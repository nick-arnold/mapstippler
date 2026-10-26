"""Unit value advisor.

Previews, before any dots are placed, what a candidate unit value would do at a
given scale and set of dot sizes: how many dots, and how much of the value sits
where each dot class would coalesce.

How coalescence is judged
-------------------------
Dots are spread uniformly inside their polygons, so the map's density is the
polygon value spread over the polygon area. What the eye judges is density over
a neighborhood a few dots wide, not per polygon: at small scales many polygons
(urban tracts, small counties) are smaller than a single dot, and their dots
necessarily mix with their neighbors'.

So density is measured on a grid whose cells are `window` dot diameters wide
(3 by default), separately for each class. It is estimated by placing a large
sample of points in the polygons in proportion to their values and counting them
per cell, averaged over two grid offsets to reduce alignment effects.

For a cell of area A holding value V, a class with dot value u and search radius
d (the distance at which two of its dots count as merged, by default the dot
diameter) would show n = V / u dots at mean spacing s = sqrt(A / n). In a CCVT blue-noise
pattern at k0 = 16, nearest-neighbor distances are about 0.5 s at the smallest
and 0.82 s on average (measured with this library), so:

    d / s < 0.5          dots stay clear of each other
    0.5 <= d / s < 0.8   some dots touch
    d / s >= 0.8         a typical dot touches its neighbor: the class coalesces

Each class is judged as if all of the value there were shown in that class,
which is the question DBSCAN promotion asks: is this class too crowded here?
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .allocate import allocate_points, largest_remainder
from .sampling import random_points_in_polygons
from .scale import MapScale

__all__ = ["TOUCH", "COALESCE", "density_sample", "coalescing_share",
           "candidate_unit_values", "advise", "ink_report"]

TOUCH = 0.5       # d / s at which some dots start to touch
COALESCE = 0.8    # d / s at which a typical dot touches its neighbor


def density_sample(geoms, values, n_sample: int = 1_000_000, rng=0):
    """Sample points in proportion to value. Returns (xy, value_per_point)."""
    values = np.asarray(values, dtype=np.float64)
    n_sample = int(min(n_sample, max(values.sum(), 1)))
    counts = largest_remainder(values, n_sample)
    xy, _ = random_points_in_polygons(np.asarray(geoms, dtype=object), counts, rng=rng)
    return xy, values.sum() / max(n_sample, 1)


def _cell_values(xy, w, cell, offset):
    ij = np.floor((xy - xy.min(axis=0) + offset) / cell).astype(np.int64)
    key = ij[:, 0] * (ij[:, 1].max() + 1) + ij[:, 1]
    _, cnt = np.unique(key, return_counts=True)
    return cnt * w


def coalescing_share(xy, w, unit: float, diameter: float, window: float = 3.0,
                     threshold: float = COALESCE) -> float:
    """Share of total value in neighborhoods where dots of this value and size
    reach `threshold` (d / s)."""
    cell = window * diameter
    area = cell * cell
    total = xy.shape[0] * w
    shares = []
    for off in (0.0, 0.5 * cell):
        v = _cell_values(xy, w, cell, off)
        ratio = diameter * np.sqrt((v / unit) / area)
        shares.append(v[ratio >= threshold].sum() / total)
    return float(np.mean(shares))


def ink_report(scale: MapScale) -> list:
    """Per class: what share of ink promotion keeps, and how much more density it absorbs.

    Promoting k dots of diameter d0 into one of diameter d1 keeps
    (d1 / d0)^2 / k of the ink. Much below 1 and the densest places lighten as
    they promote, leaving holes; at 1 the class absorbs no extra density at all,
    which is the area-proportional case the paper argues against.
    """
    out = []
    for c in range(1, scale.n_classes):
        ink = scale.ink_kept(c)
        out.append({"class": c + 1, "from_px": scale.sizes_px[c - 1], "to_px": scale.sizes_px[c],
                    "diameter_ratio": scale.sizes_px[c] / scale.sizes_px[c - 1],
                    "ink_kept": ink, "density_absorbed": 1.0 / ink if ink else float("inf")})
    return out


def candidate_unit_values(total: float, min_dots: int = 1_000, max_dots: int = 2_000_000) -> list:
    """Round unit values (1, 2, 5 x 10^k) giving between min_dots and max_dots dots."""
    lo, hi = total / max_dots, total / min_dots
    out = []
    for exp in range(int(np.floor(np.log10(max(lo, 1e-12)))) - 1, int(np.ceil(np.log10(hi))) + 1):
        for m in (1, 2, 5):
            u = m * 10.0 ** exp
            if lo <= u <= hi and u >= 1:
                out.append(int(u) if float(u).is_integer() else u)
    return out


def advise(geoms, values, scale: MapScale, unit_values=None, k0: int = 16,
           max_top_class_coalesce: float = 0.05, window: float = 3.0,
           n_sample: int = 1_000_000, rng=0) -> pd.DataFrame:
    """One row per candidate unit value.

    Columns
    -------
    unit            value of one smallest-class dot
    dots            number of smallest-class dots if nothing were promoted
    cN_touch        share of total value where class N starts touching
    cN_coal         share of total value where class N coalesces
    max_round       largest per-polygon rounding error, in value units
    suggested       True for the smallest unit value whose largest class
                    coalesces in no more than `max_top_class_coalesce` of the
                    total value: the most detail the classes can still resolve.
    """
    values = np.asarray(values, dtype=np.float64)
    total = values.sum()
    if unit_values is None:
        unit_values = candidate_unit_values(total)
    xy, w = density_sample(geoms, values, n_sample=n_sample, rng=rng)

    rows = []
    for u in unit_values:
        row = {"unit": u, "dots": int(round(total / u))}
        for c in range(scale.n_classes):
            uc, d = scale.value(u, c), scale.search_radius(c)
            row[f"c{c + 1}_touch"] = coalescing_share(xy, w, uc, d, window, TOUCH)
            row[f"c{c + 1}_coal"] = coalescing_share(xy, w, uc, d, window, COALESCE)
        row["max_round"] = float(np.abs(allocate_points(values, u, k0).error).max())
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("unit").reset_index(drop=True)

    top = f"c{scale.n_classes}_coal"
    ok = df.index[df[top] <= max_top_class_coalesce]
    df["suggested"] = False
    if len(ok):
        df.loc[ok[0], "suggested"] = True
    return df
