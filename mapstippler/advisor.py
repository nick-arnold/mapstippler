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
           "candidate_unit_values", "advise", "ink_report", "touching_share",
           "estimate_dots_by_class"]

TOUCH = 0.5       # d / s at which some dots start to touch
COALESCE = 0.8    # d / s at which a typical dot touches its neighbor

# What share of dots in a CCVT blue-noise pattern have a neighbour closer than
# d, as a function of d / s, where s is the mean spacing sqrt(area / n).
# Measured from this library's own output: 4,000 dots at k0 = 16, interior only.
# The transition is sharp, which is the point of blue noise: dots keep their
# distance until the pattern is genuinely too crowded to.
_NN_RATIO = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                      0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                      1.05, 1.1, 1.15, 1.2, 1.25, 1.3])
_NN_SHARE = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0007, 0.0007, 0.0014,
                      0.0064, 0.0233, 0.0676, 0.1496, 0.2851, 0.4719, 0.6671, 0.8207,
                      0.9172, 0.9699, 0.9894, 0.9972, 0.9993, 1.0, 1.0, 1.0])


def touching_share(ratio):
    """Share of blue-noise dots with a neighbour closer than `ratio` x mean spacing."""
    return np.interp(ratio, _NN_RATIO, _NN_SHARE, left=0.0, right=1.0)


def estimate_dots_by_class(xy, point_value: float, unit: float, scale: MapScale,
                           window: float = 3.0):
    """Estimate how many dots of each class a unit value would produce.

    Walks the promotion chain the way the pipeline does, but over densities
    rather than placed dots. `xy` is a sample of points drawn in proportion to
    value (see `density_sample`), each standing for `point_value`.

    For each class in turn, the sample is gridded into neighbourhoods a few dots
    wide. In a neighbourhood of area A holding value V, a class with dot value u
    places n = V / u dots at mean spacing s = sqrt(A / n), and the share of them
    with a neighbour inside the class's search radius is what gets promoted. That
    share of the value is carried up to the next class, and the rest stays. Each
    class is gridded at its own size, since a larger dot's neighbourhood is
    larger.

    An estimate, not a simulation: it assumes dots spread evenly within a
    neighbourhood and ignores how merged dots are pooled into groups before
    promotion, so counts are approximate. It is for choosing a unit value, not
    for reporting.
    """
    xy = np.asarray(xy, dtype=np.float64)
    weight = np.full(xy.shape[0], float(point_value))   # value each point still carries
    origin = xy.min(axis=0)
    counts = []
    for c in range(scale.n_classes):
        u = scale.value(unit, c)
        if c == scale.n_classes - 1:
            counts.append(weight.sum() / u)             # nothing above to promote into
            break
        cell = window * scale.diameter(c)
        ij = np.floor((xy - origin) / cell).astype(np.int64)
        key = ij[:, 0] * (ij[:, 1].max() + 2) + ij[:, 1]
        cell_id, inverse = np.unique(key, return_inverse=True)
        value_in_cell = np.bincount(inverse, weight, cell_id.size)

        n = value_in_cell / u
        spacing = np.sqrt(np.divide(cell * cell, n, out=np.full_like(n, np.inf), where=n > 0))
        share = touching_share(scale.search_radius(c) / spacing)
        counts.append(float((n * (1.0 - share)).sum()))
        weight = weight * share[inverse]                # carry the promoted value up
    return [int(round(x)) for x in counts]


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
    estN            estimated dots of class N that the map would actually show
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
        for c, n in enumerate(estimate_dots_by_class(xy, w, u, scale)):
            row[f"est{c + 1}"] = n
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
