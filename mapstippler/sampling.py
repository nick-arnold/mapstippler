"""Pseudo-random points inside polygons.

Rejection sampling: draw uniform points in each polygon's bounding box and keep
the ones inside. Works for any shape, including multipart polygons, holes and
detailed coastlines. Each polygon is prepared once so the inside test is fast.
"""

from __future__ import annotations

import numpy as np
import shapely

__all__ = ["random_points_in_polygons"]


def random_points_in_polygons(geoms, counts, rng=None, max_rounds: int = 1000):
    """Place counts[i] uniform random points inside geoms[i].

    Returns (xy, owner): xy is an (m, 2) float array and owner[j] is the index
    of the polygon that point j was placed in. Points are grouped by polygon.
    """
    geoms = np.asarray(geoms, dtype=object)
    counts = np.asarray(counts, dtype=np.int64)
    if geoms.shape[0] != counts.shape[0]:
        raise ValueError("need one count per polygon")
    rng = np.random.default_rng(rng)

    total = int(counts.sum())
    xy = np.empty((total, 2), dtype=np.float64)
    owner = np.repeat(np.arange(counts.size), counts)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])

    todo = np.flatnonzero(counts > 0)
    bounds = shapely.bounds(geoms[todo])
    areas = shapely.area(geoms[todo])
    for k, i in enumerate(todo):
        need = int(counts[i])
        g = geoms[i]
        if g is None or shapely.is_empty(g) or areas[k] <= 0:
            raise ValueError(f"polygon {i} has no area but was given {need} points")
        shapely.prepare(g)
        x0, y0, x1, y1 = bounds[k]
        hit_rate = areas[k] / max((x1 - x0) * (y1 - y0), 1e-300)
        got = 0
        out = xy[starts[i]:starts[i] + need]
        for _ in range(max_rounds):
            n_try = int(min(max((need - got) / hit_rate * 1.2 + 16, 64), 5_000_000))
            px = rng.uniform(x0, x1, n_try)
            py = rng.uniform(y0, y1, n_try)
            inside = shapely.contains_xy(g, px, py)
            take = min(int(inside.sum()), need - got)
            if take:
                idx = np.flatnonzero(inside)[:take]
                out[got:got + take, 0] = px[idx]
                out[got:got + take, 1] = py[idx]
                got += take
            if got == need:
                break
        else:
            raise RuntimeError(f"could not place {need} points in polygon {i}")
        shapely.destroy_prepared(g)
    return xy, owner
