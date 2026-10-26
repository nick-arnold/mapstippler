"""Dot map pipeline: pseudo-random, blue-noise and graduated dot maps.

Follows Arnold, Jenny & White (2017):

1. Pseudo-random points inside each polygon, k0 per final dot.
2. CCVT with capacity k0: each final small dot owns exactly k0 points, so every
   small dot represents exactly the unit value.
3. For each class but the largest: dots with a neighbor of the same class closer
   than that class's search radius count as merged (DBSCAN with min cluster
   size 2). The radius defaults to the dot diameter, as in the paper, but can be
   set independently. Merged dots are pooled into groups and reduced with CCVT
   again, at capacity `class_ratio`, into the next larger class.
4. Optionally, overlapping dots are nudged apart.
5. Every dot is measured: where its points came from, how far it drifted.

`make_dot_map` at the bottom runs those steps in order, and each is a `_step_*`
function that can be read, changed or tested on its own. Two small objects carry
state between the steps: `_Input` (what was read from the polygon layer) and
`_Dots` (the dots as they are built up).
"""

from __future__ import annotations

import time
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import shapely
from scipy.spatial import cKDTree

from .allocate import allocate_points, largest_remainder
from .ccvt_fast import FastCCVT
from .nudge import nudge as _nudge
from .sampling import random_points_in_polygons
from .scale import MapScale

__all__ = ["DotMapResult", "make_dot_map"]

METHODS = ("random", "blue", "graduated")
NUDGE_MODES = ("none", "top", "all")

_EQUAL_AREA_METHODS = ("albers", "equal area", "equal-area", "mollweide", "sinusoidal",
                       "eckert iv", "equal earth", "goode")


@dataclass
class DotMapResult:
    dots: "gpd.GeoDataFrame"          # one row per dot
    reconciliation: pd.DataFrame      # one row per polygon
    summary: dict
    params: dict
    timings: dict = field(default_factory=dict)


# --------------------------------------------------------------- shared state

@dataclass
class _Input:
    """What was read out of the polygon layer, as plain arrays."""
    values: np.ndarray                # total value per polygon (sum over categories)
    ids: np.ndarray                   # one id per polygon
    geoms: np.ndarray                 # shapely geometries, repaired if needed
    crs: object
    categories: tuple = ()            # category names, empty for a single-value map
    category_values: object = None    # (n_poly, n_cat) values, or None

    @property
    def n_poly(self) -> int:
        return self.geoms.size

    @property
    def n_cat(self) -> int:
        return len(self.categories)


@dataclass
class _Dots:
    """The dots as they are built up.

    `point_dot` maps each pseudo-random point to the dot that now stands for it.
    Keeping it up to date through promotion is how a large dot still knows which
    polygons its value came from.
    """
    xy: np.ndarray                    # (n, 2) positions
    cls: np.ndarray                   # (n,) class index, 0 = smallest
    value: np.ndarray                 # (n,) value each dot represents
    point_dot: np.ndarray             # (m,) dot index for each input point
    category: object = None           # (n,) category index, or None if uncategorised

    @property
    def n(self) -> int:
        return self.xy.shape[0]

    def add(self, xy, cls_index: int, value: float, category_index: int = -1) -> int:
        """Append new dots of one class. Returns the index they start at."""
        base = self.n
        n_new = xy.shape[0]
        self.xy = np.vstack([self.xy, xy])
        self.cls = np.concatenate([self.cls, np.full(n_new, cls_index, dtype=np.int64)])
        self.value = np.concatenate([self.value, np.full(n_new, float(value))])
        if self.category is not None:
            self.category = np.concatenate(
                [self.category, np.full(n_new, category_index, dtype=np.int64)])
        return base

    def reassign(self, old_dots, new_dots) -> None:
        """Point every input point that belonged to `old_dots` at `new_dots`."""
        lookup = np.full(self.n, -1, dtype=np.int64)
        lookup[old_dots] = new_dots
        moved = lookup[self.point_dot] >= 0
        self.point_dot[moved] = lookup[self.point_dot[moved]]

    def keep_only(self, keep) -> None:
        """Drop every dot not in `keep`, renumbering point_dot to match."""
        remap = np.full(self.n, -1, dtype=np.int64)
        remap[keep] = np.arange(keep.size)
        self.xy, self.cls, self.value = self.xy[keep], self.cls[keep], self.value[keep]
        if self.category is not None:
            self.category = self.category[keep]
        self.point_dot = remap[self.point_dot]

    def sizes_px(self, scale: MapScale) -> np.ndarray:
        return np.asarray(scale.sizes_px)[np.minimum(self.cls, scale.n_classes - 1)]


class _Progress:
    """Collects timings and prints step messages (or stays quiet)."""

    def __init__(self, progress):
        self.say = progress or (lambda *a, **k: None)
        self.timings: dict = {}
        self._t0 = time.perf_counter()

    @contextmanager
    def step(self, name):
        t = time.perf_counter()
        yield
        self.timings[name] = time.perf_counter() - t

    def elapsed(self, name) -> float:
        return self.timings.get(name, 0.0)

    def finish(self) -> dict:
        self.timings["total"] = time.perf_counter() - self._t0
        return self.timings


# ------------------------------------------------------------------- helpers

def _check_crs(crs):
    if crs is None:
        raise ValueError("the polygons have no coordinate system; set one first")
    if not crs.is_projected:
        raise ValueError("the polygons must be in a projected (ideally equal-area) coordinate system")
    method = ""
    try:
        method = crs.coordinate_operation.method_name.lower()
    except AttributeError:
        pass
    if not any(m in method for m in _EQUAL_AREA_METHODS):
        warnings.warn(f"'{crs.name}' may not be equal-area; dot density will be distorted "
                      "where the projection distorts area", stacklevel=3)


def _components(xy, radius):
    """Label connected groups of points linked by distances under `radius`.
    Every point gets a label; isolated points form groups of one."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n = xy.shape[0]
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    if not np.isfinite(radius):
        return np.zeros(n, dtype=np.int64)
    pairs = cKDTree(xy).query_pairs(radius * (1 - 1e-12), output_type="ndarray")
    g = coo_matrix((np.ones(pairs.shape[0]), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    return connected_components(g, directed=False)[1].astype(np.int64)


def _touching(xy, eps):
    """DBSCAN with min cluster size 2 flags every dot with a neighbor closer than eps."""
    if xy.shape[0] < 2:
        return np.zeros(xy.shape[0], dtype=bool)
    d, _ = cKDTree(xy).query(xy, k=2)
    return d[:, 1] < eps


def _pick_largest(score, labels, count):
    """Per label, a mask of the `count[label]` entries with the largest score."""
    order = np.lexsort((-score, labels))
    lab = labels[order]
    first = np.r_[0, np.flatnonzero(lab[1:] != lab[:-1]) + 1]
    rank = np.arange(order.size) - np.repeat(first, np.diff(np.r_[first, order.size]))
    mask = np.zeros(score.size, dtype=bool)
    mask[order] = rank < count[lab]
    return mask


def _nearest_same_label(xy, labels, sites, site_labels):
    """Distance from each dot to the nearest site with the same label."""
    tree = cKDTree(sites)
    k = min(8, sites.shape[0])
    d, j = tree.query(xy, k=k)
    d, j = d.reshape(xy.shape[0], k), j.reshape(xy.shape[0], k)
    same = site_labels[j] == labels[:, None]
    out = np.where(same, d, np.inf).min(axis=1)
    for i in np.flatnonzero(~np.isfinite(out)):         # rare: no same-cluster site among the k nearest
        own = sites[site_labels == labels[i]]
        out[i] = np.sqrt(((own - xy[i]) ** 2).sum(axis=1)).min()
    return out


def _majority_source(point_dot, point_poly, n_dots, n_poly):
    """For each dot: the polygon contributing most of its points, and that share."""
    codes = point_dot.astype(np.int64) * n_poly + point_poly
    uniq, cnt = np.unique(codes, return_counts=True)
    dot, poly = uniq // n_poly, uniq % n_poly
    order = np.lexsort((-cnt, dot))                 # per dot, largest count first
    dot, poly, cnt = dot[order], poly[order], cnt[order]
    first = np.ones(dot.size, dtype=bool)
    first[1:] = dot[1:] != dot[:-1]
    source = np.full(n_dots, -1, dtype=np.int64)
    top = np.zeros(n_dots, dtype=np.int64)
    source[dot[first]] = poly[first]
    top[dot[first]] = cnt[first]
    total = np.bincount(point_dot, minlength=n_dots)
    return source, top / np.maximum(total, 1)

def _overlap_counts(xy, cls, scale):
    """Number of dots of class i that overlap at least one dot of class j."""
    n = scale.n_classes
    out = np.zeros((n, n), dtype=np.int64)
    trees = [cKDTree(xy[cls == j]) if np.any(cls == j) else None for j in range(n)]
    for i in range(n):
        pts = xy[cls == i]
        if pts.shape[0] == 0:
            continue
        for j in range(n):
            if trees[j] is None:
                continue
            r = 0.5 * (scale.diameter(i) + scale.diameter(j))
            k = trees[j].query_ball_point(pts, r * (1 - 1e-9), return_length=True)
            if i == j:
                k = k - 1                           # the dot itself
            out[i, j] = int(np.count_nonzero(k > 0))
    return out


def _promote_group(xy, labels, ratio, rng, max_iterations=1000, max_rounds=5):
    """Turn each group of merged dots into larger dots of exactly `ratio` dots each.

    A group of m dots becomes m // ratio larger dots; the remaining m % ratio
    stay as they are. Which dots stay is chosen with concentration steps: start
    with the dots farthest from the group centre, build the larger dots with a
    CCVT constrained to each group, then re-pick the leftovers as the dots
    farthest from their nearest larger dot, until the choice stops changing.

    Returns (sites, owner, leftover, rounds): owner[i] is the new site for dot i,
    or -1 if dot i is a leftover.
    """
    k = labels.max() + 1
    sizes = np.bincount(labels, minlength=k)
    n_big = sizes // ratio
    n_left = sizes - n_big * ratio
    cx = np.bincount(labels, xy[:, 0], k) / sizes
    cy = np.bincount(labels, xy[:, 1], k) / sizes
    leftover = _pick_largest(np.hypot(xy[:, 0] - cx[labels], xy[:, 1] - cy[labels]), labels, n_left)
    site_group = np.repeat(np.arange(k), n_big)
    caps = np.full(site_group.size, ratio, dtype=np.int64)

    for rounds in range(1, max_rounds + 1):
        keep = np.flatnonzero(~leftover)
        opt = FastCCVT(xy[keep], caps, rng=rng, site_groups=site_group, point_groups=labels[keep])
        opt.run(max_iterations=max_iterations)
        sites = opt.site_locations()
        if not n_left.any():
            break
        dist = _nearest_same_label(xy, labels, sites, site_group)
        new_leftover = _pick_largest(dist, labels, n_left)
        if np.array_equal(new_leftover, leftover):
            break
        leftover = new_leftover

    owner = np.full(xy.shape[0], -1, dtype=np.int64)
    owner[keep] = opt.assignment()
    return sites, owner, leftover, rounds


# --------------------------------------------------------------------- steps

def _read_column(polygons, name: str, run: _Progress) -> np.ndarray:
    """One numeric column as a float array; missing values become 0."""
    col = pd.to_numeric(polygons[name], errors="coerce")
    n_missing = int(col.isna().sum())
    if n_missing:
        warnings.warn(f"{n_missing} polygons have no {name}; treated as 0", stacklevel=4)
    col = col.fillna(0).to_numpy(dtype=np.float64)
    if np.any(col < 0):
        raise ValueError(f"{name} has negative values")
    return col


def _step_read(polygons, value, id_column, scale: MapScale, run: _Progress) -> _Input:
    """Validate the layer and pull out values, ids and geometries.

    `value` is one column name, or several for a categorised map (for example
    two parties), in which case each dot is coloured by one category and the
    totals per polygon are matched exactly.
    """
    _check_crs(polygons.crs)
    if scale.radii_px is not None:
        under = [i + 1 for i in range(scale.n_classes) if scale.radius_px(i) < scale.sizes_px[i]]
        if under:
            warnings.warn(f"search radius is smaller than the dot diameter for class(es) {under}, "
                          "so dots that visibly overlap will not be promoted", stacklevel=3)

    if isinstance(value, str):
        categories, category_values = (), None
        values = _read_column(polygons, value, run)
    else:
        categories = tuple(value)
        if len(categories) < 2:
            raise ValueError("give one column name, or two or more for a categorised map")
        category_values = np.column_stack([_read_column(polygons, c, run) for c in categories])
        values = category_values.sum(axis=1)

    ids = polygons[id_column].to_numpy() if id_column else polygons.index.to_numpy()
    geoms = np.asarray(polygons.geometry.values, dtype=object)
    bad = ~shapely.is_valid(geoms)
    if bad.any():
        run.say(f"  repairing {int(bad.sum())} invalid polygons")
        geoms = geoms.copy()
        geoms[bad] = shapely.make_valid(geoms[bad])
    return _Input(values=values, ids=ids, geoms=geoms, crs=polygons.crs,
                  categories=categories, category_values=category_values)


def _step_allocate(data: _Input, unit_value: float, k0: int, method: str, run: _Progress):
    """Decide how many dots there are and how many points go in each polygon."""
    with run.step("allocate"):
        if method == "random":
            n_dots = int(round(data.values.sum() / unit_value))
            counts = largest_remainder(data.values, n_dots)
        else:
            alloc = allocate_points(data.values, unit_value, k0)
            n_dots, counts = alloc.n_dots, alloc.points
    run.say(f"[1] allocation: {data.values.sum():,.0f} total -> {n_dots:,} dots"
            + ("" if method == "random" else f", {counts.sum():,} points (k0 = {k0})"))
    return n_dots, counts


def _step_place_points(data: _Input, counts, rng, run: _Progress):
    """Scatter the pseudo-random points inside the polygons."""
    with run.step("random_points"):
        xy, point_poly = random_points_in_polygons(data.geoms, counts, rng=rng)
    run.say(f"[2] pseudo-random points: {xy.shape[0]:,} placed in "
            f"{np.count_nonzero(counts):,} polygons ({run.elapsed('random_points'):.1f} s)")
    return xy, point_poly


def _step_blue_noise(xy, n_dots: int, unit_value: float, k0: int, method: str,
                     rng, max_iterations: int, run: _Progress):
    """Reduce the points to dots with the CCVT. Returns (_Dots, number of passes)."""
    if method == "random":
        return _Dots(xy=xy.copy(), cls=np.zeros(xy.shape[0], dtype=np.int64),
                     value=np.full(xy.shape[0], float(unit_value)),
                     point_dot=np.arange(xy.shape[0])), 0

    with run.step("ccvt"):
        opt = FastCCVT(xy, np.full(n_dots, k0, dtype=np.int64), rng=rng)
        last = [time.perf_counter()]

        def report(o):
            now = time.perf_counter()
            run.say(f"      pass {o.iterations:3d}: {o.moves_last:,} points moved "
                    f"({now - last[0]:.1f} s)")
            last[0] = now

        converged = opt.run(max_iterations=max_iterations, callback=report)
    run.say(f"[3] blue noise: {n_dots:,} dots, {opt.iterations} passes, "
            f"{'converged' if converged else 'stopped at max_iterations'} "
            f"({run.elapsed('ccvt'):.1f} s)")
    return _Dots(xy=opt.site_locations(), cls=np.zeros(n_dots, dtype=np.int64),
                 value=np.full(n_dots, float(unit_value)),
                 point_dot=opt.assignment()), opt.iterations


def _step_assign_categories(dots: _Dots, point_poly, data: _Input, rng, run: _Progress):
    """Colour each dot, so every polygon's category totals come out right.

    The category has to be chosen per dot, not per point. If each point were
    given a category and each dot took the majority of its points, the larger
    category would win nearly every dot: in a 60/40 polygon with 16 points per
    dot, about 72% of dots would come out as the majority rather than 60%.

    So each polygon's dots are split by largest remainder in proportion to its
    category values, and which of its dots get which category is chosen at
    random. The dots are already a blue-noise pattern, so a random choice reads
    as evenly interleaved colours rather than clumps.
    """
    if not data.n_cat:
        return None
    with run.step("categories"):
        source, _ = _majority_source(dots.point_dot, point_poly, dots.n, data.n_poly)
        category = np.full(dots.n, -1, dtype=np.int64)
        order = np.argsort(source, kind="stable")
        starts = np.flatnonzero(np.r_[True, np.diff(source[order]) != 0])
        for start, stop in zip(starts, np.r_[starts[1:], order.size]):
            idx = order[start:stop]
            counts = largest_remainder(data.category_values[source[idx[0]]], idx.size)
            labels = np.repeat(np.arange(counts.size), counts)
            rng.shuffle(labels)
            category[idx] = labels
        dots.category = category

    shares = np.bincount(category, minlength=data.n_cat) / dots.n
    wanted = data.category_values.sum(axis=0) / data.values.sum()
    run.say("[3b] categories: " + ", ".join(
        f"{name} {got:.1%} of dots (data {want:.1%})"
        for name, got, want in zip(data.categories, shares, wanted))
        + f" ({run.elapsed('categories'):.1f} s)")
    return shares


def _promote_one_class(dots: _Dots, cls_index: int, alive, scale: MapScale, unit_value: float,
                       group_radius, rng, max_iterations: int, run: _Progress,
                       category_index: int = -1, category_name: str = ""):
    """One promotion pass: merged dots of `cls_index` become dots of the next class.

    `alive` marks the dots that have not been absorbed. Returns
    (info, stop, alive): `info` records what happened (None if there was nothing
    to consider), `stop` is True when no further class can be built, and `alive`
    is extended to cover any dots that were added.
    """
    ratio = int(round(scale.class_ratio))
    eligible_dots = alive & (dots.cls == cls_index)
    if category_index >= 0:
        # dots of different categories are promoted separately: a larger dot of
        # one colour must be worth exactly `ratio` smaller dots of that colour
        eligible_dots &= dots.category == category_index
    cand = np.flatnonzero(eligible_dots)
    if cand.size < 2:
        return None, True, alive

    eps = scale.search_radius(cls_index)
    link = (np.inf if group_radius is None
            else max(eps, group_radius * scale.diameter(cls_index + 1)))
    flagged = np.flatnonzero(_touching(dots.xy[cand], eps))
    labels = _components(dots.xy[cand[flagged]], link)
    sizes = np.bincount(labels) if flagged.size else np.zeros(0, dtype=np.int64)
    eligible = np.flatnonzero(sizes >= ratio)

    label = f" ({category_name})" if category_name else ""
    info = {"class": cls_index + 1, "category": category_name or None,
            "dots": int(cand.size), "touching": int(flagged.size),
            "link_radius_m": float(link), "clusters": int(sizes.size),
            "clusters_promoted": int(eligible.size),
            "dots_in_small_clusters": int(sizes[sizes < ratio].sum()),
            "new_dots": 0, "leftovers": 0, "rounds": 0}
    head = (f"[4] class {cls_index + 1}{label}: {info['touching']:,} of {cand.size:,} dots touch "
            f"(radius {eps:,.0f} m = {scale.radius_px(cls_index)} px); pooled within "
            + ("whole map" if not np.isfinite(link) else f"{link / scale.m_per_px:.1f} px") + ": ")
    if eligible.size == 0:
        run.say(head + f"{sizes.size:,} groups, none with {ratio} or more dots "
                       "-> nothing to promote")
        return info, True, alive

    relabel = np.full(sizes.size, -1, dtype=np.int64)
    relabel[eligible] = np.arange(eligible.size)
    keep = relabel[labels] >= 0
    members = cand[flagged[keep]]
    sites, owner, leftover, rounds = _promote_group(
        dots.xy[members], relabel[labels[keep]], ratio, rng, max_iterations)
    info.update(new_dots=int(sites.shape[0]), leftovers=int(leftover.sum()), rounds=int(rounds))
    run.say(head + f"{eligible.size:,} groups -> {sites.shape[0]:,} class {cls_index + 2} dots, "
            f"{info['leftovers']:,} leftovers, {info['dots_in_small_clusters']:,} dots in "
            f"groups under {ratio} stay")

    base = dots.add(sites, cls_index + 1, scale.value(unit_value, cls_index + 1), category_index)
    absorbed = members[owner >= 0]
    dots.reassign(absorbed, base + owner[owner >= 0])
    alive = np.concatenate([alive, np.ones(sites.shape[0], dtype=bool)])
    alive[absorbed] = False
    return info, False, alive


def _step_promote(dots: _Dots, scale: MapScale, unit_value: float, group_radius,
                  method: str, rng, max_iterations: int, run: _Progress,
                  data_categories=()) -> list:
    """Build the larger classes, one at a time, then drop the dots they absorbed."""
    if method != "graduated":
        return []
    promotions = []
    groups = ([(-1, "")] if not data_categories
              else list(enumerate(data_categories)))
    with run.step("graduated"):
        alive = np.ones(dots.n, dtype=bool)
        for c in range(scale.n_classes - 1):
            stopped = 0
            for cat_index, cat_name in groups:
                info, stop, alive = _promote_one_class(
                    dots, c, alive, scale, unit_value, group_radius, rng,
                    max_iterations, run, cat_index, cat_name)
                if info is not None:
                    promotions.append(info)
                stopped += bool(stop)
            if stopped == len(groups):
                break
        dots.keep_only(np.flatnonzero(alive))
    return promotions


def _step_nudge(dots: _Dots, scale: MapScale, nudge: str, nudge_max: float,
                rng, run: _Progress):
    """Push overlapping dots apart. Returns (distance moved per dot, summary)."""
    moved = np.zeros(dots.n)
    if nudge == "none" or dots.n < 2:
        return moved, None

    with run.step("nudge"):
        radii = 0.5 * dots.sizes_px(scale) * scale.m_per_px
        idx = (np.flatnonzero(dots.cls == dots.cls.max()) if nudge == "top"
               else np.arange(dots.n))
        res = _nudge(dots.xy[idx], radii[idx], weight=radii[idx] ** 2,
                     max_disp=nudge_max * 2 * radii[idx], rng=rng)
        dots.xy = dots.xy.copy()
        dots.xy[idx] = res.xy
        moved[idx] = res.moved

    mv = moved[moved > 0] / scale.m_per_px
    info = {"mode": nudge, "dots_considered": int(idx.size),
            "overlapping_before": res.overlapping_before,
            "overlapping_after": res.overlapping_after,
            "dots_moved": int(mv.size), "at_cap": res.capped, "rounds": res.iterations,
            "move_px_median": float(np.median(mv)) if mv.size else 0.0,
            "move_px_max": float(mv.max()) if mv.size else 0.0}
    run.say(f"[4b] nudge ({'largest-class' if nudge == 'top' else 'all'} dots, cap "
            f"{nudge_max:g} diameter{'s' if nudge_max != 1 else ''}): "
            f"{res.overlapping_before:,} overlapping -> {res.overlapping_after:,}; "
            f"{mv.size:,} moved, median {info['move_px_median']:.2f} px, max "
            f"{info['move_px_max']:.2f} px, {res.capped:,} at the cap "
            f"({run.elapsed('nudge'):.1f} s)")
    return moved, info


@dataclass
class _Measurements:
    source: np.ndarray        # polygon that contributed most of each dot's points
    purity: np.ndarray        # share of its points that came from that polygon
    spread: np.ndarray        # RMS distance from a dot to the points it stands for
    container: np.ndarray     # polygon each dot actually falls in, -1 if none
    inside: np.ndarray        # container == source
    displacement: np.ndarray  # distance beyond the source polygon's edge
    overlaps: object          # (classes x classes) counts, or None


def _step_measure(dots: _Dots, points_xy, point_poly, data: _Input, scale: MapScale,
                  method: str, run: _Progress) -> _Measurements:
    """Where each dot's value came from, how far it drifted, what it overlaps."""
    with run.step("diagnostics"):
        source, purity = _majority_source(dots.point_dot, point_poly, dots.n, data.n_poly)

        off = points_xy - dots.xy[dots.point_dot]
        n_pts = np.bincount(dots.point_dot, minlength=dots.n)
        spread = np.sqrt(np.bincount(dots.point_dot, (off * off).sum(axis=1), dots.n)
                         / np.maximum(n_pts, 1))

        pts = shapely.points(dots.xy)
        pi, gi = shapely.STRtree(data.geoms).query(pts, predicate="within")
        container = np.full(dots.n, -1, dtype=np.int64)
        container[pi[::-1]] = gi[::-1]                      # first match wins
        inside = container == source

        displacement = np.zeros(dots.n)
        outside = np.flatnonzero(~inside)
        if outside.size:
            displacement[outside] = shapely.distance(data.geoms[source[outside]], pts[outside])

        overlaps = _overlap_counts(dots.xy, dots.cls, scale) if method != "random" else None
    return _Measurements(source, purity, spread, container, inside, displacement, overlaps)


# ------------------------------------------------------------------- outputs

def _build_dots_frame(dots: _Dots, m: _Measurements, data: _Input, scale: MapScale, nudge_moved):
    import geopandas as gpd
    extra = {}
    if data.n_cat:
        # ordered categorical, so colours and legends follow the order the
        # columns were given rather than whichever dot happens to come first
        extra["category"] = pd.Categorical(np.asarray(data.categories)[dots.category],
                                           categories=list(data.categories), ordered=True)
    return gpd.GeoDataFrame({
        "dot_id": np.arange(dots.n),
        "class": dots.cls + 1,
        "value": dots.value,
        "size_px": dots.sizes_px(scale),
        "source_id": data.ids[m.source],
        "purity": np.round(m.purity, 4),
        "inside_source": m.inside,
        "displacement_px": np.round(m.displacement / scale.m_per_px, 3),
        "spread_px": np.round(m.spread / scale.m_per_px, 3),
        "nudge_px": np.round(nudge_moved / scale.m_per_px, 3),
        **extra,
    }, geometry=shapely.points(dots.xy), crs=data.crs)


def _build_reconciliation(dots: _Dots, m: _Measurements, data: _Input, counts, method: str):
    """Per polygon: its value against the value of the dots that fall inside it."""
    placed = m.container >= 0
    value_inside = np.bincount(m.container[placed], weights=dots.value[placed],
                               minlength=data.n_poly)
    rec = pd.DataFrame({"id": data.ids, "value": data.values,
                        "points": counts if method != "random" else 0,
                        "value_of_dots_inside": value_inside,
                        "difference": value_inside - data.values})
    for c in range(int(dots.cls.max()) + 1):
        sel = placed & (dots.cls == c)
        rec[f"class{c + 1}_dots_inside"] = np.bincount(m.container[sel], minlength=data.n_poly)
    rec["dots_sourced"] = np.bincount(m.source, minlength=data.n_poly)
    for k, name in enumerate(data.categories):
        rec[f"{name}_value"] = data.category_values[:, k]
        sel = dots.category == k
        rec[f"{name}_value_of_dots"] = np.bincount(m.source[sel], weights=dots.value[sel],
                                                   minlength=data.n_poly)
    return rec


def _build_summary(dots: _Dots, m: _Measurements, data: _Input, rec, scale: MapScale,
                   unit_value: float, ccvt_passes: int, promotions: list, nudge_info):
    total = dots.value.sum()
    outside = np.flatnonzero(~m.inside)
    px = scale.m_per_px
    return {
        "value_total": float(data.values.sum()),
        "value_of_all_dots": float(total),
        "dots_by_class": {int(c + 1): int(np.count_nonzero(dots.cls == c))
                          for c in np.unique(dots.cls)},
        "ccvt_passes": ccvt_passes,
        "categories": ({name: {"value": float(data.category_values[:, k].sum()),
                               "value_of_dots": float(dots.value[dots.category == k].sum()),
                               "dots": int(np.count_nonzero(dots.category == k))}
                        for k, name in enumerate(data.categories)}
                       if data.n_cat else None),
        "promotions": promotions,
        "nudge": nudge_info,
        "drift": {
            "dots_outside_source": int(outside.size),
            "value_share_outside_source": float(dots.value[~m.inside].sum() / total),
            "dots_in_no_polygon": int(np.count_nonzero(m.container < 0)),
            "displacement_px_median_of_outside":
                float(np.median(m.displacement[outside]) / px) if outside.size else 0.0,
            "displacement_px_p95_of_outside":
                float(np.percentile(m.displacement[outside], 95) / px) if outside.size else 0.0,
            "mean_purity_value_weighted": float(np.sum(m.purity * dots.value) / total),
        },
        "spread_px_max_by_class": {int(c + 1): float(m.spread[dots.cls == c].max() / px)
                                   for c in np.unique(dots.cls)},
        "overlapping_dots": ({f"class{i + 1}_with_class{j + 1}": int(m.overlaps[i, j])
                              for i in range(m.overlaps.shape[0])
                              for j in range(m.overlaps.shape[1])}
                             if m.overlaps is not None else None),
        "polygons_abs_difference_over_half_unit":
            int(np.count_nonzero(np.abs(rec["difference"]) > unit_value / 2)),
    }


# --------------------------------------------------------------------- entry

def make_dot_map(polygons, value, unit_value: float, scale: MapScale, k0: int = 16,
                 id_column: str | None = None, method: str = "graduated", rng=0,
                 max_iterations: int = 1000, group_radius: float | None = 3.0,
                 nudge: str = "none", nudge_max: float = 1.0,
                 progress=print) -> DotMapResult:
    """Make a dot map from a polygon GeoDataFrame.

    polygons must be in a projected, ideally equal-area, coordinate system.

    value: the column of counts to map, or a list of columns for a categorised
        map (for example ["dem", "rep"]). Each dot then carries one category and
        is coloured by it; every polygon's split is matched by largest remainder,
        and each category is promoted into larger dots separately, so a large dot
        is worth exactly `class_ratio` smaller dots of its own category.

    method: "random"     pseudo-random dots at the final density (no CCVT)
            "blue"       blue-noise dots of a single class
            "graduated"  blue-noise dots promoted into larger classes

    unit_value: what one dot of the smallest class represents.
    k0: pseudo-random points per dot. Higher gives a better blue-noise pattern
        and changes nothing about what the map means.

    group_radius: how merged dots are pooled before being promoted. They pool
        with other merged dots closer than this many diameters of the larger
        class, so scattered merges in one dense area pool together while
        separate towns never do. 0 pools only dots that touch directly; None
        pools every merged dot on the map, as in Arnold et al. (2017), which can
        place a dot between distant towns.

    nudge: "none", "top" (push apart overlapping dots of the largest class) or
        "all" (every class; larger dots are harder to move, in proportion to
        their area). Each dot moves at most `nudge_max` of its own diameters.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    if nudge not in NUDGE_MODES:
        raise ValueError(f"nudge must be one of {NUDGE_MODES}")
    rng = np.random.default_rng(rng)
    run = _Progress(progress)

    data = _step_read(polygons, value, id_column, scale, run)
    params = {"value": value, "categories": list(data.categories) or None, "unit_value": unit_value, "k0": k0, "method": method,
              "group_radius": group_radius, "nudge": nudge, "nudge_max": nudge_max,
              "crs": data.crs.to_string(), "polygons": data.n_poly, **scale.describe()}

    n_dots, counts = _step_allocate(data, unit_value, k0, method, run)
    points_xy, point_poly = _step_place_points(data, counts, rng, run)
    dots, ccvt_passes = _step_blue_noise(points_xy, n_dots, unit_value, k0, method,
                                         rng, max_iterations, run)
    _step_assign_categories(dots, point_poly, data, rng, run)
    promotions = _step_promote(dots, scale, unit_value, group_radius, method,
                               rng, max_iterations, run, data.categories)
    nudge_moved, nudge_info = _step_nudge(dots, scale, nudge, nudge_max, rng, run)
    measured = _step_measure(dots, points_xy, point_poly, data, scale, method, run)

    frame = _build_dots_frame(dots, measured, data, scale, nudge_moved)
    rec = _build_reconciliation(dots, measured, data, counts, method)
    summary = _build_summary(dots, measured, data, rec, scale, unit_value,
                             ccvt_passes, promotions, nudge_info)

    timings = run.finish()
    run.say(f"[5] done: {dots.n:,} dots, "
            f"{summary['drift']['value_share_outside_source']:.1%} of dot value "
            f"outside its source polygon ({timings['total']:.1f} s total)")
    return DotMapResult(dots=frame, reconciliation=rec, summary=summary,
                        params=params, timings=timings)
