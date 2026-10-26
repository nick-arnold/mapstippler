"""Capacity-constrained Voronoi tessellation (CCVT).

A faithful Python port of the discrete CCVT optimizer by Balzer, Schlömer and
Deussen (2009), "Capacity-constrained point distributions: a variant of Lloyd's
method", following the logic of Michael Balzer's reference C++ implementation
(github.com/michaelbalzer/ccvt, GPL-3.0).

This is the baseline version. It deliberately keeps the original algorithm's
behavior, including its all-pairs loop and its "stop only when a full pass makes
no swaps" rule, so it can be validated against the C++ code before the faster
methods (Li et al. 2010; Ahmed & Deussen 2017) are added.

Terminology follows the paper:
    points  the input samples (here: pseudo-random dots)
    sites   the output dots; each site owns exactly `capacity` points

Additions relative to the original:
    * every point carries its original index, so the result reports which input
      points each site absorbed (needed later for drift and source tracking);
    * `true_energy()` computes the actual CCVT energy. The original only keeps a
      running tally updated during swaps, which goes stale when sites move to
      their centroids. `reported_energy` reproduces that tally for validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

__all__ = [
    "CCVTOptimizer",
    "CCVTResult",
    "split_capacity",
    "ccvt",
    "reduce_points",
]


def split_capacity(n_points: int, n_sites: int) -> np.ndarray:
    """Divide `n_points` among `n_sites` so capacities differ by at most one.

    Mirrors the reference example program: site i receives
    remaining // (n_sites - i), so any remainder goes to the last sites.
    """
    if n_sites <= 0:
        raise ValueError("n_sites must be positive")
    if n_points < n_sites:
        raise ValueError("need at least one point per site")
    capacities = np.empty(n_sites, dtype=np.int64)
    remaining = n_points
    for i in range(n_sites):
        c = remaining // (n_sites - i)
        capacities[i] = c
        remaining -= c
    return capacities


def _sq_dist(points: np.ndarray, center: np.ndarray) -> np.ndarray:
    """Squared distances, computed in the same operation order as the C++ metric."""
    dx = points[:, 0] - center[0]
    dy = points[:, 1] - center[1]
    return dx * dx + dy * dy


def _centroid(points: np.ndarray) -> np.ndarray:
    """Centroid with sequential summation, matching the C++ accumulation order."""
    n = points.shape[0]
    return np.array(
        [np.cumsum(points[:, 0])[-1] / n, np.cumsum(points[:, 1])[-1] / n]
    )


def nearest_available_assignment(locations, capacities, points) -> list:
    """Assign each point to its nearest site that still has free capacity.

    Same greedy rule as the reference kd-tree initialization, including the
    order: the original consumes the point list from the back. Returns, for
    each site, the array of point indices it owns (in assignment order).
    """
    n = locations.shape[0]
    remaining = np.asarray(capacities, dtype=np.int64).copy()
    owned_idx = [[] for _ in range(n)]

    active = np.arange(n)                       # sites in the current search tree
    tree = cKDTree(locations)
    full_in_tree = 0

    for p_index in range(points.shape[0] - 1, -1, -1):
        p = points[p_index]
        k = min(8, active.size)
        chosen = -1
        while True:
            _, nn = tree.query(p, k=k)
            nn = np.atleast_1d(nn)
            for t in nn:
                s = active[t]
                if remaining[s] > 0:
                    chosen = s
                    break
            if chosen >= 0 or k >= active.size:
                break
            k = min(k * 4, active.size)
        owned_idx[chosen].append(p_index)
        remaining[chosen] -= 1
        if remaining[chosen] == 0:
            full_in_tree += 1
            # Rebuild the tree on available sites once half of it is full.
            if full_in_tree * 2 > active.size and remaining.any():
                active = np.flatnonzero(remaining > 0)
                tree = cKDTree(locations[active])
                full_in_tree = 0

    return [np.array(ix, dtype=np.int64) for ix in owned_idx]


class CCVTOptimizer:
    """Discrete CCVT optimizer (Balzer et al. 2009), baseline port.

    Parameters
    ----------
    sites : (n, 2) array
        Initial site locations.
    capacities : (n,) int array
        Number of points each site must own. Must sum to the number of points.
        Sites with zero capacity are ignored, as in the original.
    points : (m, 2) array
        The input points.
    """

    def __init__(self, sites, capacities, points):
        sites = np.asarray(sites, dtype=np.float64)
        capacities = np.asarray(capacities, dtype=np.int64)
        points = np.asarray(points, dtype=np.float64)
        if sites.ndim != 2 or sites.shape[1] != 2:
            raise ValueError("sites must be an (n, 2) array")
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("points must be an (m, 2) array")
        if capacities.shape != (sites.shape[0],):
            raise ValueError("need one capacity per site")
        if capacities.sum() != points.shape[0]:
            raise ValueError("the sum of the site capacities must equal the number of points")

        keep = capacities > 0
        self.site_ids = np.flatnonzero(keep)          # original index of each active site
        self.locations = sites[keep].copy()
        self.capacities = capacities[keep].copy()
        self.n = self.locations.shape[0]
        self._points_in = points

        self._assign_initial()

        # Per-site state, as in the reference Entry struct.
        self.sq_radius = np.array([_sq_dist(p, c).max() for p, c in zip(self.pts, self.locations)])
        self.radius = np.sqrt(self.sq_radius)
        self.stable = np.zeros(self.n, dtype=bool)
        self.entry_energy = np.array([_sq_dist(p, c).sum() for p, c in zip(self.pts, self.locations)])
        self.iterations = 0

    # ------------------------------------------------------------------ setup

    def _assign_initial(self) -> None:
        self.idx = nearest_available_assignment(self.locations, self.capacities, self._points_in)
        self.pts = [self._points_in[ix].copy() for ix in self.idx]

    # ----------------------------------------------------------- optimization

    def _row_candidates(self, i: int, start: int) -> np.ndarray:
        """Sites j >= start that the reference loop would not skip when paired with i."""
        js = np.arange(start, self.n)
        if js.size == 0:
            return js
        d = self.locations[js] - self.locations[i]
        dist = np.sqrt(d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1])
        skip = (self.stable[i] & self.stable[js]) | (dist > self.radius[i] + self.radius[js])
        return js[~skip]

    def _update_entry(self, e: int) -> None:
        pts = self.pts[e]
        self.sq_radius[e] = _sq_dist(pts, self.locations[e]).max()
        self.radius[e] = np.sqrt(self.sq_radius[e])
        self.stable[e] = False

    def _process_pair(self, i: int, j: int, centroidal: bool) -> bool:
        """Exchange points between sites i and j. Returns True if any swap happened."""
        a, b = i, j
        if self.pts[a].shape[0] > self.pts[b].shape[0]:
            a, b = b, a
        loc_a, loc_b = self.locations[a], self.locations[b]
        pts_a, pts_b = self.pts[a], self.pts[b]
        max_sq = max(self.sq_radius[a], self.sq_radius[b])

        # Candidates from a: points within the bounding radius of b.
        d_a_to_b = _sq_dist(pts_a, loc_b)
        cand_a = np.flatnonzero(d_a_to_b <= max_sq)
        if cand_a.size == 0:
            return False
        self_a = _sq_dist(pts_a[cand_a], loc_a)
        other_a = d_a_to_b[cand_a]
        gain_a = self_a - other_a
        min_energy = -gain_a.max()

        # Candidates from b, pre-filtered as in the reference.
        d_b_to_a = _sq_dist(pts_b, loc_a)
        near_b = np.flatnonzero(d_b_to_a <= max_sq)
        if near_b.size == 0:
            return False
        self_b = _sq_dist(pts_b[near_b], loc_b)
        other_b = d_b_to_a[near_b]
        gain_b = self_b - other_b
        keep = gain_b > min_energy
        if not keep.any():
            return False
        cand_b, self_b, other_b, gain_b = near_b[keep], self_b[keep], other_b[keep], gain_b[keep]

        # The reference pops both heaps in descending order of gain and swaps
        # the s-th best of each while their combined gain is positive.
        order_a = np.argsort(-gain_a, kind="stable")
        order_b = np.argsort(-gain_b, kind="stable")
        max_swaps = min(order_a.size, order_b.size)
        ok = (gain_a[order_a[:max_swaps]] + gain_b[order_b[:max_swaps]]) > 0
        swaps = max_swaps if ok.all() else int(np.argmin(ok))
        if swaps == 0:
            return False

        sa, sb = order_a[:swaps], order_b[:swaps]
        rows_a, rows_b = cand_a[sa], cand_b[sb]

        # Running energy tally, accumulated in the reference's order.
        for s in range(swaps):
            self.entry_energy[a] += other_b[sb[s]] - self_a[sa[s]]
            self.entry_energy[b] += other_a[sa[s]] - self_b[sb[s]]

        tmp_p = pts_a[rows_a].copy()
        pts_a[rows_a] = pts_b[rows_b]
        pts_b[rows_b] = tmp_p
        idx_a, idx_b = self.idx[a], self.idx[b]
        tmp_i = idx_a[rows_a].copy()
        idx_a[rows_a] = idx_b[rows_b]
        idx_b[rows_b] = tmp_i

        if centroidal:
            self.locations[a] = _centroid(pts_a)
            self.locations[b] = _centroid(pts_b)
        self._update_entry(a)
        self._update_entry(b)
        return True

    def optimize(self, centroidal: bool = True) -> bool:
        """Run one full pass over all site pairs. Returns True if nothing changed."""
        stability = np.ones(self.n, dtype=bool)
        for i in range(self.n - 1):
            js = self._row_candidates(i, i + 1)
            pos = 0
            while pos < js.size:
                j = int(js[pos])
                pos += 1
                if self._process_pair(i, j, centroidal):
                    stability[i] = False
                    stability[j] = False
                    # Site i moved, so re-check the rest of this row against its new state.
                    js = self._row_candidates(i, j + 1)
                    pos = 0
        self.stable[:] = stability
        self.iterations += 1
        return bool(stability.all())

    def run(self, max_iterations: int | None = None, centroidal: bool = True, callback=None) -> bool:
        """Optimize until stable (or until `max_iterations`). Returns True if converged."""
        while max_iterations is None or self.iterations < max_iterations:
            stable = self.optimize(centroidal)
            if callback is not None:
                callback(self)
            if stable:
                return True
        return False

    # ---------------------------------------------------------------- results

    @property
    def reported_energy(self) -> float:
        """The reference implementation's running energy tally (for validation)."""
        return float(self.entry_energy.sum())

    def true_energy(self) -> float:
        """Sum of squared distances from every point to the site that owns it."""
        return float(sum(_sq_dist(p, c).sum() for p, c in zip(self.pts, self.locations)))

    def assignment(self) -> np.ndarray:
        """For each input point, the index (into the original `sites`) of its owner."""
        owner = np.empty(self._points_in.shape[0], dtype=np.int64)
        for e, ix in enumerate(self.idx):
            owner[ix] = self.site_ids[e]
        return owner


@dataclass
class CCVTResult:
    sites: np.ndarray            # (n, 2) final site locations
    capacities: np.ndarray       # (n,) points owned by each site
    assignment: np.ndarray       # (m,) owning site for each input point
    iterations: int
    converged: bool
    energy: float                # true CCVT energy of the result
    history: list = field(default_factory=list)   # true energy after each pass


def ccvt(points, n_sites: int, initial_sites=None, rng=None,
         max_iterations: int | None = None, record_history: bool = False,
         method: str = "fast", tol: float | None = None) -> CCVTResult:
    """Reduce `points` to `n_sites` blue-noise sites with the CCVT method.

    method="fast" (default) uses the accelerated optimizer in `ccvt_fast.py`,
    starting from the recursive median-cut initialization unless
    `initial_sites` is given.

    method="baseline" uses the faithful port of the original algorithm. By
    default its initial sites are a random subset of the input points, which
    already follows the input density.
    """
    if method == "fast":
        from .ccvt_fast import ccvt_fast
        return ccvt_fast(points, n_sites, initial_sites=initial_sites, rng=rng,
                         max_iterations=1000 if max_iterations is None else max_iterations,
                         tol=tol, record_history=record_history)
    if method != "baseline":
        raise ValueError("method must be 'fast' or 'baseline'")

    points = np.asarray(points, dtype=np.float64)
    rng = np.random.default_rng(rng)
    if initial_sites is None:
        initial_sites = points[rng.choice(points.shape[0], size=n_sites, replace=False)]
    capacities = split_capacity(points.shape[0], n_sites)
    opt = CCVTOptimizer(initial_sites, capacities, points)
    history = []
    cb = (lambda o: history.append(o.true_energy())) if record_history else None
    converged = opt.run(max_iterations=max_iterations, callback=cb)

    sites = np.asarray(initial_sites, dtype=np.float64).copy()
    sites[opt.site_ids] = opt.locations
    return CCVTResult(
        sites=sites,
        capacities=capacities,
        assignment=opt.assignment(),
        iterations=opt.iterations,
        converged=converged,
        energy=opt.true_energy(),
        history=history,
    )


def reduce_points(points, factor: float, **kwargs) -> CCVTResult:
    """Reduce the number of points by `factor` (the paper's reduction factor k).

    The number of sites is round(m / factor). When m does not divide evenly,
    capacities differ by at most one, as in the reference implementation.
    """
    points = np.asarray(points, dtype=np.float64)
    n_sites = max(1, int(round(points.shape[0] / factor)))
    return ccvt(points, n_sites, **kwargs)
