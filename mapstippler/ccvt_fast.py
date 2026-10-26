"""Fast capacity-constrained Voronoi tessellation (CCVT).

Same optimization problem as the baseline port in `ccvt.py` (Balzer et al. 2009),
solved with the accelerations from:

    Li, Nehab, Wei, Sander & Fu (2010). Fast capacity constrained Voronoi
    tessellation. I3D '10.
    Ahmed & Deussen (2017). Capacity-constrained Voronoi tessellation revisited.
    Computer Graphics & Visual Computing (CGVC).

What changes relative to the baseline:

1. Neighbor-only exchanges. Sites only trade points with their Delaunay
   neighbors, instead of every pair of sites being considered.
2. Linear keys and capacity cuts. For a pair of sites (i, j), the energy
   difference for a point p is |p - sj|^2 - |p - si|^2 = 2 (si - sj) . (p - mid),
   a linear function of p (its position along the line between the sites).
   Pooling the pair's points and giving site i the c_i points with the largest
   keys is the optimal exchange between the two sites.
3. Recursive median-cut initialization. The point set is recursively sliced
   perpendicular to its major axis into pieces of exactly the right sizes; each
   site starts at the centroid of its slice.
4. Batched pairs. Pairs that share no site are independent, so a whole batch of
   them (a matching in the neighbor graph) is processed at once with numpy. This
   gives the same result as processing those pairs one after another.
5. Settled pairs are skipped. A pair is only re-examined if one of its sites
   changed in the previous pass, or if the two sites have just become neighbors.

Storage follows Ahmed & Deussen: capacities never change, so site s always owns
the fixed block of slots offsets[s] : offsets[s] + capacities[s]. Exchanging
points only permutes which point index sits in which slot.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import Delaunay, cKDTree

try:  # location of QhullError differs across scipy versions
    from scipy.spatial import QhullError
except ImportError:  # pragma: no cover
    from scipy.spatial.qhull import QhullError

from .ccvt import CCVTResult, nearest_available_assignment, split_capacity

__all__ = ["FastCCVT", "ccvt_fast", "recursive_split"]


def _ranges(starts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Concatenate arange(s, s + n) for each (s, n), vectorized."""
    lengths = np.asarray(lengths, dtype=np.int64)
    total = int(lengths.sum())
    if total == 0:
        return np.empty(0, dtype=np.int64)
    before = np.cumsum(lengths) - lengths
    return np.arange(total, dtype=np.int64) + np.repeat(np.asarray(starts, dtype=np.int64) - before, lengths)


def recursive_split(points: np.ndarray, capacities: np.ndarray,
                    order: np.ndarray | None = None, bounds: np.ndarray | None = None) -> np.ndarray:
    """Recursive median-cut initialization (Ahmed & Deussen 2017, section 3.4).

    Returns `slot_points`: a permutation of point indices such that site s owns
    slot_points[offsets[s] : offsets[s] + capacities[s]].

    Each slice is cut perpendicular to the major axis of its points, which keeps
    slices compact. All slices at the same depth are processed together.

    With groups: `bounds` gives the site index where each group starts (plus the
    total at the end), and `order` lists the point indices arranged so that each
    group's points fill that group's slots. Each group then starts as its own
    slice, so no site receives points from another group.
    """
    m, n = points.shape[0], capacities.shape[0]
    cum = np.concatenate([[0], np.cumsum(capacities)])
    order = np.arange(m, dtype=np.int64) if order is None else np.array(order, dtype=np.int64)
    if bounds is None:
        s0, s1 = np.array([0]), np.array([n])   # site ranges of the slices still to split
    else:
        bounds = np.asarray(bounds, dtype=np.int64)
        s0, s1 = bounds[:-1], bounds[1:]
        keep = (s1 - s0) > 1
        s0, s1 = s0[keep], s1[keep]

    while s0.size:
        seg_start, seg_len = cum[s0], cum[s1] - cum[s0]
        pos = _ranges(seg_start, seg_len)
        label = np.repeat(np.arange(s0.size), seg_len)
        xy = points[order[pos]]

        cnt = seg_len.astype(np.float64)
        mx = np.bincount(label, xy[:, 0], s0.size) / cnt
        my = np.bincount(label, xy[:, 1], s0.size) / cnt
        dx, dy = xy[:, 0] - mx[label], xy[:, 1] - my[label]
        sxx = np.bincount(label, dx * dx, s0.size)
        syy = np.bincount(label, dy * dy, s0.size)
        sxy = np.bincount(label, dx * dy, s0.size)
        theta = 0.5 * np.arctan2(2.0 * sxy, sxx - syy)     # major axis direction
        t = dx * np.cos(theta)[label] + dy * np.sin(theta)[label]

        perm = np.lexsort((t, label))                     # sort along the axis within each slice
        order[pos] = order[pos[perm]]

        mid = (s0 + s1) // 2
        c0 = np.stack([s0, mid], axis=1).ravel()
        c1 = np.stack([mid, s1], axis=1).ravel()
        keep = (c1 - c0) > 1
        s0, s1 = c0[keep], c1[keep]
    return order


class FastCCVT:
    """Fast discrete CCVT optimizer.

    Parameters
    ----------
    points : (m, 2) array
    capacities : (n,) int array, all positive, summing to m
    initial_sites : (n, 2) array, optional
        If given, points are first assigned greedily to their nearest site with
        free capacity (the baseline's rule). Otherwise the recursive median-cut
        initialization is used.
    rng : seed or Generator, used to order the batches of pairs
    site_groups, point_groups : (n,) and (m,) arrays of labels, optional
        Group constraint: a site only ever owns points with its own label. Each
        group's site capacities must add up to its number of points. Sites in
        different groups never exchange points.
    """

    SMALL_GROUP = 12      # groups this small get every pair of sites as neighbors

    def __init__(self, points, capacities, initial_sites=None, rng=None,
                 site_groups=None, point_groups=None):
        points = np.asarray(points, dtype=np.float64)
        capacities = np.asarray(capacities, dtype=np.int64)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("points must be an (m, 2) array")
        if capacities.sum() != points.shape[0]:
            raise ValueError("the sum of the site capacities must equal the number of points")
        if np.any(capacities <= 0):
            raise ValueError("all capacities must be positive")
        if (site_groups is None) != (point_groups is None):
            raise ValueError("give both site_groups and point_groups, or neither")

        # Work in coordinates centered on the data, for numerical stability with
        # large projected coordinates. Results are shifted back.
        self.origin = points.mean(axis=0)
        self.P = points - self.origin
        self.n = capacities.shape[0]
        self.rng = np.random.default_rng(rng)

        # Internally, sites are ordered by group; site_perm maps back.
        order = bounds = None
        self.site_group = None
        if site_groups is None:
            self.site_perm = np.arange(self.n)
        else:
            if initial_sites is not None:
                raise ValueError("initial_sites cannot be combined with groups")
            site_groups = np.asarray(site_groups)
            point_groups = np.asarray(point_groups)
            if site_groups.shape != (self.n,) or point_groups.shape != (points.shape[0],):
                raise ValueError("need one group label per site and per point")
            labels, sg = np.unique(site_groups, return_inverse=True)
            pg = np.searchsorted(labels, point_groups)
            if np.any(pg >= labels.size) or np.any(labels[np.minimum(pg, labels.size - 1)] != point_groups):
                raise ValueError("every point group must also have sites")
            need = np.bincount(sg, weights=capacities, minlength=labels.size)
            have = np.bincount(pg, minlength=labels.size)
            if not np.array_equal(need.astype(np.int64), have):
                raise ValueError("in every group, the site capacities must add up to its number of points")
            self.site_perm = np.argsort(sg, kind="stable")
            self.site_group = sg[self.site_perm]
            order = np.argsort(pg, kind="stable")
            bounds = np.concatenate([[0], np.cumsum(np.bincount(sg, minlength=labels.size))])

        capacities = capacities[self.site_perm]
        self.capacities = capacities
        self.offsets = np.concatenate([[0], np.cumsum(capacities)[:-1]])
        self.slot_site = np.repeat(np.arange(self.n), capacities)

        if initial_sites is None:
            self.slot_points = recursive_split(self.P, capacities, order=order, bounds=bounds)
        else:
            sites = np.asarray(initial_sites, dtype=np.float64) - self.origin
            if sites.shape != (self.n, 2):
                raise ValueError("need one initial site per capacity")
            owned = nearest_available_assignment(sites, capacities, self.P)
            self.slot_points = np.concatenate(owned)

        self.sites = self._centroids()
        self.active = np.ones(self.n, dtype=bool)
        self._last_edge_codes = np.empty(0, dtype=np.int64)
        self.iterations = 0
        self.moves_last = None

    # ---------------------------------------------------------------- helpers

    def _centroids(self) -> np.ndarray:
        xy = self.P[self.slot_points]
        c = self.capacities.astype(np.float64)
        return np.column_stack([
            np.bincount(self.slot_site, xy[:, 0], self.n) / c,
            np.bincount(self.slot_site, xy[:, 1], self.n) / c,
        ])

    def _neighbor_edges(self) -> np.ndarray:
        """Delaunay edges between sites, as (E, 2) with i < j."""
        n, S = self.n, self.sites
        if n < 2:
            return np.empty((0, 2), dtype=np.int64)
        if n == 2:
            e = np.array([[0, 1]], dtype=np.int64)
            return e if self.site_group is None or self.site_group[0] == self.site_group[1] \
                else np.empty((0, 2), dtype=np.int64)
        try:
            tri = Delaunay(S)
            simp = tri.simplices
            e = np.concatenate([simp[:, [0, 1]], simp[:, [1, 2]], simp[:, [0, 2]]])
            # Sites Qhull left out (e.g. coincident sites) are joined to the
            # vertices of the triangle nearest to them.
            if tri.coplanar.size:
                cp = tri.coplanar
                extra = np.concatenate([
                    np.column_stack([cp[:, 0], simp[cp[:, 1], k]]) for k in range(3)
                ])
                e = np.concatenate([e, extra])
        except QhullError:
            # Degenerate layouts (e.g. all sites collinear): use nearest neighbors.
            k = min(6, n - 1)
            _, nn = cKDTree(S).query(S, k=k + 1)
            e = np.column_stack([np.repeat(np.arange(n), k), nn[:, 1:].ravel()])
        # Qhull returns 32-bit indices; i * n + j below would overflow them
        # once n reaches the tens of thousands.
        e = np.sort(e.astype(np.int64), axis=1)
        e = e[e[:, 0] != e[:, 1]]
        if self.site_group is not None:
            e = e[self.site_group[e[:, 0]] == self.site_group[e[:, 1]]]
            e = np.concatenate([e, self._small_group_pairs()])
        codes = np.unique(e[:, 0] * n + e[:, 1])
        return np.column_stack([codes // n, codes % n])

    def _small_group_pairs(self) -> np.ndarray:
        """Every pair of sites within each small group. Delaunay edges alone can
        miss pairs when other groups' sites sit between them."""
        g = self.site_group
        starts = np.flatnonzero(np.r_[True, g[1:] != g[:-1]])
        sizes = np.diff(np.r_[starts, g.size])
        out = [np.empty((0, 2), dtype=np.int64)]
        for size in np.unique(sizes[(sizes >= 2) & (sizes <= self.SMALL_GROUP)]):
            a, b = np.triu_indices(size, k=1)
            st = starts[sizes == size][:, None]
            out.append(np.column_stack([(st + a).ravel(), (st + b).ravel()]))
        return np.concatenate(out).astype(np.int64)

    def _matchings(self, edges: np.ndarray):
        """Split edges into batches in which no site appears twice.

        Each round picks the edges whose random priority is the highest among
        all remaining edges at both of their endpoints (a Luby-style matching).
        """
        remaining = edges
        while remaining.shape[0]:
            pri = self.rng.random(remaining.shape[0])
            best = np.full(self.n, -1.0)
            np.maximum.at(best, remaining[:, 0], pri)
            np.maximum.at(best, remaining[:, 1], pri)
            sel = (pri == best[remaining[:, 0]]) & (pri == best[remaining[:, 1]])
            yield remaining[sel]
            remaining = remaining[~sel]

    def _process_batch(self, pairs: np.ndarray):
        """Optimally re-split the points of each (disjoint) pair of sites."""
        i, j = pairs[:, 0], pairs[:, 1]
        ci, cj = self.capacities[i], self.capacities[j]
        size = ci + cj
        b = pairs.shape[0]

        # Gather each pair's slots: site i's block followed by site j's block.
        slots = _ranges(np.stack([self.offsets[i], self.offsets[j]], 1).ravel(),
                        np.stack([ci, cj], 1).ravel())
        grp = np.repeat(np.arange(b), size)
        pos = np.arange(slots.size) - np.repeat(np.cumsum(size) - size, size)
        in_i = pos < ci[grp]                                  # slot belongs to site i

        pts = self.slot_points[slots]
        xy = self.P[pts]
        si, sj = self.sites[i], self.sites[j]
        d = si - sj
        mid = 0.5 * (si + sj)
        # key = |p - sj|^2 - |p - si|^2 : how much p prefers site i over site j
        key = 2.0 * ((xy[:, 0] - mid[grp, 0]) * d[grp, 0] + (xy[:, 1] - mid[grp, 1]) * d[grp, 1])

        # Within each pair, sort by descending key: the first c_i go to site i.
        perm = np.lexsort((-key, grp))
        now_in_i = np.empty_like(in_i)
        now_in_i[perm] = in_i                                 # new membership of each gathered point

        # Energy decrease for each pair; only apply strictly positive changes.
        gain = np.bincount(grp, key * (now_in_i.astype(np.float64) - in_i), b)
        scale = np.bincount(grp, np.abs(key), b)
        apply = gain > 1e-12 * np.maximum(scale, 1e-300)
        if not apply.any():
            return 0, np.empty(0, dtype=np.int64)

        new_pts = pts[perm]
        unchanged = ~apply[grp]
        new_pts[unchanged] = pts[unchanged]
        self.slot_points[slots] = new_pts

        moved = int(np.count_nonzero((now_in_i != in_i) & apply[grp]))

        # Move the affected sites to the centroids of their new points.
        touched = np.concatenate([i[apply], j[apply]])
        mask = apply[grp]
        site_of_slot = np.where(in_i, i[grp], j[grp])[mask]
        nxy = self.P[new_pts[mask]]
        sx = np.bincount(site_of_slot, nxy[:, 0], self.n)[touched]
        sy = np.bincount(site_of_slot, nxy[:, 1], self.n)[touched]
        c = self.capacities[touched].astype(np.float64)
        self.sites[touched, 0] = sx / c
        self.sites[touched, 1] = sy / c
        return moved, touched

    # ----------------------------------------------------------- optimization

    def iterate(self) -> int:
        """One pass over the neighbor pairs that need checking. Returns points moved."""
        edges = self._neighbor_edges()
        codes = edges[:, 0] * self.n + edges[:, 1]
        is_new = ~np.isin(codes, self._last_edge_codes, assume_unique=True)
        check = self.active[edges[:, 0]] | self.active[edges[:, 1]] | is_new
        self._last_edge_codes = codes

        next_active = np.zeros(self.n, dtype=bool)
        moved_total = 0
        for batch in self._matchings(edges[check]):
            moved, touched = self._process_batch(batch)
            moved_total += moved
            next_active[touched] = True
        self.active = next_active
        self.iterations += 1
        self.moves_last = moved_total
        return moved_total

    def run(self, max_iterations: int | None = 1000, tol: float | None = None, callback=None) -> bool:
        """Iterate until no point moves (converged), `tol` is met, or `max_iterations`.

        `tol` optionally stops early once the relative energy decrease of a
        pass falls below it. Returns True if fully converged (no moves).
        """
        prev = self.energy() if tol is not None else None
        while max_iterations is None or self.iterations < max_iterations:
            moved = self.iterate()
            if callback is not None:
                callback(self)
            if moved == 0:
                return True
            if tol is not None:
                e = self.energy()
                if prev - e <= tol * prev:
                    return False
                prev = e
        return False

    # ---------------------------------------------------------------- results

    def energy(self) -> float:
        """Sum of squared distances from every point to the site that owns it."""
        d = self.P[self.slot_points] - self.sites[self.slot_site]
        return float(np.einsum("ij,ij->", d, d))

    def site_locations(self) -> np.ndarray:
        """Site locations, in the order the sites were given."""
        out = np.empty_like(self.sites)
        out[self.site_perm] = self.sites + self.origin
        return out

    def assignment(self) -> np.ndarray:
        """For each point, the index of its site, in the order the sites were given."""
        owner = np.empty(self.P.shape[0], dtype=np.int64)
        owner[self.slot_points] = self.site_perm[self.slot_site]
        return owner


def ccvt_fast(points, n_sites: int, initial_sites=None, rng=None,
              max_iterations: int | None = 1000, tol: float | None = None,
              record_history: bool = False) -> CCVTResult:
    """Reduce `points` to `n_sites` blue-noise sites with the fast CCVT.

    Without `initial_sites`, the recursive median-cut initialization is used.
    """
    points = np.asarray(points, dtype=np.float64)
    capacities = split_capacity(points.shape[0], n_sites)
    opt = FastCCVT(points, capacities, initial_sites=initial_sites, rng=rng)
    history = []
    cb = (lambda o: history.append(o.energy())) if record_history else None
    converged = opt.run(max_iterations=max_iterations, tol=tol, callback=cb)
    return CCVTResult(
        sites=opt.site_locations(),
        capacities=capacities,
        assignment=opt.assignment(),
        iterations=opt.iterations,
        converged=converged,
        energy=opt.energy(),
        history=history,
    )
