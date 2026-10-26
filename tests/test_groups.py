import numpy as np
import pytest

from mapstippler.ccvt_fast import FastCCVT


def test_sites_only_own_points_from_their_group():
    rng = np.random.default_rng(0)
    # two interleaved groups over the same square: the hard case for Delaunay edges
    pts = rng.uniform(0, 1, (3000, 2))
    pg = (pts[:, 0] * 10).astype(int) % 2               # vertical stripes
    n0, n1 = np.bincount(pg)
    caps = np.r_[np.full(n0 // 10, 10), np.full(n1 // 10, 10)]
    extra0, extra1 = n0 % 10, n1 % 10
    caps = np.r_[caps, np.ones(extra0 + extra1, dtype=int)]
    sg = np.r_[np.zeros(n0 // 10), np.ones(n1 // 10), np.zeros(extra0), np.ones(extra1)].astype(int)
    opt = FastCCVT(pts, caps, rng=1, site_groups=sg, point_groups=pg)
    assert opt.run()
    owner = opt.assignment()
    assert np.array_equal(sg[owner], pg)                           # never crosses groups
    assert np.array_equal(np.bincount(owner, minlength=caps.size), caps)


def test_group_capacity_mismatch_is_rejected():
    pts = np.zeros((5, 2))
    with pytest.raises(ValueError, match="add up"):
        FastCCVT(pts, [3, 2], site_groups=[0, 1], point_groups=[0, 0, 0, 0, 1])


def test_many_small_groups_are_independent():
    rng = np.random.default_rng(2)
    centers = rng.uniform(0, 1000, (300, 2))
    sizes = rng.integers(10, 40, 300)
    pts = np.vstack([c + rng.normal(0, 2, (k, 2)) for c, k in zip(centers, sizes)])
    pg = np.repeat(np.arange(300), sizes)
    caps, sg = [], []
    for g, k in enumerate(sizes):
        caps += [10] * (k // 10) + [1] * (k % 10)
        sg += [g] * (k // 10 + k % 10)
    opt = FastCCVT(pts, caps, rng=3, site_groups=sg, point_groups=pg)
    assert opt.run()
    sites = opt.site_locations()
    big = np.asarray(caps) == 10
    # every large site stays inside its own cluster, never between clusters
    d = np.linalg.norm(sites[big] - centers[np.asarray(sg)[big]], axis=1)
    assert d.max() < 10
