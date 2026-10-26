import numpy as np
import pytest

from mapstippler.ccvt import CCVTOptimizer, ccvt, reduce_points, split_capacity


def test_split_capacity_differs_by_at_most_one():
    caps = split_capacity(47, 5)
    assert caps.sum() == 47
    assert caps.max() - caps.min() <= 1
    assert list(caps) == [9, 9, 9, 10, 10]


def test_capacities_must_match_point_count():
    with pytest.raises(ValueError):
        CCVTOptimizer(np.zeros((2, 2)), [1, 1], np.zeros((3, 2)))


def test_every_point_assigned_once_and_capacities_respected():
    rng = np.random.default_rng(0)
    pts = rng.uniform(0, 100, (403, 2))
    r = ccvt(pts, 40, rng=1)
    counts = np.bincount(r.assignment, minlength=40)
    assert np.array_equal(counts, r.capacities)
    assert r.converged


def test_sites_sit_at_centroids_of_their_points():
    rng = np.random.default_rng(1)
    pts = rng.uniform(0, 100, (300, 2))
    r = ccvt(pts, 30, rng=2)
    for s in range(30):
        own = pts[r.assignment == s]
        assert np.allclose(r.sites[s], own.mean(axis=0))


def test_reproducible_with_seed():
    rng = np.random.default_rng(2)
    pts = rng.uniform(0, 1, (500, 2))
    a = reduce_points(pts, 2, rng=5)
    b = reduce_points(pts, 2, rng=5)
    assert np.array_equal(a.sites, b.sites)


def test_energy_does_not_increase():
    rng = np.random.default_rng(3)
    pts = rng.uniform(0, 1, (800, 2))
    r = ccvt(pts, 50, rng=4, record_history=True)
    h = np.array(r.history)
    assert np.all(np.diff(h) <= 1e-9 * h[0])
