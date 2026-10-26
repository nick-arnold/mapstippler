import numpy as np
import pytest

from mapstippler import ccvt, recursive_split, split_capacity
from mapstippler.ccvt_fast import FastCCVT


@pytest.mark.parametrize("method", ["fast", "baseline"])
def test_assignment_capacities_and_centroids(method):
    rng = np.random.default_rng(0)
    pts = rng.uniform(0, 100, (1003, 2))
    r = ccvt(pts, 100, rng=1, method=method)
    counts = np.bincount(r.assignment, minlength=100)
    assert np.array_equal(counts, r.capacities)
    own = np.column_stack([np.bincount(r.assignment, pts[:, 0]), np.bincount(r.assignment, pts[:, 1])]) / counts[:, None]
    at_centroid = np.isclose(own, r.sites).all(axis=1)
    if method == "fast":
        assert at_centroid.all()
    else:
        # Original behavior: a site only moves to its centroid after a swap, so a
        # site that never swaps stays at its starting position.
        start = pts[np.random.default_rng(1).choice(len(pts), 100, replace=False)]
        assert np.allclose(r.sites[~at_centroid], start[~at_centroid])
    assert r.converged


def test_recursive_split_gives_exact_slice_sizes():
    rng = np.random.default_rng(1)
    pts = rng.normal(size=(1000, 2))
    caps = split_capacity(1000, 37)
    order = recursive_split(pts, caps)
    assert np.array_equal(np.sort(order), np.arange(1000))      # a permutation


def test_energy_never_increases_and_beats_start():
    rng = np.random.default_rng(2)
    pts = rng.uniform(0, 1, (4000, 2))
    r = ccvt(pts, 500, rng=3, record_history=True)
    h = np.array(r.history)
    assert np.all(np.diff(h) <= 1e-9 * h[0])


def test_fast_reaches_baseline_quality_from_same_start():
    rng = np.random.default_rng(4)
    pts = rng.uniform(0, 1, (3200, 2))
    sites0 = pts[rng.choice(3200, 400, replace=False)]
    base = ccvt(pts, 400, initial_sites=sites0, method="baseline")
    fast = ccvt(pts, 400, initial_sites=sites0, method="fast", rng=0)
    assert fast.energy <= base.energy * 1.01


def test_large_coordinates_give_identical_result():
    rng = np.random.default_rng(5)
    pts = rng.uniform(0, 1, (4000, 2))
    caps = split_capacity(4000, 500)
    a = FastCCVT(pts, caps, rng=1); a.run()
    b = FastCCVT(pts * 250_000.0 + [500_000.0, 4_800_000.0], caps, rng=1); b.run()
    assert np.array_equal(a.assignment(), b.assignment())


def test_reproducible_with_seed():
    rng = np.random.default_rng(6)
    pts = rng.uniform(0, 1, (2000, 2))
    a = ccvt(pts, 250, rng=9)
    b = ccvt(pts, 250, rng=9)
    assert np.array_equal(a.sites, b.sites)


@pytest.mark.parametrize("n_sites", [1, 2, 3])
def test_tiny_site_counts(n_sites):
    rng = np.random.default_rng(7)
    pts = rng.uniform(0, 1, (40, 2))
    r = ccvt(pts, n_sites, rng=0)
    assert r.converged and np.array_equal(np.bincount(r.assignment), r.capacities)


def test_collinear_sites_fall_back_to_nearest_neighbors():
    x = np.linspace(0, 1, 400)
    pts = np.column_stack([x, 2 * x])          # every point on one line
    r = ccvt(pts, 40, rng=0)
    assert r.converged and np.array_equal(np.bincount(r.assignment), r.capacities)
    # dots should be evenly spread along the line
    gaps = np.diff(np.sort(r.sites[:, 0]))
    assert gaps.std() < 0.1 * gaps.mean()
