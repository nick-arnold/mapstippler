import numpy as np

from mapstippler.nudge import nudge


def min_gap(X, r):
    d = np.hypot(*(X[:, None, :] - X[None, :, :]).transpose(2, 0, 1))
    np.fill_diagonal(d, np.inf)
    return (d - 2 * r).min()


def test_two_overlapping_dots_split_the_difference():
    res = nudge([[0, 0], [1.0, 0]], radius=1.0)            # need 2 apart, 1 apart now
    assert res.overlapping_before == 2 and res.overlapping_after == 0
    assert np.allclose(res.xy, [[-0.5, 0], [1.5, 0]], atol=1e-6)


def test_fixed_dots_never_move():
    res = nudge([[0, 0], [1.0, 0]], radius=1.0, movable=[False, True])
    assert np.allclose(res.xy[0], [0, 0]) and np.allclose(res.xy[1], [2, 0], atol=1e-6)


def test_heavier_dots_move_less():
    res = nudge([[0, 0], [1.0, 0]], radius=1.0, weight=[3.0, 1.0])
    assert np.allclose(res.moved, [0.25, 0.75], atol=1e-6)


def test_crowded_cluster_resolves_when_there_is_room():
    rng = np.random.default_rng(0)
    X = rng.normal(0, 3.0, (40, 2))                        # 40 dots of radius 1, badly overlapping
    res = nudge(X, radius=1.0)
    assert res.overlapping_before > 20 and res.overlapping_after == 0
    assert min_gap(res.xy, 1.0) > -1e-5            # within the 1e-6 relative tolerance


def test_cap_is_respected_and_leftover_overlap_is_reported():
    rng = np.random.default_rng(1)
    X = rng.normal(0, 1.0, (40, 2))                        # far too crowded to fix within the cap
    res = nudge(X, radius=1.0, max_disp=1.0)
    assert res.moved.max() <= 1.0 + 1e-9
    assert res.overlapping_after > 0 and res.capped > 0


def test_nothing_moves_without_overlap():
    X = np.array([[0, 0], [5, 0], [0, 5]], dtype=float)
    res = nudge(X, radius=1.0)
    assert res.iterations == 0 and np.array_equal(res.xy, X)
