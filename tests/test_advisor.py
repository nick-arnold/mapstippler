import numpy as np
import pytest
import shapely

from mapstippler.advisor import advise, coalescing_share, density_sample
from mapstippler.scale import MapScale


def test_uniform_square_matches_hand_calculation():
    # 100 km square, 1,000,000 people, dot diameter 1,000 m.
    # d / s = d * sqrt(density / u) = 1000 * sqrt(1e-4 / u) = 10 / sqrt(u),
    # so a class coalesces (d / s >= 0.8) exactly when u <= 156.25.
    sq = np.array([shapely.box(0, 0, 100_000, 100_000)], dtype=object)
    xy, w = density_sample(sq, [1_000_000], n_sample=400_000, rng=0)
    assert coalescing_share(xy, w, unit=100, diameter=1000) > 0.9
    assert coalescing_share(xy, w, unit=250, diameter=1000) < 0.01


def test_tiny_dense_polygons_are_judged_with_their_neighbors():
    # A grid of 400 m cells, each far smaller than a 3 km dot. Per polygon they
    # look extremely dense; over dot-sized neighborhoods the density is modest.
    cells = [shapely.box(x, y, x + 400, y + 400) for x in range(0, 40_000, 400) for y in range(0, 40_000, 400)]
    vals = np.full(len(cells), 10.0)                 # 100,000 people over 1,600 km2
    scale = MapScale(1500.0)                         # class 1 = 3 km, 2 = 5.25 km, 3 = 9 km
    df = advise(np.array(cells, dtype=object), vals, scale, unit_values=[1000, 2000, 5000],
                n_sample=200_000)
    # class 1 at u = 1000: d / s = 3000 * sqrt(62.5e-6 / 1000) = 0.75 -> not coalescing
    assert df.loc[df["unit"] == 1000, "c1_coal"].iloc[0] < 0.05
    assert df["suggested"].sum() == 1


def test_ink_kept_matches_the_formula_and_is_reciprocal_to_density():
    from mapstippler.advisor import ink_report

    s = MapScale(100.0, sizes_px=(2.0, 4.5, 10.0), class_ratio=10.0)
    rep = ink_report(s)
    assert rep[0]["ink_kept"] == pytest.approx((4.5 / 2.0) ** 2 / 10)
    for r in rep:
        assert r["ink_kept"] * r["density_absorbed"] == pytest.approx(1.0)


def test_area_proportional_sizes_keep_all_the_ink_and_absorb_nothing():
    # diameter ratio sqrt(10) is the area-proportional case
    s = MapScale(100.0, sizes_px=(2.0, 2.0 * 10 ** 0.5), class_ratio=10.0)
    assert s.ink_kept(1) == pytest.approx(1.0)


def test_first_class_has_nothing_to_promote_from():
    assert MapScale(100.0).ink_kept(0) == 1.0


def test_single_class_keeps_the_search_radius():
    s = MapScale(100.0, sizes_px=(2.0, 4.5, 10.0), radii_px=(3.0, 6.0, 12.0)).single_class()
    assert s.n_classes == 1
    assert s.sizes_px == (2.0,)
    assert s.radius_px(0) == 3.0            # not silently reset to the diameter


def test_single_class_without_radii_still_uses_diameters():
    s = MapScale(100.0).single_class()
    assert s.radii_px is None and s.radius_px(0) == s.sizes_px[0]
