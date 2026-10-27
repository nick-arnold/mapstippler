import json

import numpy as np
import pytest
import shapely

gpd = pytest.importorskip("geopandas")

from mapstippler.graduated import make_dot_map
from mapstippler.output import render_preview, write_outputs
from mapstippler.sampling import random_points_in_polygons
from mapstippler.scale import MapScale
from mapstippler.__main__ import main


def synthetic_grid(n=12, cell=10_000.0):
    """n x n squares in North America Albers, a dense center fading outward."""
    polys, vals, ids = [], [], []
    c = (n - 1) / 2
    for i in range(n):
        for j in range(n):
            x, y = (i - c) * cell, (j - c) * cell
            polys.append(shapely.box(x, y, x + cell, y + cell))
            r2 = ((i - c) ** 2 + (j - c) ** 2) / (n / 4) ** 2
            vals.append(int(400_000 * np.exp(-r2) + 800))
            ids.append(f"{i:02d}{j:02d}")
    return gpd.GeoDataFrame({"gid": ids, "pop": vals}, geometry=polys, crs="ESRI:102008")


SCALE = MapScale(250.0, sizes_px=(2.0, 3.5, 6.0))   # pinned, not the library defaults:
# small dot = 500 m, so the center coalesces at 1,000 per dot


def test_random_points_land_inside_their_polygons():
    ring = shapely.Point(0, 0).buffer(1000).difference(shapely.Point(0, 0).buffer(500))
    geoms = np.array([ring, shapely.box(3000, 0, 4000, 5000)], dtype=object)
    xy, owner = random_points_in_polygons(geoms, [500, 300], rng=0)
    assert xy.shape == (800, 2) and np.array_equal(np.bincount(owner), [500, 300])
    for k, g in enumerate(geoms):
        assert shapely.contains_xy(g, xy[owner == k, 0], xy[owner == k, 1]).all()


@pytest.mark.parametrize("method", ["random", "blue", "graduated"])
def test_pipeline_conserves_value_and_accounts_for_every_dot(method):
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", method=method, rng=1, progress=None)
    n_dots = int(round(g["pop"].sum() / 1000))
    assert r.dots["value"].sum() == pytest.approx(n_dots * 1000)
    assert r.summary["value_of_all_dots"] == pytest.approx(n_dots * 1000)
    assert set(r.dots["source_id"]) <= set(g["gid"])
    assert 0 <= r.summary["per_dot"]["purity"]["min"] <= 1
    if method == "graduated":
        assert r.dots["class"].max() >= 2                     # something got promoted
        for cls, v in ((1, 1_000), (2, 10_000), (3, 100_000)):
            assert (r.dots.loc[r.dots["class"] == cls, "value"] == v).all()   # exact, no variation
    if method != "graduated":
        assert (r.dots["class"] == 1).all()


def test_touching_small_dots_are_only_leftovers_or_small_clusters():
    # Per-cluster promotion keeps values exact, so a cluster of 14 touching dots
    # becomes one larger dot and 4 leftovers, and clusters under 10 stay as they
    # are. Any small dots still touching must come from those two groups.
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, method="graduated", rng=2, progress=None)
    p1 = r.summary["promotions"][0]
    still_touching = r.summary["overlapping_dots"]["class1_with_class1"]
    assert still_touching <= p1["leftovers"] + p1["dots_in_small_clusters"]
    assert p1["new_dots"] > 0


def test_outputs_and_preview(tmp_path):
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=3, progress=None)
    dots_path, rec_path, info_path = write_outputs(r, tmp_path)
    back = gpd.read_file(dots_path, layer="dots")
    assert len(back) == len(r.dots) and back.crs == g.crs
    info = json.loads(info_path.read_text())
    assert info["params"]["search_radii_m"][0] == pytest.approx(500.0)
    png, (w, h) = render_preview(r.dots, g, SCALE, tmp_path / "p.png")
    assert png.exists() and w == pytest.approx(120_000 / 250 + 20, abs=1)


def test_cli_runs_end_to_end(tmp_path, capsys):
    g = synthetic_grid()
    src = tmp_path / "grid.parquet"
    g.to_parquet(src)
    main(["advise", str(src), "--value", "pop", "--mpp", "250"])
    assert "suggested" in capsys.readouterr().out
    out = tmp_path / "run"
    main(["dots", str(src), "--value", "pop", "--unit", "1000", "--mpp", "250",
          "--id", "gid", "--out", str(out), "--preview"])
    assert (out / "dots.gpkg").exists() and (out / "preview.png").exists()


def test_rejects_geographic_coordinates():
    g = synthetic_grid().to_crs(4326)
    with pytest.raises(ValueError, match="projected"):
        make_dot_map(g, "pop", 1000, SCALE, progress=None)


def test_distant_clusters_never_merge_into_a_dot_between_them():
    # Two small, dense towns 200 km apart, with nothing in between. Each has
    # enough touching small dots for one or two medium dots plus leftovers.
    # Global pooling used to combine their remainders into a dot in the empty
    # space between them; per-cluster promotion must never do that.
    town_a = shapely.box(0, 0, 3_000, 3_000)
    town_b = shapely.box(200_000, 0, 203_000, 3_000)
    empty = shapely.box(-20_000, -20_000, 223_000, 23_000).difference(
        shapely.MultiPolygon([town_a, town_b]))
    g = gpd.GeoDataFrame({"gid": ["a", "b", "empty"], "pop": [17_000, 14_000, 0]},
                         geometry=[town_a, town_b, empty], crs="ESRI:102008")
    r = make_dot_map(g, "pop", 1_000, SCALE, id_column="gid", rng=4, progress=None)
    x = r.dots.geometry.x.to_numpy()
    assert np.all((x < 10_000) | (x > 193_000))                     # nothing in the gap
    assert (r.dots.loc[r.dots["class"] == 2, "value"] == 10_000).all()
    assert r.dots["value"].sum() == 31_000
    assert r.summary["per_dot"]["spread_px"]["max_px"] < 20          # no dot stands for far-away points


def test_nudging_the_top_class_resolves_its_overlaps():
    g = synthetic_grid()
    before = make_dot_map(g, "pop", 150, SCALE, rng=1, progress=None)
    after = make_dot_map(g, "pop", 150, SCALE, rng=1, nudge="top", progress=None)
    assert before.summary["overlapping_dots"]["class3_with_class3"] > 0
    assert after.summary["overlapping_dots"]["class3_with_class3"] == 0
    d = after.dots
    moved = after.summary["per_dot"]["nudge_px"]
    assert moved["max_px"] <= 6.0 + 1e-6                                 # at most one 6 px diameter
    assert d["value"].sum() == before.dots["value"].sum()               # values untouched
    assert after.summary["nudge"]["overlapping_after"] == 0


def test_preview_can_draw_another_outline_layer(tmp_path):
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, rng=3, progress=None)
    outline = gpd.GeoDataFrame(geometry=[g.union_all()], crs=g.crs)
    png, _ = render_preview(r.dots, g, SCALE, tmp_path / "p.png", outlines=outline)
    assert png.exists()


def test_single_class_blue_noise_map_and_advice(tmp_path, capsys):
    g = synthetic_grid()
    src = tmp_path / "grid.parquet"
    g.to_parquet(src)
    main(["advise", str(src), "--value", "pop", "--mpp", "250", "--method", "blue"])
    out = capsys.readouterr().out
    assert "dots merge" in out and "suggested" in out and "class 2" not in out
    run = tmp_path / "blue"
    main(["dots", str(src), "--value", "pop", "--unit", "2000", "--mpp", "250", "--method", "blue",
          "--sizes", "2", "--out", str(run)])
    d = gpd.read_file(run / "dots.gpkg", layer="dots")
    assert (d["class"] == 1).all() and (d["value"] == 2000).all() and (d["size_px"] == 2).all()


def test_search_radius_can_be_set_independently_of_dot_size():
    g = synthetic_grid()
    small = MapScale(250.0, sizes_px=(2.0, 3.5, 6.0))
    wide = MapScale(250.0, sizes_px=(2.0, 3.5, 6.0), radii_px=(3.0, 5.0, 8.0))
    a = make_dot_map(g, "pop", 1000, small, rng=1, progress=None)
    b = make_dot_map(g, "pop", 1000, wide, rng=1, progress=None)
    # a wider radius counts more dots as merged, so more get promoted
    assert b.summary["promotions"][0]["touching"] > a.summary["promotions"][0]["touching"]
    assert (b.dots["size_px"].unique() == a.dots["size_px"].unique()).all()   # drawn the same
    assert b.params["radii_px"] == [3.0, 5.0, 8.0]
    assert b.params["sizes_px"] == [2.0, 3.5, 6.0]


def test_single_radius_value_acts_as_a_multiplier():
    s = MapScale(100.0, sizes_px=(2.0, 4.0), radii_px=(1.5,))
    assert [s.radius_px(i) for i in (0, 1)] == [3.0, 6.0]


def test_radius_under_dot_size_warns():
    g = synthetic_grid()
    tight = MapScale(250.0, sizes_px=(2.0, 3.5, 6.0), radii_px=(1.0, 3.5, 6.0))
    with pytest.warns(UserWarning, match="smaller than the dot diameter"):
        make_dot_map(g, "pop", 1000, tight, rng=1, progress=None)


def test_dots_state_object_tracks_points_through_promotion():
    # every input point must still map to a real dot after promotion, and the
    # values of the dots it maps to must add up to the total
    from mapstippler.graduated import _Dots
    import numpy as np

    d = _Dots(xy=np.zeros((4, 2)), cls=np.zeros(4, dtype=np.int64),
              value=np.full(4, 10.0), point_dot=np.array([0, 0, 1, 2, 3]))
    base = d.add(np.zeros((1, 2)), 1, 40.0)          # one big dot worth four small
    assert base == 4 and d.n == 5
    d.reassign(np.array([0, 1, 2, 3]), np.array([4, 4, 4, 4]))
    assert np.all(d.point_dot == 4)
    d.keep_only(np.array([4]))
    assert d.n == 1 and np.all(d.point_dot == 0) and d.value[0] == 40.0


def test_placement_accuracy_is_reported_separately_from_promotion():
    """Before promotion every dot is worth one unit, so a dot crossing a boundary
    moves one unit. After promotion a single dot can carry 100 units across, which
    is a property of the method, not of how well the dots were placed. The summary
    must separate the two."""
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=1, progress=None)
    acc = r.summary["accuracy"]
    assert acc["placement"] is not None and acc["final"] is not None

    rec = r.reconciliation
    assert "difference_before_promotion" in rec
    # placement error is bounded: a polygon can only lose or gain whole unit-value
    # dots, so the worst case stays modest
    assert acc["placement"]["max_dots"] < acc["final"]["max_dots"]
    assert acc["placement"]["median_dots"] <= acc["final"]["median_dots"] + 1e-9


def test_before_promotion_column_accounts_for_the_whole_total():
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=2, progress=None)
    rec = r.reconciliation
    # every dot is inside some polygon in this synthetic grid, so the
    # pre-promotion values must add up to the value of all the dots
    assert rec["value_inside_before_promotion"].sum() == pytest.approx(
        r.summary["value_of_all_dots"], rel=1e-9)


def test_single_class_maps_have_no_separate_placement_figure():
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, method="blue", rng=3, progress=None)
    # nothing is promoted, so placement and final are the same thing
    assert r.summary["accuracy"]["placement"] is None
    assert "difference_before_promotion" not in r.reconciliation


def test_geopackage_is_the_default_and_round_trips(tmp_path):
    """QGIS reads GeoPackage everywhere; it reads GeoParquet only when its GDAL
    was built with Arrow, which many builds are not."""
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=1, progress=None)
    dots_path, _, _ = write_outputs(r, tmp_path)
    assert dots_path.name == "dots.gpkg"
    back = gpd.read_file(dots_path, layer="dots")
    assert len(back) == len(r.dots)
    assert back.crs == g.crs
    for col in ("class", "value", "size_px", "source_id"):
        assert col in back.columns


def test_categories_survive_the_geopackage_write(tmp_path):
    """GeoPackage has no categorical type, so the labels are written as text."""
    g = synthetic_grid()
    g["a"] = (g["pop"] * 0.6).astype(int)
    g["b"] = g["pop"] - g["a"]
    r = make_dot_map(g, ["a", "b"], 1000, SCALE, id_column="gid", rng=2, progress=None)
    dots_path, _, _ = write_outputs(r, tmp_path)
    back = gpd.read_file(dots_path, layer="dots")
    assert set(back["category"].unique()) == {"a", "b"}


def test_geopackage_carries_a_qgis_style(tmp_path):
    """The layer should open already symbolised, with dots at the size they were
    designed for. QGIS reads styles from a layer_styles table in the package."""
    import sqlite3
    import xml.etree.ElementTree as ET

    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=1, progress=None)
    dots_path, _, _ = write_outputs(r, tmp_path, scale=SCALE)

    con = sqlite3.connect(dots_path)
    rows = con.execute("SELECT f_table_name, useAsDefault, styleQML FROM layer_styles").fetchall()
    con.close()
    assert len(rows) == 1 and rows[0][0] == "dots" and rows[0][1] == 1

    qml = rows[0][2]
    root = ET.fromstring(qml[qml.index("<qgis"):])          # must be valid XML
    rules = root.findall(".//rule")
    assert len(rules) == len(r.dots["class"].unique())

    # sizes are in map units, equal to the dot's CSS pixels times meters per pixel
    sizes = sorted(float({o.get("name"): o.get("value")
                          for o in s.findall(".//Option[@name]")}["size"])
                   for s in root.findall(".//symbol"))
    expected = sorted(SCALE.sizes_px[c - 1] * SCALE.m_per_px
                      for c in sorted(r.dots["class"].unique()))
    assert sizes == pytest.approx(expected)
    assert (tmp_path / "dots.qml").exists()                  # and alongside, for other tools


def test_style_has_one_rule_per_class_and_category():
    from mapstippler.output import style_qml
    import xml.etree.ElementTree as ET

    g = synthetic_grid()
    g["a"] = (g["pop"] * 0.6).astype(int)
    g["b"] = g["pop"] - g["a"]
    r = make_dot_map(g, ["a", "b"], 1000, SCALE, id_column="gid", rng=2, progress=None)
    qml = style_qml(r, SCALE)
    root = ET.fromstring(qml[qml.index("<qgis"):])
    rules = root.findall(".//rule")
    assert len(rules) == len(r.dots["class"].unique()) * 2
    colours = {o.get("value") for s in root.findall(".//symbol")
               for o in s.findall(".//Option[@name='color']")}
    assert len(colours) == 2                                 # one per category


def test_dot_layer_carries_only_what_a_map_needs():
    """Diagnostics belong in the run record, not on every dot."""
    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=1, progress=None)
    assert set(r.dots.columns) == {"dot_id", "class", "value", "size_px",
                                   "source_id", "geometry"}
    per_dot = r.summary["per_dot"]
    for key in ("spread_px", "displacement_px", "purity"):
        assert key in per_dot


def test_geopackage_carries_a_qgis_style(tmp_path):
    """The layer should open already symbolised, with dots at the size they were
    designed for: size_px x m_per_px, in map units."""
    import sqlite3
    import xml.etree.ElementTree as ET

    g = synthetic_grid()
    r = make_dot_map(g, "pop", 1000, SCALE, id_column="gid", rng=1, progress=None)
    dots_path, _, _ = write_outputs(r, tmp_path, scale=SCALE)

    con = sqlite3.connect(dots_path)
    try:
        row = con.execute("SELECT f_table_name, useAsDefault, styleQML "
                          "FROM layer_styles").fetchone()
    finally:
        con.close()
    assert row[0] == "dots" and row[1] == 1

    qml = row[2]
    root = ET.fromstring(qml[qml.index("<qgis"):])          # must be well-formed
    sizes = sorted({float(o.get("value"))
                    for o in root.findall(".//Option[@name='size']")})
    expected = sorted({s * SCALE.m_per_px for s in SCALE.sizes_px[:r.dots["class"].max()]})
    assert sizes == pytest.approx(expected)
    units = {o.get("value") for o in root.findall(".//Option[@name='size_unit']")}
    assert units == {"MapUnit"}
    assert (tmp_path / "dots.qml").exists()
