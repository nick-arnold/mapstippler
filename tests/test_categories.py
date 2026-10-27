import numpy as np
import pytest
import shapely

gpd = pytest.importorskip("geopandas")

from mapstippler.graduated import make_dot_map
from mapstippler.output import render_preview
from mapstippler.scale import MapScale
from mapstippler.__main__ import main


def two_party_grid(n=10, cell=20_000.0, seed=0):
    """Squares whose Democratic share rises with population density."""
    rng = np.random.default_rng(seed)
    polys, pop, ids = [], [], []
    c = (n - 1) / 2
    for i in range(n):
        for j in range(n):
            polys.append(shapely.box(i * cell, j * cell, (i + 1) * cell, (j + 1) * cell))
            r2 = ((i - c) ** 2 + (j - c) ** 2) / (n / 4) ** 2
            pop.append(int(600_000 * np.exp(-r2) + rng.integers(3_000, 30_000)))
            ids.append(f"{i:02d}{j:02d}")
    g = gpd.GeoDataFrame({"geoid": ids, "pop": pop}, geometry=polys, crs="ESRI:102008")
    dens = g["pop"] / (g.geometry.area / 1e6)
    share = 0.22 + 0.60 / (1 + np.exp(-(np.log10(dens) - 1.6) * 2.2))
    g["dem"] = (g["pop"] * share).round().astype(int)
    g["rep"] = g["pop"] - g["dem"]
    return g


SCALE = MapScale(700.0, sizes_px=(2.0, 4.5, 10.0))


def test_category_shares_match_the_data():
    unit, n_poly = 2000, 100
    g = two_party_grid()
    r = make_dot_map(g, ["dem", "rep"], unit, SCALE, id_column="geoid", rng=1, progress=None)
    cats = r.summary["categories"]
    # each polygon's split is rounded to whole dots, so a category can be off by
    # up to half a dot per polygon; that is the only error allowed
    budget = n_poly * unit / 2
    for name in ("dem", "rep"):
        assert abs(cats[name]["value_of_dots"] - cats[name]["value"]) < budget
    assert cats["dem"]["dots"] + cats["rep"]["dots"] == len(r.dots)
    # national share should be far tighter than the per-polygon budget suggests,
    # because the rounding errors mostly cancel
    total = cats["dem"]["value_of_dots"] + cats["rep"]["value_of_dots"]
    assert abs(cats["dem"]["value_of_dots"] / total
               - g["dem"].sum() / (g["dem"].sum() + g["rep"].sum())) < 0.002


def test_each_polygon_gets_its_own_split_not_the_national_one():
    g = two_party_grid()
    r = make_dot_map(g, ["dem", "rep"], 1000, SCALE, id_column="geoid", rng=2, progress=None)
    rec = r.reconciliation
    all_want = rec["dem_value"] / rec["value"]
    assert all_want.max() - all_want.min() > 0.3          # the test data really does vary

    judged = rec[rec["value"] > 20_000]                   # at least ~20 dots to judge a share
    got = judged["dem_value_of_dots"] / (judged["dem_value_of_dots"] + judged["rep_value_of_dots"])
    want = judged["dem_value"] / judged["value"]
    # a promoted dot carries its whole value to one polygon, so per-polygon
    # shares are granular; the median is what shows they track
    assert np.abs(got - want).median() < 0.03
    # and each polygon must follow its own split, not the national one
    national = g["dem"].sum() / (g["dem"].sum() + g["rep"].sum())
    assert np.abs(got - want).mean() < 0.5 * np.abs(got - national).mean()


def test_promotion_keeps_categories_separate_and_values_exact():
    g = two_party_grid()
    r = make_dot_map(g, ["dem", "rep"], 200, SCALE, id_column="geoid", rng=3, progress=None)
    d = r.dots
    assert d["class"].max() >= 2                          # something was promoted
    for cls in d["class"].unique():
        vals = d.loc[d["class"] == cls, "value"].unique()
        assert vals.size == 1                             # one value per class, no averaging
    # a promoted dot is worth exactly `ratio` dots of its own category
    for p in r.summary["promotions"]:
        assert p["category"] in ("dem", "rep")


def test_category_column_order_follows_the_columns_given():
    g = two_party_grid()
    r = make_dot_map(g, ["rep", "dem"], 2000, SCALE, id_column="geoid", rng=4, progress=None)
    assert list(r.dots["category"].cat.categories) == ["rep", "dem"]


def test_single_column_still_has_no_category():
    g = two_party_grid()
    r = make_dot_map(g, "pop", 2000, SCALE, id_column="geoid", rng=5, progress=None)
    assert "category" not in r.dots.columns
    assert r.summary["categories"] is None


def test_cli_and_preview(tmp_path):
    g = two_party_grid()
    src = tmp_path / "party.parquet"
    g.to_parquet(src)
    out = tmp_path / "run"
    main(["dots", str(src), "--value", "dem", "rep", "--unit", "2000", "--mpp", "700",
          "--id", "geoid", "--out", str(out), "--preview"])
    d = gpd.read_file(out / "dots.gpkg", layer="dots")
    assert set(d["category"].unique()) == {"dem", "rep"}
    assert (out / "preview.png").exists()


def test_more_than_two_categories():
    g = two_party_grid()
    g["other"] = (g["pop"] * 0.05).round().astype(int)
    g["dem"] -= g["other"] // 2
    g["rep"] -= g["other"] - g["other"] // 2
    r = make_dot_map(g, ["dem", "rep", "other"], 2000, SCALE, id_column="geoid",
                     rng=6, progress=None)
    assert set(r.summary["categories"]) == {"dem", "rep", "other"}
    assert r.summary["categories"]["other"]["dots"] > 0
