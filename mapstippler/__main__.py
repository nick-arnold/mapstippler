"""Command line for mapstippler.

    mapstippler advise county_pop.parquet --value population
    mapstippler dots county_pop.parquet --value population --unit 10000 --id geoid --preview

Also runnable as `python -m mapstippler ...`.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path


def _value_arg(args):
    """One column name, or a list of them for a categorised map."""
    return args.value[0] if len(args.value) == 1 else list(args.value)


def _read(path, where):
    import geopandas as gpd
    p = Path(path)
    g = gpd.read_parquet(p) if p.suffix.lower() in (".parquet", ".geoparquet") else gpd.read_file(p)
    if where:
        g = g.query(where)
    return g.reset_index(drop=True)


def _scale(args):
    from .scale import MapScale
    kw = {"sizes_px": tuple(args.sizes), "class_ratio": args.ratio,
          "radii_px": tuple(args.radii) if args.radii else None}
    return MapScale(args.mpp, **kw) if args.mpp else MapScale.zoom(args.zoom, **kw)


def cmd_advise(args):
    from .advisor import advise, ink_report

    t = time.perf_counter()
    g = _read(args.input, args.where)
    scale = _scale(args)
    single = args.method == "blue"
    if single:
        scale = scale.single_class()
    value = _value_arg(args)
    names = [value] if isinstance(value, str) else value
    vals = g[names].sum(axis=1).to_numpy(dtype=float) if len(names) > 1 \
        else g[names[0]].to_numpy(dtype=float)
    label = names[0] if len(names) == 1 else " + ".join(names)
    print(f"{len(g):,} polygons, total {label} {vals.sum():,.0f}  (read in {time.perf_counter() - t:.1f} s)")
    d = scale.describe()
    print(f"scale: {scale.m_per_px:,.0f} m per CSS pixel (about 1:{d['approx_scale_denominator']:,})")
    print(("dot size: " if single else "dot classes: ")
          + ", ".join(f"{s} px = {m:,.0f} m" for s, m in zip(d["sizes_px"], d["diameters_m"])))
    if not d["radii_are_diameters"]:
        print("search radii: " + ", ".join(f"{r} px = {m:,.0f} m" for r, m in zip(d["radii_px"], d["search_radii_m"])))
    for r in ink_report(scale):
        note = ("holes where it promotes" if r["ink_kept"] < 0.35 else
                "absorbs little extra density" if r["ink_kept"] > 0.8 else "balanced")
        print(f"promotion into class {r['class']}: {r['from_px']}→{r['to_px']} px keeps "
              f"{r['ink_kept']:.0%} of the ink, absorbs {r['density_absorbed']:.1f}x the density ({note})")
    t = time.perf_counter()
    df = advise(g.geometry.values, vals, scale, unit_values=args.units, k0=args.k0,
                max_top_class_coalesce=args.max_top)
    print(f"density measured over neighborhoods 3 dot diameters wide  ({time.perf_counter() - t:.1f} s)")
    print()

    n = scale.n_classes
    if single:
        print(f"{'unit':>9} {'dots':>10}  {'dots merge':>10}")
        for _, r in df.iterrows():
            mark = "  <- suggested" if r["suggested"] else ""
            print(f"{r['unit']:>9,.0f} {int(r['dots']):>10,}  {r['c1_coal']:>10.1%}{mark}")
    else:
        print(f"{'unit':>9}  " + "  ".join(f"{'class ' + str(c):>10}" for c in range(1, n + 1))
              + f"  {'still':>8}")
        print(f"{'':>9}  " + "  ".join(f"{'dots':>10}" for _ in range(1, n + 1))
              + f"  {'crowded':>8}")
        for _, r in df.iterrows():
            cells = "  ".join(f"{int(r[f'est{c}']):>10,}" for c in range(1, n + 1))
            mark = "  <- suggested" if r["suggested"] else ""
            print(f"{r['unit']:>9,.0f}  {cells}  {r[f'c{n}_coal']:>7.1%}{mark}")
    print()
    if single:
        print("Percentages are shares of the total value in places where a typical dot would")
        print("touch its neighbor. In a single-class dot map some merging is expected: it is how")
        print("the densest areas read as dense. The suggestion is the smallest unit where dots")
        print(f"merge in at most {args.max_top:.0%} of the value, so they just begin to coalesce in")
        print("the densest areas (Mackay's rule of thumb).")
    else:
        print("Dot counts are estimates of what the map would actually show, after promotion;")
        print("they are typically within a few percent. A class with almost no dots is a class")
        print("you do not need. 'Still crowded' is the share of the value where even the largest")
        print("class would coalesce, so promotion cannot resolve it; the suggestion is the")
        print(f"smallest unit that keeps it at or under {args.max_top:.0%}.")


def cmd_dots(args):
    from .graduated import make_dot_map
    from .output import render_preview, write_outputs

    t = time.perf_counter()
    g = _read(args.input, args.where)
    print(f"read {len(g):,} polygons in {time.perf_counter() - t:.1f} s")
    scale = _scale(args)
    res = make_dot_map(g, _value_arg(args), args.unit, scale, k0=args.k0, id_column=args.id,
                       method=args.method, rng=args.seed,
                       group_radius=None if str(args.group_radius).lower() == "global"
                       else float(args.group_radius),
                       nudge=args.nudge, nudge_max=args.nudge_max)
    out = Path(args.out)
    paths = write_outputs(res, out, fmt=args.format, scale=scale)
    print("wrote " + ", ".join(str(p) for p in paths))
    s = res.summary
    print(f"dots by class: {s['dots_by_class']}")
    if s.get("categories"):
        print("by category: " + ", ".join(
            f"{k} {v['dots']:,} dots ({v['value_of_dots']:,.0f} of {v['value']:,.0f})"
            for k, v in s["categories"].items()))
    dr = s["drift"]
    print(f"drift: {dr['dots_outside_source']:,} dots ({dr['value_share_outside_source']:.1%} of value) "
          f"outside their source polygon; median {dr['displacement_px_median_of_outside']:.2f} px, "
          f"95th percentile {dr['displacement_px_p95_of_outside']:.2f} px beyond its edge; "
          f"{dr['dots_in_no_polygon']:,} in no polygon")
    acc = s.get("accuracy") or {}
    if acc.get("final"):
        def fmt(a):
            return (f"median {a['median_dots']:.2f}, 90th {a['p90_dots']:.2f}, "
                    f"max {a['max_dots']:.1f} dots; {a['within_one_dot']:.0%} within one dot")
        if acc.get("placement"):
            print(f"accuracy, placement only: {fmt(acc['placement'])}")
        print(f"accuracy, as mapped      : {fmt(acc['final'])}")
    if s["overlapping_dots"]:
        nonzero = {k: v for k, v in s["overlapping_dots"].items() if v}
        print(f"overlapping dots: {nonzero or 'none'}")
    if args.preview:
        outl = {"input": True, "none": False}.get(args.outlines.lower()) if args.outlines.lower() in ("input", "none") \
            else _read(args.outlines, None)
        p, (w, h) = render_preview(res.dots, g, scale, out / "preview.png", outlines=outl)
        print(f"preview: {p} ({w} x {h} CSS px, drawn at 2x)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mapstippler", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("input", help="polygon file (GeoParquet, GeoPackage, shapefile, ...)")
        p.add_argument("--value", required=True, nargs="+",
                       help="column of counts to map, or several columns for a categorised "
                            "map (e.g. --value dem rep), which colours each dot by category")
        p.add_argument("--where", help="pandas query to filter polygons, e.g. \"statefp == '41'\"")
        p.add_argument("--zoom", type=float, default=0, help="zoom level (0 = lower 48 ~1,200 px wide)")
        p.add_argument("--mpp", type=float, help="meters per CSS pixel (overrides --zoom)")
        p.add_argument("--sizes", type=float, nargs="+", default=[2.0, 4.5, 10.0], help="dot diameters in CSS px")
        p.add_argument("--radii", type=float, nargs="+",
                       help="search radii in CSS px, one per class, or a single multiplier of the "
                            "dot sizes. Two dots closer than this count as merged and get promoted. "
                            "Default: the dot diameters")
        p.add_argument("--ratio", type=float, default=10.0, help="value ratio between classes")
        p.add_argument("--k0", type=int, default=16, help="pseudo-random points per dot")

    a = sub.add_parser("advise", help="preview candidate unit values")
    common(a)
    a.add_argument("--method", choices=["graduated", "blue"], default="graduated",
                   help="map style to advise for: graduated classes, or single-class blue noise")
    a.add_argument("--units", type=float, nargs="+",
                   help="unit values to compare (default: round values 1, 2, 5 x 10^k)")
    a.add_argument("--max-top", type=float, default=0.05,
                   help="largest share of value where the largest class may coalesce (default 0.05)")
    a.set_defaults(func=cmd_advise)

    d = sub.add_parser("dots", help="make the dot map")
    common(d)
    d.add_argument("--unit", type=float, required=True, help="value of one smallest dot")
    d.add_argument("--id", help="polygon ID column (e.g. geoid)")
    d.add_argument("--method", choices=["random", "blue", "graduated"], default="graduated")
    d.add_argument("--seed", type=int, default=0)
    d.add_argument("--group-radius", default="3",
                   help="pool touching dots within this many larger-class diameters "
                        "(0 = only directly touching; 'global' = whole map, as in the 2017 paper)")
    d.add_argument("--out", default="mapstippler_output", help="output folder")
    d.add_argument("--format", choices=["gpkg", "parquet"], default="gpkg",
                   help="output format for the dots. GeoPackage by default, which every "
                        "QGIS reads; GeoParquet is smaller but needs GDAL built with Arrow")
    d.add_argument("--preview", action="store_true", help="also render preview.png at the true zoom")
    d.add_argument("--outlines", default="input",
                   help="preview outlines: 'input' (default), 'none', or a path to another polygon file")
    d.add_argument("--nudge", choices=["none", "top", "all"], default="none",
                   help="push overlapping dots apart: largest class only, or all classes")
    d.add_argument("--nudge-max", type=float, default=1.0,
                   help="farthest a dot may be nudged, in its own diameters (default 1)")
    d.set_defaults(func=cmd_dots)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
