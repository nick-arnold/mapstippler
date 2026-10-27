"""Preview images and output files."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .scale import MapScale

__all__ = ["render_preview", "write_outputs", "style_qml"]


DEFAULT_CATEGORY_COLORS = ("#2166ac", "#b2182b", "#1b7837", "#762a83",
                           "#e08214", "#01665e", "#8c510a", "#4d4d4d")


def render_preview(dots, polygons, scale: MapScale, path, device_pixel_ratio: float = 2.0,
                   outline_color="#c8c8c8", dot_color="#1a1a1a", outlines=True,
                   category_colors=None):
    """Draw the dots exactly as a web map at this zoom would: the image is
    (extent / m_per_px) CSS pixels wide, rendered at `device_pixel_ratio`,
    and every dot is drawn at its size_px in CSS pixels.

    outlines: True draws the input polygons, False draws none, and a
    GeoDataFrame draws that layer instead (e.g. counties over tract dots).

    category_colors: for a categorised map, a {category: colour} mapping, or a
    list of colours in category order. Defaults to a colour-blind-safe set.
    Categories are drawn in an interleaved order so neither colour sits on top."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x0, y0, x1, y1 = polygons.total_bounds
    pad = 10 * scale.m_per_px
    x0, y0, x1, y1 = x0 - pad, y0 - pad, x1 + pad, y1 + pad
    w_css = (x1 - x0) / scale.m_per_px
    h_css = (y1 - y0) / scale.m_per_px
    dpi = 96.0 * device_pixel_ratio                     # 1 CSS px = 1/96 inch
    fig = plt.figure(figsize=(w_css / 96.0, h_css / 96.0), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.axis("off")
    layer = polygons if outlines is True else (outlines if outlines is not False else None)
    if layer is not None:
        if layer.crs is not None and polygons.crs is not None and layer.crs != polygons.crs:
            layer = layer.to_crs(polygons.crs)
        layer.boundary.plot(ax=ax, color=outline_color, linewidth=0.5 * 0.75)    # 0.5 CSS px
    xy = np.column_stack([dots.geometry.x, dots.geometry.y])
    has_cats = "category" in dots.columns
    if has_cats:
        cats = (list(dots["category"].cat.categories)
                if hasattr(dots["category"], "cat")
                else list(dict.fromkeys(dots["category"])))
        if category_colors is None:
            colors = {c: DEFAULT_CATEGORY_COLORS[i % len(DEFAULT_CATEGORY_COLORS)]
                      for i, c in enumerate(cats)}
        elif isinstance(category_colors, dict):
            colors = category_colors
        else:
            colors = {c: category_colors[i % len(category_colors)] for i, c in enumerate(cats)}
    for cls in sorted(dots["class"].unique()):
        in_class = (dots["class"] == cls).to_numpy()
        size_pt = float(dots.loc[in_class, "size_px"].iloc[0]) * 0.75            # 1 CSS px = 0.75 pt
        if not has_cats:
            ax.scatter(xy[in_class, 0], xy[in_class, 1], s=size_pt ** 2, c=dot_color,
                       linewidths=0, zorder=2 + cls)
            continue
        # draw the categories of a class in one pass, shuffled, so neither
        # colour is systematically painted over the other
        sel = np.flatnonzero(in_class)
        order = np.random.default_rng(0).permutation(sel.size)
        sel = sel[order]
        ax.scatter(xy[sel, 0], xy[sel, 1], s=size_pt ** 2,
                   c=[colors[c] for c in dots["category"].to_numpy()[sel]],
                   linewidths=0, zorder=2 + cls)
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)
    return Path(path), (round(w_css), round(h_css))


_MARKER_OPTIONS = (
    ("angle", "0"), ("cap_style", "square"), ("horizontal_anchor_point", "1"),
    ("joinstyle", "bevel"), ("name", "circle"), ("offset", "0,0"),
    ("offset_map_unit_scale", "3x:0,0,0,0,0,0"), ("offset_unit", "MapUnit"),
    ("outline_style", "no"), ("outline_width", "0"), ("outline_width_unit", "MapUnit"),
    ("scale_method", "diameter"), ("size_map_unit_scale", "3x:0,0,0,0,0,0"),
    ("size_unit", "MapUnit"), ("vertical_anchor_point", "1"),
)


def _marker(size_map_units: float, color: str) -> str:
    opts = dict(_MARKER_OPTIONS)
    opts["color"] = color
    opts["size"] = f"{size_map_units:g}"
    body = "".join(f'<Option name="{k}" type="QString" value="{v}"/>'
                   for k, v in sorted(opts.items()))
    return ('<layer class="SimpleMarker" enabled="1" locked="0" pass="0">'
            f'<Option type="Map">{body}</Option></layer>')


def _rgba(color: str) -> str:
    """'#2166ac' -> '33,102,172,255', which is how QML writes colours."""
    c = color.lstrip("#")
    return ",".join(str(int(c[i:i + 2], 16)) for i in (0, 2, 4)) + ",255"


def style_qml(result, scale: MapScale, dot_color: str = "#1a1a1a",
              category_colors=None) -> str:
    """A QGIS style for the dot layer, as QML.

    Dot sizes are in **map units**, set to the size the dots were designed for
    (size in CSS pixels x meters per pixel). The layer then draws at the right
    size for the scale it was made for, whatever the zoom.

    One symbol per class, or per class and category together when the map is
    categorised.
    """
    dots = result.dots
    classes = sorted(int(c) for c in dots["class"].unique())
    unit = float(result.params.get("unit_value", 1))
    ratio = float(result.params.get("class_ratio", 10))
    has_cats = "category" in dots.columns
    cats = ([str(c) for c in (dots["category"].cat.categories
                              if hasattr(dots["category"], "cat")
                              else dict.fromkeys(dots["category"]))]
            if has_cats else [None])
    if category_colors is None:
        category_colors = {c: DEFAULT_CATEGORY_COLORS[i % len(DEFAULT_CATEGORY_COLORS)]
                           for i, c in enumerate(cats) if c is not None}

    rules, symbols = [], []
    for cls in classes:
        size = scale.sizes_px[min(cls - 1, scale.n_classes - 1)] * scale.m_per_px
        for cat in cats:
            i = len(symbols)
            colour = _rgba(dot_color if cat is None else category_colors.get(cat, dot_color))
            value = unit * ratio ** (cls - 1)
            label = f"{value:,.0f}" + (f" ({cat})" if cat else "")
            # field names are quoted in QGIS expressions; the quotes have to be
            # XML-escaped because the expression lives in an attribute
            filt = "&quot;class&quot; = %d" % cls
            if cat:
                filt += " AND &quot;category&quot; = '%s'" % cat
            rules.append(f'<rule key="{{{i:08x}-0000-0000-0000-000000000000}}" '
                         f'filter="{filt}" label="{label}" symbol="{i}"/>')
            symbols.append(f'<symbol type="marker" name="{i}" alpha="1" '
                           f'clip_to_extent="1" force_rhr="0" frame_rate="10">'
                           f'{_marker(size, colour)}</symbol>')

    return ('<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>'
            '<qgis version="3.34.0" styleCategories="Symbology">'
            '<renderer-v2 type="RuleRenderer" forceraster="0" symbollevels="0" '
            'enableorderby="0" referencescale="-1">'
            f'<rules key="{{00000000-0000-0000-0000-000000000001}}">{"".join(rules)}</rules>'
            f'<symbols>{"".join(symbols)}</symbols>'
            '</renderer-v2></qgis>')


def _embed_style(gpkg_path, layer: str, qml: str, style_name: str = "mapstippler") -> None:
    """Put the style in the GeoPackage's layer_styles table, as QGIS does, so the
    layer opens already symbolised."""
    import sqlite3

    con = sqlite3.connect(gpkg_path)
    try:
        con.execute("""CREATE TABLE IF NOT EXISTS layer_styles(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            f_table_catalog TEXT, f_table_schema TEXT, f_table_name TEXT,
            f_geometry_column TEXT, styleName TEXT, styleQML TEXT, styleSLD TEXT,
            useAsDefault BOOLEAN, description TEXT, owner TEXT, ui TEXT,
            update_time DATETIME DEFAULT (strftime('%Y-%m-%dT%H:%M:%fZ','now')))""")
        con.execute("DELETE FROM layer_styles WHERE f_table_name = ? AND styleName = ?",
                    (layer, style_name))
        con.execute("""INSERT INTO layer_styles
            (f_table_catalog, f_table_schema, f_table_name, f_geometry_column,
             styleName, styleQML, styleSLD, useAsDefault, description, owner, ui)
            VALUES ('', '', ?, 'geom', ?, ?, NULL, 1, ?, '', NULL)""",
                    (layer, style_name, qml, "Dot sizes in map units, for the scale "
                                             "the map was made for"))
        con.commit()
    finally:
        con.close()


def write_outputs(result, out_dir, fmt: str = "gpkg", scale: MapScale | None = None):
    """Write the dots, the reconciliation table (CSV) and a JSON file with
    parameters, summary and timings.

    With `scale`, a GeoPackage also carries a QGIS style, so the layer opens
    already symbolised with the dots at their intended size. The same style is
    written alongside as `dots.qml`.

    GeoPackage by default: every QGIS reads it. GeoParquet is smaller and faster
    for large outputs, but QGIS only reads it when its GDAL was built with Arrow
    support, which many builds leave out.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    dots = result.dots
    if fmt == "gpkg":
        dots_path = out / "dots.gpkg"
        # GeoPackage has no categorical type; write the labels as plain text
        cats = [c for c in dots.columns if str(dots[c].dtype) == "category"]
        if cats:
            dots = dots.copy()
            for c in cats:
                dots[c] = dots[c].astype(str)
        dots.to_file(dots_path, layer="dots", driver="GPKG")
        if scale is not None:
            qml = style_qml(result, scale)
            _embed_style(dots_path, "dots", qml)
            (out / "dots.qml").write_text(qml)
    elif fmt == "parquet":
        dots_path = out / "dots.parquet"
        dots.to_parquet(dots_path)
    else:
        raise ValueError("fmt must be 'gpkg' or 'parquet'")
    rec_path = out / "reconciliation.csv"
    result.reconciliation.to_csv(rec_path, index=False)
    info_path = out / "run.json"
    info_path.write_text(json.dumps({"params": result.params, "summary": result.summary,
                                     "timings_s": {k: round(v, 2) for k, v in result.timings.items()}},
                                    indent=2, default=str))
    return dots_path, rec_path, info_path
