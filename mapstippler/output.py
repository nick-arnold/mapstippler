"""Preview images and output files."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .scale import MapScale

__all__ = ["render_preview", "write_outputs"]


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


def write_outputs(result, out_dir, fmt: str = "parquet"):
    """Write dots (GeoParquet or GeoPackage), the reconciliation table (CSV) and a
    JSON file with parameters, summary and timings."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if fmt == "parquet":
        dots_path = out / "dots.parquet"
        result.dots.to_parquet(dots_path)
    elif fmt == "gpkg":
        dots_path = out / "dots.gpkg"
        result.dots.to_file(dots_path, layer="dots", driver="GPKG")
    else:
        raise ValueError("fmt must be 'parquet' or 'gpkg'")
    rec_path = out / "reconciliation.csv"
    result.reconciliation.to_csv(rec_path, index=False)
    info_path = out / "run.json"
    info_path.write_text(json.dumps({"params": result.params, "summary": result.summary,
                                     "timings_s": {k: round(v, 2) for k, v in result.timings.items()}},
                                    indent=2, default=str))
    return dots_path, rec_path, info_path
