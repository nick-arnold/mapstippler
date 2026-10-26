"""A two-colour dot map from synthetic data, where the Democratic share rises
with population density, roughly as it does in the United States.

Adapt the `share_d` line to your own data, or replace the whole first half with
a real table of votes by county.

    python examples/two_color_synthetic.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import geopandas as gpd
import numpy as np
import shapely

from mapstippler.graduated import make_dot_map
from mapstippler.output import render_preview
from mapstippler.scale import MapScale

rng = np.random.default_rng(1)
n, cell = 14, 30_000.0
polys, pop, ids = [], [], []
centre = (n - 1) / 2
for i in range(n):
    for j in range(n):
        polys.append(shapely.box(i * cell, j * cell, (i + 1) * cell, (j + 1) * cell))
        r2 = ((i - centre) ** 2 + (j - centre) ** 2) / (n / 4) ** 2
        pop.append(int(900_000 * np.exp(-r2) + rng.integers(4_000, 40_000)))
        ids.append(f"{i:02d}{j:02d}")

g = gpd.GeoDataFrame({"geoid": ids, "pop": pop}, geometry=polys, crs="ESRI:102008")

# Democratic share as a smooth function of density: about 22% in the sparsest
# places, about 82% in the densest.
density = g["pop"] / (g.geometry.area / 1e6)
share_d = 0.22 + 0.60 / (1 + np.exp(-(np.log10(density) - 1.6) * 2.2))
g["dem"] = (g["pop"] * share_d).round().astype(int)
g["rep"] = g["pop"] - g["dem"]
print(f"synthetic data: {g['dem'].sum() / g['pop'].sum():.1%} Democratic overall, "
      f"{share_d.min():.0%} to {share_d.max():.0%} by polygon")

scale = MapScale(700.0)
result = make_dot_map(g, ["dem", "rep"], unit_value=2000, scale=scale,
                      id_column="geoid", nudge="all", rng=2)

out = Path(__file__).resolve().parent / "two_color_synthetic.png"
render_preview(result.dots, g, scale, out, outlines=False,
               category_colors={"dem": "#2166ac", "rep": "#b2182b"})
print("wrote", out)
