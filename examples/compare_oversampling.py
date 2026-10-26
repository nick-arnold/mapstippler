"""Visual comparison: pseudo-random dots vs CCVT at k0 = 2 and k0 = 16, same number of dots."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from mapstippler.ccvt import ccvt

def density_points(n, rng):
    """Non-constant test density from the reference example, on the unit square."""
    out = np.empty((0, 2))
    while len(out) < n:
        xy = rng.uniform(-1, 1, (4 * n, 2))
        x, y = xy[:, 0], xy[:, 1]
        p = np.exp(-20 * x * x - 20 * y * y) + 0.2 * np.sin(np.pi * x) ** 2 * np.sin(np.pi * y) ** 2
        out = np.vstack([out, (xy[p >= rng.uniform(size=len(p))] + 1) / 2])
    return out[:n]

def touching_share(p, diameter):
    d, _ = cKDTree(p).query(p, k=2)
    return (d[:, 1] < diameter).mean()

rng = np.random.default_rng(3)
n_dots = 1500
dot_diameter = 0.012          # in map units (unit square)

panels = [("Pseudo-random", density_points(n_dots, rng))]
for k in (2, 16):
    r = ccvt(density_points(n_dots * k, rng), n_dots, rng=rng)
    panels.append((f"CCVT, k\u2080 = {k}  ({r.iterations} passes)", r.sites))

fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
for ax, (title, p) in zip(axes, panels):
    # marker size in points^2: convert map-unit diameter to points for this axes
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    fig.canvas.draw()
    bbox = ax.get_window_extent()
    pts_per_unit = bbox.width / fig.dpi * 72
    ax.scatter(p[:, 0], p[:, 1], s=(dot_diameter * pts_per_unit) ** 2, c="black", linewidths=0)
    share = touching_share(p, dot_diameter)
    ax.set_title(f"{title}\n{share:.0%} of dots touch or overlap a neighbor", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"{n_dots} dots on the same density surface", fontsize=12, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = Path(__file__).resolve().parent / "oversampling_comparison.png"
fig.savefig(out, dpi=130)
print(out)
