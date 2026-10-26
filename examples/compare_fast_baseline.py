"""Visual check: baseline vs fast CCVT at k0 = 16, plus a large map with the fast version."""
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


rng = np.random.default_rng(3)
k0 = 16
runs = []
pts = density_points(1500 * k0, rng)
for method in ("baseline", "fast"):
    t = time.perf_counter()
    r = ccvt(pts, 1500, rng=1, method=method)
    runs.append((f"{method.capitalize()} CCVT: 1,500 dots", r, time.perf_counter() - t, 0.012))
big = density_points(15000 * k0, rng)
t = time.perf_counter()
r = ccvt(big, 15000, rng=1)
runs.append(("Fast CCVT: 15,000 dots", r, time.perf_counter() - t, 0.0038))

fig, axes = plt.subplots(1, 3, figsize=(15, 5.6))
for ax, (title, r, secs, diam) in zip(axes, runs):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    fig.canvas.draw()
    pts_per_unit = ax.get_window_extent().width / fig.dpi * 72
    ax.scatter(r.sites[:, 0], r.sites[:, 1], s=(diam * pts_per_unit) ** 2, c="black", linewidths=0)
    ax.set_title(f"{title}\n{secs:.2f} s, {r.iterations} passes", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"k\u2080 = {k0}, same density surface", fontsize=12, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = Path(__file__).resolve().parent / "fast_vs_baseline.png"
fig.savefig(out, dpi=130)
print(out)
