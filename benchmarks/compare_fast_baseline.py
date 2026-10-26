"""Fast CCVT vs the verified baseline port: quality, speed, scaling.

Usage: python benchmarks/compare_fast_baseline.py
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mapstippler.ccvt import CCVTOptimizer, split_capacity  # noqa: E402
from mapstippler.ccvt_fast import FastCCVT  # noqa: E402


def spacing_stats(sites, n):
    """Nearest-neighbor distance / mean spacing, interior dots of the unit square."""
    d, _ = cKDTree(sites).query(sites, k=2)
    inner = np.all((sites > 0.05) & (sites < 0.95), axis=1)
    dn = d[inner, 1] / np.sqrt(1.0 / n)
    return dn.mean(), dn.min(), (dn < 0.5).mean()


def run_baseline(points, caps, sites0):
    t = time.perf_counter()
    opt = CCVTOptimizer(sites0, caps, points)
    opt.run()
    return opt.locations, opt.true_energy(), opt.iterations, time.perf_counter() - t


def run_fast(points, caps, sites0=None, seed=0):
    t = time.perf_counter()
    opt = FastCCVT(points, caps, initial_sites=sites0, rng=seed)
    opt.run()
    return opt.site_locations(), opt.energy(), opt.iterations, time.perf_counter() - t


def row(label, sites, energy, iters, secs, n, ref_energy=None):
    m, mn, close = spacing_stats(sites, n)
    rel = "" if ref_energy is None else f"{energy / ref_energy:10.4f}"
    print(f"  {label:<28}{m:8.3f}{mn:8.3f}{close:9.1%}{iters:7d}{secs:9.2f}s {rel}")


def main():
    rng = np.random.default_rng(11)
    n = 2000
    head = f"  {'':<28}{'meanNN':>8}{'minNN':>8}{'close':>9}{'passes':>7}{'time':>10}  energy vs baseline"

    print("A. Same starting state (greedy assignment to the same random sites)")
    for k in (2, 4, 16):
        pts = rng.uniform(0, 1, (n * k, 2))
        caps = split_capacity(len(pts), n)
        sites0 = pts[rng.choice(len(pts), n, replace=False)]
        print(f" k0 = {k}"); print(head)
        s, e_b, it, t = run_baseline(pts, caps, sites0)
        row("baseline (all pairs)", s, e_b, it, t, n, e_b)
        s, e, it, t = run_fast(pts, caps, sites0)
        row("fast (Delaunay neighbors)", s, e, it, t, n, e_b)
        s, e, it, t = run_fast(pts, caps, None)
        row("fast + recursive init", s, e, it, t, n, e_b)

    print("\nB. Scaling at k0 = 16 (fast, recursive init)")
    print(f"  {'dots':>8}{'points':>11}{'passes':>8}{'time':>10}")
    for n_dots in (1_000, 10_000, 100_000):
        pts = rng.uniform(0, 1, (n_dots * 16, 2))
        caps = split_capacity(len(pts), n_dots)
        _, _, it, t = run_fast(pts, caps)
        print(f"  {n_dots:>8,}{len(pts):>11,}{it:>8}{t:>9.2f}s")

    print("\nC. Baseline scaling for comparison (k0 = 16)")
    for n_dots in (1_000, 4_000):
        pts = rng.uniform(0, 1, (n_dots * 16, 2))
        caps = split_capacity(len(pts), n_dots)
        sites0 = pts[rng.choice(len(pts), n_dots, replace=False)]
        _, _, it, t = run_baseline(pts, caps, sites0)
        print(f"  {n_dots:>8,}{len(pts):>11,}{it:>8}{t:>9.2f}s")

    print("\nD. Large projected coordinates")
    pts = rng.uniform(0, 1, (5000 * 8, 2))
    caps = split_capacity(len(pts), 5000)
    scale, shift = 250_000.0, np.array([500_000.0, 4_800_000.0])
    a = FastCCVT(pts, caps, rng=3); a.run()

    # D1: scaled to meters and shifted to UTM-sized values. Nothing about the
    # geometry changes, so a numerically stable implementation gives the same result.
    b = FastCCVT(pts * scale + shift, caps, rng=3); b.run()
    same = (a.assignment() == b.assignment()).mean()
    diff = np.abs((b.site_locations() - shift) / scale - a.site_locations()).max()
    print(f"  scaled + shifted: same assignment {same:.1%}, max site difference {diff:.1e} "
          f"(unit-square units), passes {a.iterations} vs {b.iterations}")

    # D2: also rotated. Rotation can flip the direction of the first cut in the
    # recursive initialization, which renumbers the sites and changes the batch
    # order, so the run takes a different path. Quality should be the same.
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    c = FastCCVT((pts @ R.T) * scale + shift, caps, rng=3); c.run()
    back = ((c.site_locations() - shift) / scale) @ R           # undo the transform
    print(f"  rotated too: energy ratio {c.energy() / scale**2 / a.energy():.4f}, "
          f"mean NN {spacing_stats(a.site_locations(), 5000)[0]:.3f} vs {spacing_stats(back, 5000)[0]:.3f}")


if __name__ == "__main__":
    main()
