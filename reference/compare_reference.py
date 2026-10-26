"""Compare the Python CCVT port against Balzer's original C++ on identical inputs.

Usage: python reference/compare_reference.py
Requires the compiled driver: g++ -O2 -o reference/ccvt_reference reference/ccvt_reference_driver.cpp
"""

import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from mapstippler.ccvt import CCVTOptimizer, split_capacity  # noqa: E402

DRIVER = HERE / "ccvt_reference"


def density_points(n, rng, size=1000.0):
    """Rejection-sample the non-constant density used in the reference example."""
    out = []
    while len(out) < n:
        x, y = rng.uniform(-1, 1, 2)
        p = np.exp(-20 * x * x - 20 * y * y) + 0.2 * np.sin(np.pi * x) ** 2 * np.sin(np.pi * y) ** 2
        if p >= rng.uniform():
            out.append(((x + 1) / 2 * size, (y + 1) / 2 * size))
    return np.array(out)


def run_reference(sites, capacities, points):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(f"{len(sites)} {len(points)}\n")
        for (x, y), c in zip(sites, capacities):
            f.write(f"{x:.17g} {y:.17g} {int(c)}\n")
        for x, y in points:
            f.write(f"{x:.17g} {y:.17g}\n")
        path = f.name
    out = subprocess.run([str(DRIVER), path], capture_output=True, text=True, check=True).stdout
    energies, locs, info = [], [], {}
    for line in out.splitlines():
        parts = line.split()
        if parts[0] == "iteration":
            energies.append(float(parts[2]))
        elif parts[0] == "site":
            locs.append((float(parts[1]), float(parts[2])))
        else:
            info[parts[0]] = float(parts[1])
    return np.array(locs), energies, info


def compare(name, sites, capacities, points):
    ref_sites, ref_energy, info = run_reference(sites, capacities, points)

    t0 = time.perf_counter()
    opt = CCVTOptimizer(sites, capacities, points)
    py_energy = []
    opt.run(callback=lambda o: py_energy.append(o.reported_energy))
    py_seconds = time.perf_counter() - t0

    diff = np.abs(opt.locations - ref_sites).max()
    same_iters = len(py_energy) == len(ref_energy)
    rel_e = abs(py_energy[-1] - ref_energy[-1]) / ref_energy[-1]
    print(f"{name}")
    print(f"  sites={len(sites)} points={len(points)} capacity={capacities.min()}-{capacities.max()}")
    print(f"  passes: C++ {len(ref_energy)}, Python {len(py_energy)}  (match: {same_iters})")
    print(f"  max site position difference: {diff:.3e}")
    print(f"  final reported energy relative difference: {rel_e:.3e}")
    print(f"  time: C++ {info['seconds']:.3f}s, Python {py_seconds:.3f}s")
    return diff


def main():
    rng = np.random.default_rng(42)
    cases = []

    # Dot-map regime from the paper: reduction factor k0 = 2.
    pts = rng.uniform(0, 1000, (600, 2))
    cases.append(("uniform, k=2", pts, 300))

    pts = density_points(1000, rng)
    cases.append(("non-uniform density, k=2", pts, 500))

    # Promotion regime: capacity 10, with a remainder (capacities 10 and 11).
    pts = density_points(2003, rng)
    cases.append(("non-uniform density, k=10, uneven", pts, 200))

    # Graphics regime: many points per site.
    pts = density_points(128 * 64, rng)
    cases.append(("non-uniform density, 64 points per site", pts, 128))

    worst = 0.0
    for name, pts, n_sites in cases:
        sites = pts[rng.choice(len(pts), n_sites, replace=False)]
        caps = split_capacity(len(pts), n_sites)
        worst = max(worst, compare(name, sites, caps, pts))
    print(f"\nlargest site difference across all cases: {worst:.3e}")


if __name__ == "__main__":
    main()
