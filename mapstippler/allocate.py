"""Allocating dots and pseudo-random points to polygons.

The CCVT reduces the pseudo-random points by the factor k0, so the number of
points must be exactly k0 times the number of final dots for every dot to
represent exactly the unit value. This module makes that allocation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Allocation", "largest_remainder", "allocate_points"]


def largest_remainder(values, total: int) -> np.ndarray:
    """Split the integer `total` in proportion to `values` (Hamilton's method).

    Each share is rounded down, and the leftover units go to the largest
    fractional parts (ties go to the earlier entry). The result always sums
    exactly to `total`, and no share differs from its exact quota by 1 or more.
    """
    values = np.asarray(values, dtype=np.float64)
    if np.any(values < 0):
        raise ValueError("values must be non-negative")
    s = values.sum()
    if s == 0:
        if total != 0:
            raise ValueError("cannot allocate a positive total to all-zero values")
        return np.zeros(values.shape, dtype=np.int64)
    quota = values * (total / s)
    base = np.floor(quota).astype(np.int64)
    leftover = int(total - base.sum())
    if leftover:
        frac = quota - base
        order = np.lexsort((np.arange(values.size), -frac))   # largest fraction first
        base[order[:leftover]] += 1
    return base


@dataclass
class Allocation:
    unit_value: float          # value represented by one final dot
    k0: int                    # pseudo-random points per final dot
    n_dots: int                # total final dots
    points: np.ndarray         # pseudo-random points per polygon
    represented: np.ndarray    # value those points stand for, per polygon
    error: np.ndarray          # represented - actual, per polygon

    @property
    def n_points(self) -> int:
        return int(self.points.sum())

    @property
    def point_value(self) -> float:
        return self.unit_value / self.k0


def allocate_points(values, unit_value: float, k0: int = 16) -> Allocation:
    """Decide how many pseudo-random points go in each polygon.

    The number of final dots is round(total / unit_value). The points total is
    exactly k0 times that, split across polygons in proportion to their values,
    so the CCVT gives every final dot exactly k0 points, and each dot represents
    exactly `unit_value`. Rounding happens at the point level, which is k0 times
    finer than rounding whole dots.
    """
    values = np.asarray(values, dtype=np.float64)
    if unit_value <= 0:
        raise ValueError("unit_value must be positive")
    if k0 < 1:
        raise ValueError("k0 must be at least 1")
    n_dots = int(round(values.sum() / unit_value))
    points = largest_remainder(values, n_dots * k0)
    represented = points * (unit_value / k0)
    return Allocation(unit_value, k0, n_dots, points, represented, represented - values)
