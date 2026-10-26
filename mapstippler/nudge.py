"""Nudging overlapping dots apart.

The idea is the one behind Dorling cartograms and D3's forceCollide, and it is
simple: any two dots closer together than they should be are pushed apart along
the line between their centres, each moving a share of the overlap. Pushing one
pair apart can push a dot into a third, so this repeats until nothing overlaps.

Three things make it more than a loop over pairs:

* **Shares.** A pair splits the overlap in inverse proportion to `weight`, so a
  heavy dot barely moves while a light one gives way. The pipeline passes each
  dot's area as its weight, so small dots move around large ones. A dot marked
  immovable pushes others but never moves itself.
* **Batches.** Each round, the overlapping pairs are split into batches that
  share no dot, so a whole batch can be moved at once with array arithmetic.
  Positions update after each batch rather than at the end of the round, which
  spreads crowding outward much faster in a jammed cluster.
* **A cap.** Each dot may move only so far from where it started, so a dot can
  never drift away from the value it represents. Where there is not enough room
  within the cap, some overlap remains: it is reported, not hidden.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

__all__ = ["NudgeResult", "nudge"]

# Overlaps smaller than this share of the required distance count as resolved.
# Without it, floating-point noise keeps pairs "overlapping" forever.
TOL = 1e-6


@dataclass
class NudgeResult:
    xy: np.ndarray            # new positions
    moved: np.ndarray         # distance each dot moved from its start
    iterations: int
    overlapping_before: int   # dots overlapping at least one other, at the start
    overlapping_after: int    # ... at the end
    capped: int               # dots that reached their displacement cap


def _find_overlapping_pairs(positions, radius, gap):
    """Every pair of dots closer than their radii (plus `gap`) allow.

    The KD-tree first finds all pairs within the largest possible distance, then
    each pair is checked against its own required distance.
    """
    widest = 2 * radius.max() + gap
    pairs = cKDTree(positions).query_pairs(widest, output_type="ndarray")
    if pairs.shape[0] == 0:
        return pairs
    a, b = pairs[:, 0], pairs[:, 1]
    required = radius[a] + radius[b] + gap
    apart = np.hypot(positions[b, 0] - positions[a, 0], positions[b, 1] - positions[a, 1])
    return pairs[apart < required * (1 - TOL)]


def _count_overlapping_dots(positions, radius, gap) -> int:
    pairs = _find_overlapping_pairs(positions, radius, gap)
    return int(np.unique(pairs).size) if pairs.size else 0


def _independent_batches(pairs, n_dots, rng):
    """Split pairs into batches in which no dot appears twice.

    Two pairs sharing a dot cannot be moved in the same array operation, because
    each would overwrite the other's move. Each round gives every remaining pair
    a random priority and takes the pairs that win at both of their dots, which
    is Luby's algorithm for a maximal matching.
    """
    remaining = pairs
    while remaining.shape[0]:
        priority = rng.random(remaining.shape[0])
        best_at_dot = np.full(n_dots, -1.0)
        np.maximum.at(best_at_dot, remaining[:, 0], priority)
        np.maximum.at(best_at_dot, remaining[:, 1], priority)
        wins_both = ((priority == best_at_dot[remaining[:, 0]])
                     & (priority == best_at_dot[remaining[:, 1]]))
        yield remaining[wins_both]
        remaining = remaining[~wins_both]


def _apply_cap(proposed, start, cap, dots):
    """Pull any dot that would move further than its cap back onto the cap.

    `proposed` and `start` are the new and original positions of `dots`.
    """
    if cap is None:
        return proposed
    offset = proposed - start[dots]
    distance = np.hypot(offset[:, 0], offset[:, 1])
    too_far = distance > cap[dots]
    if too_far.any():
        shrink = (cap[dots][too_far] / distance[too_far])[:, None]
        proposed[too_far] = start[dots][too_far] + offset[too_far] * shrink
    return proposed


def _separate_batch(positions, start, batch, radius, gap, ease, cap, rng):
    """Push one batch of pairs apart. Returns the largest distance any dot moved.

    The batch shares no dots, so every pair is handled at once.
    """
    a, b = batch[:, 0], batch[:, 1]
    offset = positions[b] - positions[a]
    apart = np.hypot(offset[:, 0], offset[:, 1])
    required = radius[a] + radius[b] + gap

    still_overlapping = apart < required * (1 - TOL)
    if not still_overlapping.any():
        return 0.0
    a, b = a[still_overlapping], b[still_overlapping]
    offset, apart = offset[still_overlapping], apart[still_overlapping]
    required = required[still_overlapping]

    # Dots at exactly the same spot have no direction to separate along, so pick one.
    coincident = apart == 0
    if coincident.any():
        angle = rng.uniform(0, 2 * np.pi, int(coincident.sum()))
        offset[coincident] = np.column_stack([np.cos(angle), np.sin(angle)])
        apart = np.where(coincident, 1.0, apart)

    direction = offset / apart[:, None]              # unit vector from a towards b
    overlap = (required - apart)[:, None]            # how far apart they still need to get
    # Who gives way: each side moves in proportion to how readily it gives.
    # Both shares are computed from `give` rather than one as 1 - the other, so
    # a pair splits the overlap the same way whichever order it arrives in.
    give = ease[a] + ease[b]
    share_a = (ease[a] / give)[:, None]
    share_b = (ease[b] / give)[:, None]

    moved_a = _apply_cap(positions[a] - overlap * direction * share_a, start, cap, a)
    moved_b = _apply_cap(positions[b] + overlap * direction * share_b, start, cap, b)
    largest_step = max(np.abs(moved_a - positions[a]).max(), np.abs(moved_b - positions[b]).max())
    positions[a], positions[b] = moved_a, moved_b
    return float(largest_step)



def nudge(xy, radius, movable=None, weight=None, gap: float = 0.0, max_disp=None,
          iterations: int = 1000, rng=0) -> NudgeResult:
    """Push overlapping circles apart, moving them as little as possible.

    xy        (n, 2) circle centres
    radius    (n,) radii, in the same units as xy
    movable   (n,) bool; an immovable circle pushes others but never moves itself
    weight    (n,) how hard a circle is to move; a pair splits the overlap in
              inverse proportion, so heavier circles give way less
    gap       extra empty space required between circles
    max_disp  (n,) or one number: how far a circle may move from where it started
    iterations  maximum rounds before giving up
    """
    start = np.asarray(xy, dtype=np.float64)
    positions = start.copy()
    n_dots = start.shape[0]
    radius = np.broadcast_to(np.asarray(radius, dtype=np.float64), (n_dots,))
    movable = np.ones(n_dots, bool) if movable is None else np.asarray(movable, bool)
    weight = np.ones(n_dots) if weight is None else np.asarray(weight, dtype=np.float64)
    cap = (None if max_disp is None
           else np.broadcast_to(np.asarray(max_disp, dtype=np.float64), (n_dots,)))
    rng = np.random.default_rng(rng)

    if n_dots < 2:
        return NudgeResult(positions, np.zeros(n_dots), 0, 0, 0, 0)

    # How readily each dot gives way: the inverse of its weight, and zero if it
    # cannot move at all. A pair with no give on either side is left alone.
    ease = np.where(movable, 1.0 / weight, 0.0)
    before = _count_overlapping_dots(positions, radius, gap)

    rounds = 0
    for rounds in range(1, iterations + 1):
        pairs = _find_overlapping_pairs(positions, radius, gap)
        if pairs.size:
            pairs = pairs[(ease[pairs[:, 0]] + ease[pairs[:, 1]]) > 0]
        if pairs.shape[0] == 0:
            rounds -= 1                     # this round found nothing to do
            break

        largest_step = 0.0
        for batch in _independent_batches(pairs, n_dots, rng):
            largest_step = max(largest_step,
                               _separate_batch(positions, start, batch, radius, gap,
                                               ease, cap, rng))
        if largest_step < TOL * 1e-3 * max(radius.max(), 1e-300):
            break                           # everything that can move is at its cap

    moved = np.hypot(*(positions - start).T)
    return NudgeResult(
        xy=positions,
        moved=moved,
        iterations=rounds,
        overlapping_before=before,
        overlapping_after=_count_overlapping_dots(positions, radius, gap),
        capped=0 if cap is None else int(np.count_nonzero(moved >= cap * (1 - 1e-9))),
    )
