# Design notes

Why this implementation differs from Arnold, Jenny & White (2017), and what the
evidence was for each change. Read alongside the paper; section numbers refer to it.

Contents:
1. What stayed the same
2. Speed: the fast CCVT
3. k0: from 2 to 16
4. Promotion: per-group instead of global
5. Exact dot values and leftovers
6. Scale, dot size and the DBSCAN radius
7. Nudging overlapping dots
8. The unit value advisor
9. Diagnostics
10. Things deliberately not done
11. Open questions

---

## 1. What stayed the same

The method is still the paper's, and each of these was a deliberate choice to keep:

- **CCVT does the reduction.** Pseudo-random points are reduced to dots by a
  capacity-constrained Voronoi tessellation, so each dot owns the same number of
  points and therefore represents the same value, and the dots follow the input
  density (section 3.2).
- **DBSCAN finds coalescence,** with the search radius set to the dot diameter
  and the minimum cluster size set to 2, so any two dots that touch are found
  (section 3.3).
- **Classes are built by running CCVT again** on the dots that touch, at a
  capacity equal to the ratio of consecutive unit values (section 3.4).
- **Three classes with a ratio of 10** as the default (section 3.4, 4.1).
- **Dots may cross polygon boundaries.** The paper reported this as a limitation;
  see section 9 below for why it is better treated as an input choice.

---

## 2. Speed: the fast CCVT

**What changed.** The optimizer was rewritten following Li et al. (2010) and
Ahmed & Deussen (2017): sites only exchange points with their Delaunay
neighbours; each pair is re-split optimally with a linear key and a cut at the
capacity; sites start from a recursive median cut; independent pairs are
processed together in numpy.

**Why.** The original checks every pair of sites, so its run time grows with the
square of the number of dots. Ahmed & Deussen note that at 100,000 sites this
means about 10 billion neighbour queries per iteration, over 100 times the useful
work. That is what limited the 2017 tool to modest areas.

**Evidence.** `reference/compare_reference.py` runs the original C++ and the
Python port on identical inputs: identical results, bit for bit, in four regimes.
`benchmarks/compare_fast_baseline.py` then compares the fast version against that
verified baseline: from the same starting state it reaches equal or slightly
lower energy (each pair split is optimal, where the original's bounding-circle
shortcut misses some swaps), and run time grows roughly in proportion to the dot
count. 100,000 dots from 1.6 million points take about 31 s.

**Trade-off.** More code, and a second implementation to keep correct. Mitigated
by keeping the faithful port available as `method="baseline"` and testing the two
against each other.

---

## 3. k0: from 2 to 16

**What changed.** The default number of pseudo-random points per dot went from 2
(section 3.4) to 16.

**Why.** With k0 = 2, each dot sits at the midpoint of two random points, and the
optimizer converges after about 4 passes. There is too little to average over,
so clumps in the random input survive into the result.

**Evidence.** 2,000 dots on a uniform density, nearest-neighbour distance as a
share of mean spacing (`examples/compare_oversampling.py`):

| method | mean NN | dots closer than half the spacing |
|---|---|---|
| pseudo-random | 0.49 | 55% |
| CCVT, k0 = 2 | 0.61 | 30% |
| CCVT, k0 = 4 | 0.72 | 9% |
| CCVT, k0 = 8 | 0.78 | 1% |
| CCVT, k0 = 16 | 0.82 | 0% |

At k0 = 2, 30% of dots still have a neighbour closer than half the mean spacing,
which is only a little better than random. This matters beyond appearance: those
leftover clumps are what DBSCAN flags, so at k0 = 2 some dots are promoted
because the sampling was uneven, not because the area is dense.

**Why it is safe.** k0 changes neither the number of dots nor their value. It
only changes how many pseudo-random points the method starts from and what each
point is worth (d0 = d1 / k0). The map's meaning is unchanged.

**Trade-off.** More points to optimize, so slower. At national scale this is
seconds, not minutes, so the quality is worth it.

---

## 4. Promotion: per-group instead of global

**What changed.** The paper pools every dot DBSCAN flags into one CCVT run
(section 3.4). Here, flagged dots are first split into groups of dots within a
linking radius of each other (default 3 diameters of the next class up), and each
group gets its own constrained CCVT run.

**Why.** A dot always lands at the centroid of the points it owns. Under global
pooling, one larger dot can be built from touching dots in two distant towns, and
it then lands somewhere between them, on empty ground. You saw this on the first
national map: a stray large dot near the Oregon–Idaho–Nevada corner, roughly
halfway between Boise and Reno, and more in the Midwest.

**Evidence.** On a synthetic scatter of 220 towns, counting larger dots that fall
on empty land between towns:

| pooling | larger dots on empty land | max spread (px) |
|---|---|---|
| whole map (2017 method) | 119 | 230 |
| linked within 3 diameters | 1 | 7 |
| directly touching only (0) | 1 | 7 |

Dense areas promote almost identically either way; only isolated clusters differ.
`tests/test_pipeline.py::test_distant_clusters_never_merge_into_a_dot_between_them`
locks this in: two towns 200 km apart, and no dot may land in the gap.

**Why a linking radius rather than "directly touching only".** With 0, a dense
area breaks into many small chains of touching dots, each too small to promote,
so much of the crowding is never resolved (7,435 small dots still touching in the
dense-core test, against 127 at radius 3). The linking radius lets scattered
touches across one dense area pool together, while keeping separate towns apart.
3 is a middle setting; `--group-radius global` restores the paper's behaviour.

**Trade-off.** Another parameter. Small clusters (under 10 touching dots) now
stay as they are instead of being merged into a distant dot, so a few honest
overlaps remain where the paper's method would have hidden them elsewhere. That
is the right trade: an overlap in the right place beats a dot in the wrong one.

---

## 5. Exact dot values and leftovers

**What changed.** Every dot of a class is worth exactly its class value. Where a
group of touching dots does not divide evenly, the remainder stays as dots of
their own class.

**Why.** Balzer's implementation allows capacities to differ by one, and pooled
nationally that is harmless. Per group it is not: a group of 14 promoted as a
single dot would make that dot worth 1.4 times its stated value. A dot map's
contract with the reader is that every dot of a size means the same number.

**How leftovers are chosen.** The CCVT is run on the group with as many
capacity-`ratio` sites as fit; which dots are left over is decided by
concentration steps: start with the dots farthest from the group's centre,
optimize, then re-pick the leftovers as the dots farthest from their nearest new
larger dot, and repeat until the choice stops changing (usually 2 to 4 rounds).
The larger dots end up on the dense core of the group and the leftovers on its
edge, which is where a cartographer would put them.

**Trade-off.** Leftover dots in a dense area may still touch. They are counted in
the overlap diagnostics, and nudging (section 7) can resolve them.

---

## 6. Scale, dot size and the DBSCAN radius

**What changed.** The design is explicit that a graduated dot layer is made for
one fixed display scale. `MapScale` holds meters per CSS pixel and the dot
diameters in CSS pixels; every real-world distance, above all the DBSCAN radius,
is derived from those.

**Why.** "These two dots touch" only means something at a known scale. The paper
sets the search radius to the dot diameter (section 3.3), which presumes a
display size. Making that dependency explicit means the radius can never drift
out of step with what is drawn.

**Why dot sizes are fixed and the unit value moves.** Dot sizes are in CSS
pixels, so they are constant across scales. Zooming in one level halves the
meters per pixel, so each dot covers a quarter of the ground area, and the unit
value should fall by about 4 to keep the same balance of sparse and crowded
areas. Esri's dot density renderer works the same way, recalculating the dot
value against a reference scale.

**The default sizes, 2 / 3.5 / 6 px.** Older screen guidance (Jenny, Jenny &
Räber 2008) recommended a minimum of about 6 px for point symbols, but that was
based on displays under 100 ppi. Ledermann's work at TU Wien on high-density
displays found point symbols identifiable down to about 0.6 mm, close to the
print minimum. One CSS pixel is about 0.17 mm on a phone and 0.23 mm on a laptop,
so 0.6 mm is roughly 2.5 to 3.5 CSS px. The smallest dot only has to be *seen*,
not identified by shape, so 2 px is reasonable; the largest at 6 px lands at the
old recommended minimum for a recognizable symbol.

**Why the ratios matter.** With a unit value ratio of 10 between classes, a class
absorbs 10 / (diameter ratio)² more density before it coalesces. At ratio 2 that
is 2.5x per class; at 1.7, about 3.5x. Area-proportional dots (ratio 3.16) would
give exactly 1: every class would coalesce at the same density, which is why
graduated dots should not be area-proportional.

**The opposite constraint: ink.** Testing on national cattle data turned up a
failure the paper does not discuss. Promoting k small dots into one larger dot
removes their ink and draws a single dot in its place, keeping
(d_large / d_small)² / k of it. At the original 2 / 3.5 / 6 px defaults that is
31%, so the densest areas *lightened* as they promoted and the map showed holes
exactly where it should have been darkest. The two constraints are exactly
reciprocal — ink kept × extra density absorbed = 1 — so the usable window for a
value ratio of 10 is a diameter ratio of roughly 2 to 2.5. The defaults are now
2 / 4.5 / 10 px (ratio 2.25: 51% of the ink, twice the density per class), and
`advise` reports both figures per class so the trade-off is visible before
anything is drawn.

**Search radius, and why it is now separate.** The radius that decides promotion
defaults to the dot diameter, so dots merge exactly when they overlap. But at
small sizes two dots read as a single blob before they touch, and the radius also
controls *how much* gets promoted, which is otherwise tangled with how heavy the
smallest class looks. `radii_px` sets the radii independently (one per class, or
a single multiplier), and the advisor uses the same values, so the table predicts
the run. A radius below the diameter is allowed but warns, since it leaves
visibly overlapping dots unpromoted. Nudging always uses the true drawn size,
since that is about real overlap.

**Trade-off.** Each layer serves one zoom. Multi-scale output is future work.

---

## 7. Nudging overlapping dots

**What changed.** New optional step: overlapping dots are pushed apart, `top` for
the largest class only or `all` for every class, with larger dots harder to move
in proportion to their area.

**Why.** The paper lists large-on-small overlap as unresolved (section 5): the
CCVT that places a class knows nothing about the dots of other classes already on
the map. Leftovers (section 5 above) can also touch.

**How.** Overlapping pairs are pushed apart along the line between their centres,
in batches that share no dot, with positions updated immediately after each
batch. This is the same idea as Dorling cartograms and D3's forceCollide. Each
dot is capped at `--nudge-max` of its own diameters (default 1), so where there
is not enough room, some overlap remains and is reported rather than dots being
pushed far from the people they represent.

**Trade-off.** Nudging moves dots away from where the density optimization put
them, so it trades a little positional accuracy for legibility. The cap bounds
that trade, and `nudge_px` records it per dot. For single-class maps, leave it
off: merging in dense areas is how those maps convey density.

---

## 8. The unit value advisor

**What it is.** A preview, before any dots are placed, of what candidate unit
values would do at the chosen scale and dot sizes.

**Why.** The paper notes that unit value and dot size must be chosen together
(section 2, citing Mackay 1949), but leaves the choice to the cartographer with
no way to see the consequences without making the map.

**How density is judged.** Each class is assessed as if all the value in an area
were shown in that class: for value V over area A, n = V / u dots at mean spacing
s = sqrt(A / n), and the class coalesces when the diameter reaches about 0.8 s
(measured from this library's own blue-noise output). Density is measured over a
grid of cells 3 dot diameters wide, not per polygon.

**Why neighbourhoods rather than polygons.** The first version measured density
per polygon and suggested 100,000 people per dot for the country, which would
have left 2,499 counties empty. The cause: at national scale many dense polygons
are smaller than a single dot (75 counties were smaller than one large dot), so
per-polygon density says nothing about how crowded the map looks. Their dots
necessarily mix with their neighbours'. Measuring over dot-sized neighbourhoods
is both more realistic and closer to what the eye does.

**Trade-off.** It is an estimate, not a simulation: it assumes dots spread evenly
within each polygon and ignores promotion. It is meant to narrow the choice to
one or two candidates, not to replace looking at the map.

---

## 9. Diagnostics

The paper reports drift as a known limitation but never quantifies it. Every run
now measures it:

- **`inside_source` / `displacement_px`** — whether a dot ended up outside the
  polygon that contributed most of its points, and how far beyond its edge.
- **`purity`** — what share of a dot's points came from that one polygon.
- **`spread_px`** — how far the points a dot stands for lie from it. This is the
  direct test for the misplaced-dot problem of section 4: a dot built from two
  distant towns has a huge spread.
- **`nudge_px`** — how far nudging moved the dot.
- **overlap counts** for every pair of classes.
- **`reconciliation.csv`** — per polygon: its value, the value of the dots inside
  it, and the difference.

**One result worth understanding.** On the national tract map, 85% of dot value
sits outside its source polygon, which sounds alarming and is not. Half of those
dots are within 0.56 px of their source tract's edge and 95% within 4 px. At this
scale a dot stands for 16 points drawn from several tiny tracts, so it cannot be
inside all of them; what matters is that it is within a pixel of where its people
are. This is why the fix for drift is choosing enumeration units small relative
to the dots at the display scale, not constraining dots to polygons: with county
input the dots were trapped inside counties and the map looked like a choropleth
made of dots, with Phoenix smeared across all of Maricopa.

---

## 10. Things deliberately not done

- **Constraining dots to their polygons.** The group machinery in
  `ccvt_fast.py` would support it, but with appropriately small input units drift
  is sub-pixel, and with coarse units the constraint produces the choropleth look
  described above. Choosing the right input is the better answer.
- **A fourth dot class.** Extreme cores (NYC, LA, Chicago) still crowd at the
  largest class. A fourth class worth 10 million would fix it; nudging currently
  handles it well enough.
- **Water masking.** Tract polygons include coastal water, so dots can fall in
  sounds and bays. Block group 0 is water-only and would clip the coast for free;
  inland lakes would need a water layer.
- **Multiple zoom levels.** Each layer serves one zoom.

---

## 11. Open questions

Worth testing, and candidates for new work rather than reimplementation:

1. **Is 3 classes right?** The paper picked 3 without testing it. With promotion
   now cheap, 2, 3, 4 and 5 could be compared on the same data.
2. **The subitizing idea** (paper, section 5): choose unit values so most areas
   hold 4 or 5 dots per class, which readers can take in at a glance. The
   advisor's framework could search for such values.
3. **A user study without the confound.** The 2017 study compared single-class
   maps at 200 per dot against graduated maps starting at 1,000, which the paper
   notes may explain part of the advantage. This tool can now generate all three
   map types at identical unit values and dot sizes, so the comparison can be run
   cleanly.
4. **Multi-scale graduated dots.** Dots that promote as you zoom out and split as
   you zoom in. Nothing published does this with graduated classes.
5. **Does per-group promotion change what readers estimate?** The dots are in
   better places; whether that improves accuracy is an empirical question.
