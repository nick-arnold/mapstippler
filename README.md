# mapstippler

Pseudo-random, blue-noise and graduated dot maps for Python, rebuilding and
extending the method in:

> Arnold, N.D., Jenny, B. & White, D. (2017). Creation and evaluation of graduated
> dot maps. *International Journal of Geographical Information Science*.
> doi:10.1080/13658816.2017.1359747

Developed with the assistance of AI. The algorithms are implemented from the
cited papers rather than ported from existing code.

See `DESIGN_NOTES.md` for why this implementation differs from the paper: the
evidence behind each change, its trade-offs, and what was deliberately left alone.

## Status

**Stage 1: CCVT baseline port (done).** `mapstippler/ccvt.py` is a faithful Python port of
the discrete capacity-constrained Voronoi tessellation by Balzer, Schlömer &
Deussen (2009), following Michael Balzer's reference C++ code. It deliberately
keeps the original's all-pairs loop and strict stopping rule, and it matches the
C++ bit for bit. Use it with `method="baseline"`.

**Stage 2: fast CCVT (done, now the default).** `mapstippler/ccvt_fast.py` solves the same
problem with the accelerations from Li et al. (2010) and Ahmed & Deussen (2017):

1. sites only exchange points with their Delaunay neighbors (no all-pairs loop);
2. each pair is re-split optimally with a linear key (position along the line between
   the two sites) and a cut at the capacity;
3. recursive median-cut initialization, cutting each slice across its major axis;
4. pairs that share no site are processed together in batches with numpy, which gives
   the same result as processing them one after another;
5. pairs whose sites did not change in the last pass are skipped.

Additions over the original, in both versions: every point carries its original index,
so the result reports which input points each output dot absorbed, and energy is
computed from scratch. (The original reports a running tally that goes stale when
sites move.) The fast version also moves every site to its centroid before
optimizing; the original only centers a site after it swaps, so a site that never
swaps stays at its random starting point.

**Stage 3: dot map pipeline, version 1 (done).** Polygons in, dots out:

1. allocate points so every small dot is worth exactly the unit value (`allocate.py`);
2. place pseudo-random points in each polygon (`sampling.py`);
3. fast CCVT at capacity k0;
4. DBSCAN promotion: dots of a class with a same-class neighbor closer than one dot
   diameter are flagged; flagged dots within 3 diameters of the next class are pooled
   into groups, and each group is reduced by a CCVT constrained to that group into
   larger dots worth exactly 10 dots each. Remainders stay as they are, chosen as
   the dots farthest from the new larger dots;
5. every dot records its source polygon, how mixed its sources are, and whether it
   drifted outside its source polygon.

Map scale and dot sizes live in `scale.py`. Dots are sized in CSS pixels
(defaults 2, 3.5 and 6), the scale is meters per CSS pixel for one fixed web zoom
(zoom 0 = lower 48 about 1,200 px wide, 3,800 m/px), and the DBSCAN radius of each
class is derived from the two, so they always change together.

**v0.3.2: per-group promotion.** The 2017 method pooled every touching dot on the map
into one CCVT run, which could build a larger dot from pieces of distant towns and
place it on empty land between them. Promotion now pools only touching dots within
a linking radius (`--group-radius`, in diameters of the next class; default 3). On a
synthetic scatter of towns, whole-map pooling put 119 larger dots on empty land and
this puts 1, while dense areas are promoted almost exactly as before. Every dot is
worth exactly its class value, and each dot records `spread_px`, how far the original
points it stands for lie from it.

**v0.3.3: nudging.** `--nudge top` pushes overlapping dots of the largest class apart
(`--nudge all` does every class, with larger dots harder to move in proportion to
their area). Overlapping pairs are pushed apart in batches that share no dot, with
positions updated immediately, the same idea as Dorling cartograms and D3's
forceCollide. Each dot moves at most `--nudge-max` of its own diameters (default 1),
so where there is not enough room some overlap remains and is reported. Each dot
records `nudge_px`. `--outlines` chooses the preview's outline layer: the input
polygons, none, or another file (for example counties over tract dots).

**v0.3.4: single-class blue-noise maps.** `--method blue` makes a classic dot map with
one dot size and value, laid out as blue noise; `--method random` makes the
pseudo-random version for comparison. `advise --method blue` words its table for a
single class, where some merging in the densest areas is expected and desirable
(Mackay's rule of thumb). The advisor's default rule is now 5% (`--max-top`), and
`--units` compares specific unit values.

**v0.3.5: search radius independent of dot size.** `--radii` sets the distance at
which two dots of a class count as merged, in CSS pixels, separately from how big
they are drawn. It takes one value per class, or a single number as a multiplier of
the dot sizes. The default is still the dot diameter, as in the paper. This
separates two decisions that were tangled together: how big a dot looks, and how
close is too close. At small sizes dots read as one blob before they physically
overlap, so a radius above the diameter often matches what the eye sees, and it
lets the texture of the smallest class be tuned without changing how much gets
promoted. The advisor uses the same radii, so its table predicts the run.

**v0.3.6: ink preservation, and new default sizes.** Promoting k dots into one
larger dot removes their ink and draws one dot in its place, so the ink kept is
(d_large / d_small)^2 / k. Sizes too close together empty out the places that were
densest, leaving holes exactly where the map should be darkest; sizes too far apart
mean the larger class coalesces at the same density as the one below it, which is
the paper's argument against area-proportional dots. The two are exactly
reciprocal: ink kept x extra density absorbed = 1. `advise` now reports both per
class, and the defaults changed from 2 / 3.5 / 6 px (31% of the ink, flagged as
"holes where it promotes") to **2 / 4.5 / 10 px** (51%, "balanced"), found by
testing on national cattle data.

**v0.4.1: refactor.** `make_dot_map` was one 270-line function; it is now a short
sequence of `_step_*` functions, one per stage, with two small objects carrying
state between them (`_Input` for what was read from the layer, `_Dots` for the dots
as they are built up). Verified by running the old and new versions on the same
data: identical dot positions, columns, reconciliation, summary and params in six
configurations. Also fixed: `advise --method blue` silently dropped `--radii`.

**v0.4.2: `nudge.py` rewritten for readability.** Same algorithm, same results, but
the inner loop no longer rebinds its variables mid-function and the names say what
they hold (`direction`, `overlap`, `share_a`) instead of `u`, `d` and `tot`. The
separation step, the batching and the displacement cap are each their own small
function. Verified identical against the previous version on seven configurations
(mixed radii and weights, caps, fixed dots, coincident centres, required gaps).

**v0.5.0: categorised (two-colour) dot maps.** `--value dem rep` (or a list in
Python) colours each dot by category, in proportion, rather than colouring whole
polygons by their winner. A 60/40 county gets 60% blue dots and 40% red,
interleaved.

The category has to be chosen *per dot*, not per point: if each pseudo-random
point were given a category and each dot took the majority of its points, a 60/40
polygon would come out about 72/28 at k0 = 16, because the majority wins nearly
every dot. Instead the CCVT runs once on the total, giving one well-spaced
pattern with no possibility of colours colliding, and each polygon's dots are then
split by largest remainder and shuffled. Promotion runs separately per category,
so a large dot is worth exactly `class_ratio` smaller dots *of its own colour*.

Outputs gain a `category` column (an ordered categorical, so colours follow the
order the columns were given), per-category columns in `reconciliation.csv`, and a
`categories` block in the summary. The preview colours them with a
colour-blind-safe default palette, drawing each class in one shuffled pass so
neither colour is painted over the other.

Deliberately left for later: keeping dots inside their polygons, resolving
large-on-small overlaps (counted, not fixed), water masking, multiple zooms.


## Validation against the original

`reference/compare_reference.py` runs the original C++ (with the plain
Euclidean metric) and the Python port on identical inputs. Results are
bit-for-bit identical in site positions, with the same number of passes, for:

| case | sites | points per site | passes |
|---|---|---|---|
| uniform density, k = 2 | 300 | 2 | 4 |
| non-uniform density, k = 2 | 500 | 2 | 4 |
| non-uniform density, uneven capacities | 200 | 10–11 | 8 |
| non-uniform density, graphics regime | 128 | 64 | 16 |

To rerun: copy the four `ccvt_*.h` headers from github.com/michaelbalzer/ccvt into
`reference/`, then

```
g++ -O2 -o reference/ccvt_reference reference/ccvt_reference_driver.cpp
python reference/compare_reference.py
```

## Finding: k0 = 2 gives only partial blue noise

With two points per site, each dot sits near the midpoint of two random points and
the optimizer stops after a handful of passes. Spacing statistics for 2,000 dots on
a uniform density (nearest-neighbor distance divided by the mean spacing; interior
dots only):

| method | mean NN | min NN | NN < 0.5 spacing | passes |
|---|---|---|---|---|
| pseudo-random | 0.49 | 0.02 | 55% | – |
| CCVT, k0 = 2 | 0.61 | 0.09 | 30% | 4 |
| CCVT, k0 = 4 | 0.72 | 0.26 | 9% | 6 |
| CCVT, k0 = 8 | 0.78 | 0.43 | 1% | 8 |
| CCVT, k0 = 16 | 0.82 | 0.51 | 0% | 11 |
| CCVT, k0 = 32 | 0.84 | 0.57 | 0% | 16 |

k0 does not change the number or value of the final dots, only how many
pseudo-random points the method starts from, so raising it changes no map meaning.
It also matters for the graduated pipeline: at k0 = 2, leftover random clumps in
sparse areas can be flagged by DBSCAN and promoted, so some large dots would reflect
residual randomness rather than real density. `examples/compare_oversampling.py`
shows this visually. The fast version behaves the same way at k0 = 2 (it also stops
after about 4 passes), so this is a limit of having only two points per site, not of
the original code.

## Install

```
pip install mapstippler
```

## Command line

```
mapstippler advise polygons.parquet --value population
mapstippler dots polygons.parquet --value population --unit 10000 --id geoid --preview
```

## In Python

```python
import geopandas as gpd
from mapstippler import MapScale, make_dot_map
from mapstippler.output import write_outputs

polygons = gpd.read_parquet("tracts.parquet")      # an equal-area projection
scale = MapScale.zoom(0)                            # lower 48 about 1,200 px wide
result = make_dot_map(polygons, "population", unit_value=10_000, scale=scale,
                      id_column="geoid", nudge="all")

result.dots            # GeoDataFrame: class, value, size_px, source_id, drift
result.reconciliation  # per polygon: its value vs the value of the dots inside it
result.summary         # counts by class, promotions, drift, overlaps
write_outputs(result, "out/")
```

The input must be in a projected, ideally equal-area, coordinate system: dot
density is preserved per unit area of the flat map.

`advise` previews candidate unit values at the chosen zoom and dot sizes without
placing any dots: how many dots, and how much of the value sits in polygons where
each class would coalesce. `dots` makes the map and writes `dots.parquet`,
`reconciliation.csv` (per-polygon value vs value of the dots inside it) and
`run.json` (parameters, summary, drift, overlaps, timings). `--preview` renders
`preview.png` at the true zoom with dots at their CSS pixel sizes, at 2x resolution.

Options: `--zoom`, `--mpp` (meters per pixel), `--sizes`, `--ratio`, `--k0`,
`--method random|blue|graduated`, `--where` (pandas filter), `--format parquet|gpkg`, `--seed`,
`--radii`, `--group-radius`, `--nudge none|top|all`, `--nudge-max`, `--outlines`.
`--value` takes several column names for a categorised map. For `advise`:
`--method graduated|blue`, `--units`, `--max-top`.

## Fast version vs baseline

`benchmarks/compare_fast_baseline.py`. Both versions start from the identical state
(greedy assignment to the same random sites), 2,000 dots on a uniform density.
Energy is relative to the baseline; lower is better.

| k0 | method | mean NN | points < 0.5 spacing | time | energy |
|---|---|---|---|---|---|
| 2 | baseline | 0.61 | 32% | 0.71 s | 1.00 |
| 2 | fast | 0.64 | 28% | 0.15 s | 0.86 |
| 2 | fast + recursive init | 0.63 | 26% | 0.04 s | 0.65 |
| 16 | baseline | 0.82 | 0% | 3.89 s | 1.00 |
| 16 | fast | 0.82 | 0% | 1.09 s | 0.99 |
| 16 | fast + recursive init | 0.81 | 0% | 0.25 s | 0.98 |

The fast version reaches equal or lower energy because each pair split is optimal,
while the original's bounding-circle shortcut misses some beneficial swaps.

Scaling at k0 = 16, fast version with recursive initialization:

| dots | points | time |
|---|---|---|
| 1,000 | 16,000 | 0.13 s |
| 10,000 | 160,000 | 2.3 s |
| 100,000 | 1,600,000 | 31 s |

The baseline needs 9.5 s for 4,000 dots, and its run time grows with the square of
the number of dots.

Numerical check: the same data scaled to meters and shifted to UTM-sized coordinates
gives an identical result (every point assigned to the same dot, site positions equal
to about 1e-15). Rotating the data takes a different but equally good path (energy
within 0.1%), because rotation can flip the direction of the first cut in the
recursive initialization.

## Layout

```
DESIGN_NOTES.md                        every change from the paper, and why
mapstippler/ccvt.py                         baseline CCVT (faithful port) and the ccvt() entry point
mapstippler/ccvt_fast.py                    fast CCVT
mapstippler/allocate.py                     points per polygon (largest remainder)
mapstippler/sampling.py                     pseudo-random points in polygons
mapstippler/scale.py                        map scale, dot sizes, derived search radii
mapstippler/advisor.py                      unit value advisor
mapstippler/graduated.py                    the dot map pipeline
mapstippler/output.py                       output files and preview rendering
mapstippler/nudge.py                        pushing overlapping dots apart
mapstippler/__main__.py                     command line
tests/                                 invariants, edge cases, fast-vs-baseline quality
reference/ccvt_reference_driver.cpp    runs the original C++ on an input file
reference/compare_reference.py         port-vs-original validation
benchmarks/compare_fast_baseline.py    quality, speed, scaling, large coordinates
examples/compare_oversampling.py       random vs k0 = 2 vs k0 = 16 figure
examples/compare_fast_baseline.py      baseline vs fast figure
```
