# Parameter guide

Every option, what it does, and when to change it. Most have defaults that work;
the ones you will actually set are marked **you set this**.

Two commands share most parameters:

```
mapstippler advise <input> --value <column>            # preview unit values
mapstippler dots   <input> --value <column> --unit N   # make the map
```

---

## The input

### `input` — **you set this**
A polygon file: GeoParquet, GeoPackage, shapefile, GeoJSON, or anything else GDAL
reads.

It must be in a **projected, equal-area** coordinate system. The tool refuses
geographic (latitude/longitude) coordinates and warns if the projection is not
equal-area. Dot density is preserved per unit area of the flat map, so a
projection that inflates area (Web Mercator, for example) makes those places look
sparser than they are. For the United States, ESRI:102008 or EPSG:5070 work well.

**Choose units small relative to your dots.** This matters more than any other
decision. If the polygons are large compared with the spacing between dots, the
dots spread evenly across each polygon and the map reads as a choropleth made of
dots, with visible boundaries and whole cities smeared across rural counties. At a
national scale, census tracts work; counties do not. The sign that units are too
coarse is that you can see the polygon shapes in the finished map.

### `--value` — **you set this**
The column of counts to map. Must be counts of things (people, cattle, votes), not
rates, percentages or averages.

Give two or more column names for a **categorised map**, which colours each dot by
category in proportion:

```
--value dem rep
```

Each polygon's split is matched by largest remainder, and each category promotes
into larger dots separately, so a large dot is worth exactly `--ratio` smaller dots
of its own colour. Colours follow the order you list the columns.

### `--where`
A filter, written as a pandas query, applied before anything else:

```
--where "statefp == '41'"
```

Useful for making a state map from a national file. Note that filtering changes
the density picture, so re-run `advise` after filtering.

### `--id`
The column holding each polygon's identifier (`geoid`, say). It becomes the
`source_id` on every dot and the `id` in `reconciliation.csv`. Without it the
row number is used. Set it if you want to join the dots back to your data later,
which is how you colour dots by an attribute you did not map.

---

## Scale: what the map is made for

A graduated dot layer is designed for **one fixed display scale**. The dots are
placed so they just avoid each other at that scale; zoom in or out far and the
balance is gone.

### `--mpp` — **you set this** (or `--zoom`)
Ground meters covered by one CSS pixel. This is the whole definition of scale.

| view | roughly | `--mpp` |
|---|---|---|
| lower 48, ~1,200 px wide | 1:14,000,000 | 3800 |
| a large state | 1:1,800,000 | 500 |
| a metro area | 1:210,000 | 60 |

To work it out for your own map: take the width of your area in meters and divide
by how many pixels wide you want it to be drawn.

### `--zoom`
An alternative to `--mpp`, on a ladder where zoom 0 is the lower 48 at about 1,200
px wide (3,800 m/px) and each level halves the meters per pixel. `--mpp` overrides
it. Only useful if you are working in that specific ladder.

### `--sizes`
Dot diameters in CSS pixels, smallest to largest. Default `2 4.5 10`.

**The number of values sets the number of classes.** `--sizes 2 4.5` gives a
two-class map; `--sizes 2` with `--method blue` gives one size.

The smallest dot only has to be *seen*, not identified by shape, so 2 px is
reasonable — about 0.35 mm on a phone, 0.45 mm on a laptop. Below about 1.5 px
dots start to disappear on standard-resolution screens.

**The ratio between sizes matters more than the sizes themselves.** See `--ratio`
below: the two have to be chosen together.

### `--radii`
The distance at which two dots of a class count as merged, in CSS pixels, one per
class. Defaults to the dot diameters, which means dots merge exactly when they
physically overlap (as in the 2017 paper).

A single number acts as a multiplier: `--radii 1.5` means one and a half diameters
for every class.

**Set this when the texture you want and the promotion you want disagree.** The
radius controls how much gets promoted; the size controls how heavy the smallest
class looks. Tied together, there may be no setting that gives both — small dots
with enough promotion, or enough promotion with dots too heavy. Separating them
solves it. At small sizes a radius above the diameter is also closer to what the
eye does, since two 2 px dots read as one blob before they actually touch.

A radius *below* the diameter is allowed but warns: it leaves visibly overlapping
dots unpromoted.

---

## What a dot is worth

### `--unit` — **you set this** (`dots` only)
What one dot of the smallest class represents. The single most consequential
choice you make.

Too low and the map fills in to solid black in dense areas; too high and sparse
areas go empty and the pattern disappears. Use `advise` to see the trade-off
before committing. There is no universally right answer: it depends on your data's
density range, the scale, and the dot sizes.

### `--ratio`
The value multiplier between consecutive classes. Default 10, so with
`--unit 1000` the classes are worth 1,000 / 10,000 / 100,000.

Must be a whole number, because promotion needs exactly that many dots to make one
larger dot.

**Match it to your data's density range.** This is the least obvious parameter and
one of the most useful:

* **Wide range** (population: metro tracts are thousands of times denser than
  rural ones) wants a **high ratio**, 10 or more. Each class then absorbs a lot
  more density than the one below it.
* **Narrow range** (cattle: the densest counties are perhaps 100× the sparsest)
  wants a **low ratio**, around 5. With a high ratio the upper classes are barely
  used, and you effectively have a one-class map with a few outliers.

A sign the ratio is too high for your data: the largest class has almost no dots.
A sign it is too low: the densest areas still coalesce at every class.

### `--sizes` and `--ratio` together: the ink rule

Promoting k small dots into one larger dot removes their ink and draws a single
dot in its place. The ink kept is:

```
(larger diameter / smaller diameter)² / ratio
```

Two constraints pull in opposite directions, and they are exactly reciprocal
(ink kept × extra density absorbed = 1):

* **Sizes too close together** → ink drops at promotion, and you get *holes*
  exactly where the map should be darkest. The old defaults (2 / 3.5 / 6 px at
  ratio 10) kept only 31% of the ink and did this visibly.
* **Sizes too far apart** → the larger class coalesces at the same density as the
  one below it, so the classes add nothing. Area-proportional dots (diameter ratio
  √10 ≈ 3.16 at ratio 10) are exactly this case, which is why graduated dots
  should not be area-proportional.

For a ratio of 10, a diameter ratio near **2.25** balances them — about half the
ink, twice the density per class. That is where the `2 4.5 10` default comes from.
For a ratio of 5, the balance point is around **1.6**, so `2 3.2 5.1`.

`advise` prints this for your settings:

```
promotion into class 2: 2.0→4.5 px keeps 51% of the ink, absorbs 2.0x the density (balanced)
```

Aim for "balanced". "Holes where it promotes" and "absorbs little extra density"
are the two failure modes above.

### `--k0`
Pseudo-random points generated per final dot. Default 16.

**Changes nothing about what the map means** — not the number of dots, not their
values — only how good the blue-noise pattern is. The 2017 paper used 2, which
leaves about 30% of dots with a neighbour closer than half the mean spacing,
little better than random. At 16 that is 0%.

Lower it only if you are impatient on a very large map; raise it above 16 for
marginal gains at real cost.

---

## Which kind of dot map

### `--method`
* `graduated` (default) — blue-noise dots promoted into larger classes where they
  would coalesce. The method of the paper.
* `blue` — a single class, blue-noise placed. The classic dot map, done well. Use
  `--sizes` with one value.
* `random` — pseudo-random placement, no optimisation. The conventional automated
  dot map, included for comparison: it shows the clumping that blue noise removes.

Running all three on the same data at the same unit value and dot size is the
clearest demonstration of what the method does.

### `--group-radius`
How merged dots are pooled before being promoted, in diameters of the *larger*
class. Default 3.

Dots that merge are pooled with other merged dots within this distance, and each
pool is promoted on its own. This exists because pooling every merged dot on the
map together — as the 2017 paper did — can build one larger dot out of touching
dots in two distant towns, and place it on empty ground between them.

* **0** pools only dots that touch directly. Dense areas then break into many small
  chains, each too small to promote, and much crowding goes unresolved.
* **3** (default) lets scattered merges across one dense area pool together while
  keeping separate towns apart.
* **`global`** reproduces the paper's behaviour, for comparison.

### `--nudge`
Push overlapping dots apart after placement. Default `none`.

* `none` — leave them. **Right for single-class maps**, where dots merging in dense
  areas is how the map conveys density.
* `top` — only the largest class.
* `all` — every class, with larger dots harder to move in proportion to their area,
  so small dots give way around large ones.

`all` is usually what you want for a graduated map: it resolves the large-on-small
overlaps that the 2017 paper left unsolved.

### `--nudge-max`
How far a dot may be pushed, in its own diameters. Default 1.

The cap exists so a dot can never drift far from the value it represents. Where
there is not enough room within the cap, some overlap remains and is reported
rather than hidden. Raise it (2, say) if dense cores stay crowded and you are
willing to trade a little positional accuracy; lower it (0.5) to keep dots closer
to where the density put them.

### `--seed`
Random seed. Default 0. The same seed on the same data gives the identical map,
on any machine. Change it if you want a different arrangement of the same
density — useful when a dot happens to land somewhere awkward.

---

## Output

### `--out`
Output folder. Default `mapstippler_output`. Contains:

* `dots.parquet` (or `.gpkg`) — the dot layer
* `reconciliation.csv` — per polygon: its value against the value of the dots
  inside it
* `run.json` — every parameter, the summary, timings

### `--format`
`parquet` (default) or `gpkg`. GeoPackage if you want to double-click it into
QGIS; GeoParquet is smaller and faster for large outputs.

### `--preview`
Render `preview.png` at the true scale, with dots at their real pixel sizes and
the image sized as the map would be. The fastest way to judge settings — much
faster than loading the layer and styling it.

### `--outlines`
What the preview draws underneath: `input` (default), `none`, or a path to another
polygon file. Use the third when the input units are too fine to outline — county
outlines over tract dots, for example.

---

## `advise` only

### `--units`
Specific unit values to compare, instead of the automatic round numbers:

```
--units 2000 2500 3000 4000
```

### `--max-top`
The largest share of the total value that may remain crowded — that is, where even
the biggest class would coalesce, so promotion cannot resolve it. Default 0.05.

Only affects which row is marked as suggested. Raise it if your data has a few
extreme cores (a national population map always will) and the suggestion is coming
back far too coarse.

### Reading the table

```
     unit     class 1     class 2     class 3     still
                 dots        dots        dots   crowded
      500      19,627       1,898         585    26.3%
    1,000      19,059       1,728         122     0.0%  <- suggested
    2,000      17,929         627           1     0.0%
    5,000       9,679           3           0     0.0%
```

The dot counts are **estimates of what the map would actually show**, after
promotion, not the raw value divided by the unit. They come from walking the
promotion chain over the measured density, and run within a few percent of a real
run.

They answer the question that matters: **is each class earning its place?** At
5,000 above, classes 2 and 3 are empty, so that is really a one-class map and
`--sizes` should have one value. At 1,000, class 3 has 122 dots — present but
marginal. At 500 it has 585, doing real work.

"Still crowded" is the share of the value where even the largest class coalesces.
Some is unavoidable on a national population map, where a few urban cores are
denser than any dot size can resolve.

### `--method` (advise)
`graduated` (default) or `blue`. The single-class version words its table
differently: there, some merging in the densest areas is expected and desirable,
following Mackay's rule of thumb that dots should just begin to coalesce where
density is highest.

---

## Reading the output

Columns on each dot:

| column | meaning |
|---|---|
| `class` | 1 = smallest |
| `value` | what this dot represents; always exact for its class |
| `size_px` | diameter in CSS pixels |
| `source_id` | the polygon that contributed most of its points |
| `category` | for categorised maps |
| `purity` | share of its points from that one polygon; 1.0 = all |
| `inside_source` | whether it ended up inside its source polygon |
| `displacement_px` | how far beyond that polygon's edge, if outside |
| `spread_px` | how far the points it stands for lie from it |
| `nudge_px` | how far nudging moved it |

**`spread_px` is the one to check.** A dot standing for points scattered far away
is a dot in the wrong place. Large values mean something has gone wrong with
pooling.

**Do not be alarmed by `inside_source`.** On a national tract map, 85% of dot value
can sit outside its source polygon, which sounds terrible and is not: half of those
dots are within half a pixel of the polygon's edge. A dot built from 16 points
drawn from several tiny tracts cannot be inside all of them. Look at
`displacement_px` instead, which is in pixels and tells you whether it is visible.

### Styling in QGIS

Set symbol size in **map units**, equal to `size_px × m_per_px`. With
`--sizes 2 4.5 10` at `--mpp 3800`, that is 7,600 / 17,100 / 38,000 map units. The
dots then draw at the intended size whatever the zoom, and the map looks right at
the scale it was made for.
