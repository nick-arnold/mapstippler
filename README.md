# mapstippler

Dot maps for Python: pseudo-random, blue-noise, and graduated.

Implements and extends the method in:

> Arnold, N.D., Jenny, B. & White, D. (2017). Creation and evaluation of graduated
> dot maps. *International Journal of Geographical Information Science*.
> doi:10.1080/13658816.2017.1359747

---

## How it works

1. Place pseudo-random dots in each polygon, several for every final dot.
2. Reduce them with a capacity-constrained Voronoi tessellation (CCVT). Every dot
   takes the same number of points, so every dot is worth the same and the
   pattern follows the input density without clumping.
3. Find the dots that overlap at the map's scale.
4. Reduce those with CCVT again, ten to one, into fewer dots of greater value.
5. Repeat 3 and 4 for each larger dot class.

Optionally, push apart any dots still overlapping across classes.

---

## Install

```
pip install mapstippler
```

Python 3.9 or newer.

---

## Use

Input: a polygon layer in a projected, equal-area coordinate system, with a column
of counts (not rates or percentages).

### Choose a unit value

`advise` estimates what each candidate would produce, without placing dots:

```
mapstippler advise tracts.parquet --value population --mpp 500
```

```
     unit     class 1     class 2     class 3     still
                 dots        dots        dots   crowded
      500       3,052         185          36    21.5%
    1,000       1,986         111          12     0.0%  <- suggested
    2,000       1,142          34           0     0.0%
```

`--mpp` is the scale: ground meters per screen pixel. 3,800 fits the lower 48 in
about 1,200 pixels; 500 suits a state; 60 suits a city.

### Make the map

```
mapstippler dots tracts.parquet --value population --unit 1000 --mpp 500 \
    --id geoid --nudge all --out mymap --preview
```

Writes `dots.gpkg`, a per-polygon reconciliation table, a record of the
parameters, and a PNG at the true scale. `--format parquet` writes GeoParquet
instead, which is smaller but which QGIS only reads when its GDAL was built with
Arrow support.

### Other methods

```
mapstippler dots tracts.parquet --value population --unit 1000 --mpp 500 --method blue
mapstippler dots tracts.parquet --value population --unit 1000 --mpp 500 --method random
```

`blue` is one dot size, blue-noise placed. `random` is pseudo-random placement
with no optimisation.

### Two colours

Several columns colour each dot by category, matching each polygon's split:

```
mapstippler dots counties.parquet --value dem rep --unit 10000 --mpp 3800 --preview
```

### One state

```
mapstippler dots tracts.parquet --value population --unit 1000 --mpp 500 \
    --where "statefp == '41'" --out oregon --preview
```

### Python

```python
import geopandas as gpd
from mapstippler import MapScale, make_dot_map
from mapstippler.output import write_outputs

polygons = gpd.read_parquet("tracts.parquet")
result = make_dot_map(polygons, "population", unit_value=1000,
                      scale=MapScale(500.0), id_column="geoid", nudge="all")

result.dots            # class, value, source polygon, diagnostics
result.reconciliation  # per polygon: its value against the dots inside it
write_outputs(result, "mymap/")
```

---

## Documentation

* [PARAMETERS.md](PARAMETERS.md) — every option and when to change it.
* [DESIGN_NOTES.md](DESIGN_NOTES.md) — how this differs from the 2017 paper.
* [CHANGELOG.md](CHANGELOG.md) — what changed at each version.

## License

MIT.
