# mapstippler

Dot maps for Python: pseudo-random, blue-noise, and graduated dots.

Implements and extends the method in:

> Arnold, N.D., Jenny, B. & White, D. (2017). Creation and evaluation of graduated
> dot maps. *International Journal of Geographical Information Science*.
> doi:10.1080/13658816.2017.1359747

---

## How it works

1. Pseudo-randomly place dots in each polygon.
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

Input: a polygon layer in an **equal area coordinate system**, with a column
of counts for the dots (not rates or percentages).

### Sample data
Sample data: the [latest release](https://github.com/nick-arnold/mapstippler/releases/tag/v0.7.4)
has Oregon and Indiana census tracts with population (2024 TIGER, ACS 5-year
estimates) in ESRI:102008. The commands below assume `oregon_tracts.parquet`.

### Choose a unit value

`advise` estimates a range of unit values without placing dots:

```
mapstippler advise oregon_tracts.parquet --value population --mpp 500
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
mapstippler dots oregon_tracts.parquet --value population --unit 1000 --mpp 500 \
    --id geoid --nudge all --out [file] --preview
```

Writes `dots.gpkg`, a record of the parameters, and a PNG at the true scale. `--format parquet` writes GeoParquet
instead, which is smaller.

### Other methods

```
mapstippler dots oregon_tracts.parquet --value population --unit 1000 --mpp 500 --method blue
mapstippler dots oregon_tracts.parquet --value population --unit 1000 --mpp 500 --method random
```

`blue` is one dot size in a blue noise pattern. `random` is pseudo-random placement.

### Two colors

Two columns color each dot by category, matching each polygon's split:

```
mapstippler dots counties.parquet --value [field 1] [field 2] --unit 10000 --mpp 3800 --preview
```

### A single feature

```
mapstippler dots oregon_tracts.parquet --value population --unit 1000 --mpp 500 \
    --where "statefp == '41'" --out [file] --preview
```

### Python

```python
import geopandas as gpd
from mapstippler import MapScale, make_dot_map
from mapstippler.output import write_outputs

polygons = gpd.read_parquet("oregon_tracts.parquet")
result = make_dot_map(polygons, "population", unit_value=1000,
                      scale=MapScale(500.0), id_column="geoid", nudge="all")

result.dots            # class, value, source polygon, diagnostics
result.reconciliation  # per polygon: its value against the dots inside it
write_outputs(result, "mymap/")
```

---

## License

MIT.
