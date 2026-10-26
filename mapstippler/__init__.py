"""mapstippler: pseudo-random, blue-noise and graduated dot maps.

Pseudo-random, blue-noise, and graduated dot maps, after
Arnold, Jenny & White (2017), "Creation and evaluation of graduated dot maps",
International Journal of Geographical Information Science.

Command line:
    mapstippler advise polygons.parquet --value population
    mapstippler dots polygons.parquet --value population --unit 10000 --id geoid --preview
"""

from .allocate import allocate_points, largest_remainder
from .ccvt import CCVTOptimizer, CCVTResult, ccvt, reduce_points, split_capacity
from .ccvt_fast import FastCCVT, ccvt_fast, recursive_split
from .scale import MapScale

__all__ = [
    "CCVTOptimizer", "CCVTResult", "ccvt", "reduce_points", "split_capacity",
    "FastCCVT", "ccvt_fast", "recursive_split",
    "allocate_points", "largest_remainder", "MapScale", "make_dot_map",
]


def __getattr__(name):
    # geopandas is only needed for the map pipeline, so import it lazily
    if name == "make_dot_map":
        from .graduated import make_dot_map
        return make_dot_map
    raise AttributeError(name)
