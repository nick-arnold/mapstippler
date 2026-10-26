"""Map scale and dot sizes.

A graduated dot layer is designed for one fixed web-map zoom. The zoom fixes the
ground distance covered by one CSS pixel (meters per pixel, in the units of the
equal-area projection). Dot sizes are set in CSS pixels, so everything that
depends on real distances follows from those two:

    dot diameter in map units  = size_px * m_per_px

Two dots of a class are treated as merged when they are closer than that class's
search radius, which is what DBSCAN looks for and what decides promotion. By
default the radius is the dot diameter (Arnold et al. 2017): dots merge exactly
when they physically overlap.

Dot sizes trade off against each other in a way the 2017 paper does not discuss.
Promoting k small dots into one larger dot removes their ink and draws one dot in
its place, so the ink kept is (d_large / d_small)^2 / k. Too small a size ratio
empties out the places that were densest, leaving holes where the map should be
darkest; too large a ratio means the larger class coalesces at the same density
as the one below it, which is the paper's argument against area-proportional
dots. The two are exactly reciprocal: ink kept x extra density absorbed = 1. For
a value ratio of 10, a diameter ratio near 2.25 keeps about half the ink while
absorbing twice the density, which is where the defaults sit.

`radii_px` sets the radii independently. A radius larger than the diameter
promotes dots that are close but not quite touching, which is often what the eye
sees: at small sizes two dots with a fraction of a pixel between them read as one
blob. It also separates two decisions that are otherwise tangled — how big a dot
looks, and how close is too close — so the texture of the smallest class can be
tuned without changing how much gets promoted, and vice versa.

Zoom ladder: zoom 0 shows the lower 48 about 1,200 CSS pixels wide
(3,800 m/px in a continental Albers projection). Each zoom level halves the
meters per pixel.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["DEFAULT_SIZES_PX", "CONUS_ZOOM0_M_PER_PX", "MapScale"]

DEFAULT_SIZES_PX = (2.0, 4.5, 10.0)     # smallest to largest dot class, CSS pixels
CONUS_ZOOM0_M_PER_PX = 3800.0


@dataclass(frozen=True)
class MapScale:
    m_per_px: float                                 # map units (meters) per CSS pixel
    sizes_px: tuple = field(default=DEFAULT_SIZES_PX)
    class_ratio: float = 10.0                       # value ratio between consecutive classes
    radii_px: tuple | None = None                   # search radii in CSS px; None = the diameters

    def __post_init__(self):
        if self.m_per_px <= 0:
            raise ValueError("m_per_px must be positive")
        if list(self.sizes_px) != sorted(self.sizes_px):
            raise ValueError("dot sizes must go from smallest to largest")
        if self.radii_px is not None:
            r = tuple(float(x) for x in self.radii_px)
            if len(r) == 1 and len(self.sizes_px) > 1:
                # one number given: scale every class's diameter by it
                r = tuple(r[0] * s for s in self.sizes_px)
            if len(r) != len(self.sizes_px):
                raise ValueError("need one search radius per dot class, or a single multiplier")
            if any(x <= 0 for x in r):
                raise ValueError("search radii must be positive")
            object.__setattr__(self, "radii_px", r)

    @classmethod
    def zoom(cls, level: float, **kwargs) -> "MapScale":
        """Scale at a zoom level of the CONUS ladder (0 = lower 48 at ~1,200 px wide)."""
        return cls(CONUS_ZOOM0_M_PER_PX / 2.0 ** level, **kwargs)

    def single_class(self) -> "MapScale":
        """The same scale with only the smallest dot class, for single-class maps."""
        return MapScale(self.m_per_px, sizes_px=(self.sizes_px[0],),
                        class_ratio=self.class_ratio,
                        radii_px=None if self.radii_px is None else (self.radii_px[0],))

    @property
    def n_classes(self) -> int:
        return len(self.sizes_px)

    def diameter(self, cls_index: int) -> float:
        """Diameter of a dot class in map units (0 = smallest class)."""
        return self.sizes_px[cls_index] * self.m_per_px

    def radius_px(self, cls_index: int) -> float:
        """Search radius of a class in CSS pixels (default: the dot diameter)."""
        return self.sizes_px[cls_index] if self.radii_px is None else self.radii_px[cls_index]

    def search_radius(self, cls_index: int) -> float:
        """Two dots of this class closer than this are treated as merged, in map units."""
        return self.radius_px(cls_index) * self.m_per_px

    def value(self, unit_value: float, cls_index: int) -> float:
        return unit_value * self.class_ratio ** cls_index

    def ink_kept(self, cls_index: int) -> float:
        """Share of ink kept when class `cls_index` - 1 is promoted into this class.

        1.0 means the larger dot covers exactly as much as the dots it replaced;
        below that, promotion lightens the densest places.
        """
        if cls_index == 0:
            return 1.0
        r = self.sizes_px[cls_index] / self.sizes_px[cls_index - 1]
        return r * r / self.class_ratio

    def describe(self) -> dict:
        return {
            "m_per_px": self.m_per_px,
            "approx_scale_denominator": round(self.m_per_px / 0.00028),   # OGC 0.28 mm pixel
            "sizes_px": list(self.sizes_px),
            "radii_px": [self.radius_px(i) for i in range(self.n_classes)],
            "radii_are_diameters": self.radii_px is None,
            "diameters_m": [self.diameter(i) for i in range(self.n_classes)],
            "search_radii_m": [self.search_radius(i) for i in range(self.n_classes)],
            "class_ratio": self.class_ratio,
            "ink_kept": [self.ink_kept(i) for i in range(self.n_classes)],
        }
