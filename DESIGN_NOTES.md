# Design notes


## New and improved CCVT

**What changed.** The optimizer was rewritten following Li et al. (2010) and
Ahmed & Deussen (2017): sites only exchange points with their Delaunay
neighbours rather than all points. 

## Pseudo random points per dot from 2 to 16

**What changed.** The default number of pseudo-random points per dot went from 2 to 16 to reduce clumping the initial set of dots.


## Promotion: per-group instead of global

**What changed.** The paper pools every dot DBSCAN flags into one CCVT run. Here,
 flagged dots are first split into groups of dots within a
linking radius of each other and each group gets its own constrained CCVT run. 

**Why.** A dot always lands at the centroid of the points it owns. Under global
pooling, one larger dot can be built from touching dots in two distant towns, and
it then lands somewhere between them, on empty ground. 


## Exact dot values and leftovers

**What changed.** Every dot of a class is worth exactly its class value. Where a
group of touching dots does not divide evenly, the remainder stays as dots of
their own class.

**Why.** Balzer's implementation allows capacities to differ by one, and pooled
nationally that is harmless. Per group it is not: a group of 14 promoted as a
single dot would make that dot worth 1.4 times its stated value. 


## Nudging overlapping dots

**What changed.** New optional step: overlapping dots are pushed apart, `top` for
the largest class only or `all` for every class, with larger dots harder to move
in proportion to their area.

## The unit value advisor

**What it is.** A preview, before any dots are placed, of what candidate unit
values would do at the chosen scale and dot sizes.


