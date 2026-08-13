# Hardy q-cell physical-spacing obstruction gate

Date: 2026-08-06

Status: rigorous adjacent endpoint-cell obstruction; not a proof and no continuous coverage claim

## Test

The physical adapter fixes the first 32 block-20 pivots at

```text
rae_j = 657064 + 420*(j-1).
```

For each branch and each adjacent pair, this gate compares the exact binary128
`a1` centers with the exact rational `a1` radii of the two saved q-selector
cells.  Because `a1` is reduced modulo one, the comparison uses the shorter
distance on the circle.  This is generous: each ordinary real interval is
embedded into a circular arc, so disjoint circular arcs imply disjoint original
cell projections.

For every one of the `62` branch-adjacent pairs,

```text
circular_distance(a1_j,a1_(j+1)) > radius_j + radius_(j+1).
```

The minimum certified projected gap is
`1.1216003255563527107367602509942289092491966117183336698E-1`.  The smallest ratio of
center distance to the sum of endpoint radii is
`2.0970230408097908230171012777146997371337099379833859979E+2`.

The raw `a1` formula is smooth for `rae/a>1`; modulo one it is a continuous map
to the circle.  A connected path joining two centers cannot stay in the union
of two disjoint endpoint arcs.  Therefore the two saved endpoint q cells cannot
by themselves certify the intervening continuous physical segment.  No
assumption about `fracL(rae)` is needed for this falsification because failure
in one coordinate already defeats four-dimensional box containment.

## What this changes

The local q cells are valid, but they are islands around discrete source calls,
not an adjacent interval cover.  Trying to stretch them across a physical step
of 420 is the wrong scaling strategy.

The source evaluator itself uses a discrete pivot roster, so continuous coverage
between consecutive pivots is not required merely to reproduce that finite
computation.  The economical next route is to prove the selector and analytic
margins directly on the discrete index `j`, first for all 212 pivots of block 20,
using modular/floor segmentation.  Adaptive intermediate q cells remain a valid
fallback only if a later theorem genuinely needs a continuous `rae` interval.

## Boundary

This gate does not invalidate the q-cell certificates.  It invalidates only the
proposal that adjacent endpoint cells already form a continuous physical cover.
It does not prove the remaining source pivots, special-function values, the
Riemann Hypothesis, or any prize-level conclusion.
