# Hardy block-20 discrete selector atlas

Date: 2026-08-06

Status: complete finite selector atlas for all 212 block-20 pivots; not a proof and analytic q remains open

## Source roster

The separate full-block telemetry fixture records both cubic branches at every
source pivot

```text
rae_j = 657064 + 420*(j-1),  1 <= j <= 212.
```

All 424 chains end in the same exact evaluator checkpoint state and the same
fifteen displayed Hardy values as the accepted 64-chain fixture.  The observer
therefore exposes source decisions without changing the finite calculation.

## The transition the prefix missed

The first 64 chains all had `MIT=2`.  Across the full block, however,

```text
MIT=1 direct-kernel chains : 50
MIT=2 one-step chains      : 374.
```

The 50 direct cases are scattered through the later roster rather than forming
one terminal tail.  They have one saved level, kernel length 104, and no q
recurrence.  They are recorded as their own exact selector word and are not
forced into the earlier `MIT=2` model.

For each of the 374 recursive chains, the Section 11.243 Arb constructor was
replayed at the exact saved center.  Every chain admits a nonzero-radius box
preserving the child integer/parity/orientation decisions, q sign and loop
branches, all six `PSI` paths, both `ERF` paths, frac endpoints, denominator,
and stationary discriminant.  The maximum shrink was
`12` halvings.

The complete discrete route therefore contains

```text
424 / 424 classified branch calls,
374 rigorous local q-selector cells,
 50 exact direct-kernel calls,
  0 unclassified calls.
```

Branch 1 has `211` maximal constant
selector-word segments and branch 2 has
`211`.  This high fragmentation is
why extrapolating the first 32 pivots was unsafe.

## What is and is not covered

This is the all-212 **discrete** selector atlas required by the actual evaluator.
It does not claim continuous coverage between pivots; Section 11.245 proves that
the saved endpoint cells cannot provide such a cover.  It also does not yet
replace the source `PSI`/`ERF` approximations with rigorous values or enclose the
analytic `W1`--`W5` combination and `q`.  Direct `MIT=1` classification is not an
external-error theorem.

The next target is route-specific: rigorously evaluate special functions and the
correlated q correction on the 374 occupied cells, while separately auditing the
50 direct kernels.  Recurrence accumulation and the outer Hardy representation
remain independent obligations.  No physical carrier, determinant sign,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion follows here.
