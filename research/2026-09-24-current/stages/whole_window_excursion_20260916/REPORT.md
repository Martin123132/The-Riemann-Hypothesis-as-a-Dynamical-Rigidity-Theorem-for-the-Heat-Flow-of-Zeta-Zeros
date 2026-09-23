# The entire N=2980 window has a positive paid margin

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
Private local research, 16 September 2026. No publication performed.

## Result

The previously uncomputed complement of the local gated excursion is now
included. On the complete original second window

    N=2980, r in [299/400,301/400], d=9, M=2976,

retain the same central gate X_+, positive weight W=(1+cos chi)^2, source
payment and actual retained-polynomial margin m=min_(u in [0,1]) Q/H.
Let P=int W X_+ m_+, L=int W X_+ (-m)_+, and Z=int W. Then

    (P-L)/Z > 1.24773 > 0,
    P > 0.0102391,
    0.0000281908 < L < 0.000158166,
    L/P < 0.015448.

These are outward-rounded consequences of the 8192-cell exact-Fraction
reaggregation. All negative regions, uncertain gates and uncertain minimizer
locations are retained. The signed integral is also enclosed directly; its
lower bound agrees with the separately paid positive-credit-minus-loss test.
No earlier local recovery credit has been added to this complete cover.

The negative cost over the complete window is more than eight million times
the former single-excursion upper bound 3.42839e-12. That is not a contradiction:
the former calculation concerned one local component. It shows why its small
cost could not stand in for the uncomputed complement.

The older minimum-of-Bernstein-controls score also passes, with normalized
lower bound greater than 1.24750. The scores remain separate. Transverse
refinement resolves pointwise uncertainty, but is not necessary for positivity
of this particular full-window average.

## What was calculated

- All 2976 finite-prefix terms, with their shared arithmetic phase.
- Every height in the window, using 4096 cells at 192 bits and 8192 at 256 bits.
- The whole transverse interval [0,1], by restricting both numerator and
  denominator to eight exact dyadic pieces on every height cell.
- Ten coefficients and their first two height derivatives at every centre;
  third-derivative errors bounded over the entire height window.
- Both signed integration and separate positive/negative accounting.

Cell contributions are whole-cell interval bounds, not sampled quadrature.
No uniqueness or fixed branch of the transverse minimizer is assumed.
The same inherited analytic source reserve is paid in the numerator.

## Verification

The six incoming archives verify 156/156 manifest entries and pass all 214
named tests. The old 409600-panel response quadrature and every earlier
source calculation were not repeated.

The new suite passes 36/36 tests. The two entire covers have overlapping
integral enclosures. A separate exact-Fraction implementation reaggregates
every cell and independently reconstructs five transverse restrictions per
cover. A separate mpmath.iv direct-term evaluation at five fixed new centres
agrees with all ten native finite-prefix coefficients. These five checks
are not an independent enclosure of the entire window.

Both full source covers use FLINT/Arb, at different precisions and grids;
they are not independent special-function implementations. The finite
Taylor-sum adapter also passes a 23-coefficient direct eight-term test and
agrees with both historical source-jet records in all 30 entries.

The initial exact reaggregation rejected slight outward widening of a weight
ball. Its failed receipt and original checker are retained. Rational
intersection with the analytically known [0,4] range fixed serialization
handling only; no source, score or geometry changed.

## Scope and next theorem

This removes the missing-complement obligation for ONE complete finite window.
Under the inherited source comparison, the positive average supplies a
suitable height in that window. It does not make the whole window pointwise
positive, count the exact number of negative excursions, or prove RH.

The next analytic target is a uniform aggregate inequality

    L_N < P_N

on an appropriate unbounded sequence of cutoffs, with all gates, transverse
branches and source errors retained. The observed 1.5448% upper bound is NOT
assumed at other cutoffs. Winding and the absolute integer-count closure are
still separate obligations. An additional successful finite scan would not
replace either theorem.

See DERIVATION.md, SUMMARY.json and the complete append-only cell records.
The sealed package's separate replay receipt records extraction, manifest,
test and result-reproduction checks actually performed.
