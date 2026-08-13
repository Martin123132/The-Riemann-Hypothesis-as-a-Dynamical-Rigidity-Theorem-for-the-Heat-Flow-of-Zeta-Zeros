# Newman C1 Hardy Joint Correlated-Error Ledger Gate

Date: 2026-08-05

Status: finite two-height correlated-error diagnostic and scalar-extrapolation
obstruction validated; no physical-height error theorem or interval enclosure;
not a proof of RH.

## Residual modes

The accepted `h=0.01` residual vectors are `fast_main-direct_main` on the
fifteen shifts `-0.07,...,0.07`.  A discrete orthonormal polynomial basis on
those same points gives:

| Height | max sample error | residual energy above degree 3 | above degree 5 |
|---|---:|---:|---:|
| `10^10` | `0.00724738923861` | `1.11354e-06` | `2.74694e-09` |
| `10^12` | `0.00288990976752` | `1.43085e-05` | `8.59347e-09` |

The profiles are smooth but not scalar copies.  Their cosine similarity is
`0.927285625844`, while the best scalar multiple leaves
`0.374355` of the `10^12` L2 norm.

## Physical projections

All `24` saved physical coefficient vectors are projected against both error
profiles.  `21` annihilate a constant shift error to relative sensitivity
below `1e-12`.  The observed correlations reduce the independent L1 bounds by
factors ranging from `7.0262` to
`1.77172e+06` at `10^10`, and from
`11.1023` to
`416509` at `10^12`.

This cancellation is empirical, not a bound.  More importantly, the raw
maximum sample error falls by a factor
`0.398752`, while
`21/24` physical projections increase.
The largest increase factor is
`3.48541`.  A single scalar
`T`-decay envelope therefore does not control these coefficient-weighted
errors even on the two saved calibration heights.

## Cutoff stability

The physical kernel design matrices have condition numbers near `2.5e16`.
The conclusion above was replayed with least-squares cutoffs from `1e-8` to
`1e-16`.  Across all `6` choices, exactly `21/24` projections
increase.  The maximum projected errors remain in the ranges
`[0.0093747921, 0.0094619115]`
at `10^10` and
`[0.015933506, 0.01606433]`
at `10^12`, despite substantial coefficient-norm variation.

## Missing theorem

The admissible upgrade is a modewise external-error theorem.  For the fixed
shift grid, write the error as

```text
epsilon(T,delta_j)=sum_(k=0)^5 a_k(T) q_k(j)+r_6(T,j).
```

One needs explicit physical-height bounds for every `a_k(T)` and for the
rough remainder, followed by

```text
|sum_j c_j epsilon(T,delta_j)|
 <= sum_(k=0)^5 |<c,q_k>| A_k(T)+||c||_2 R_6(T).
```

A scalar bound on `max_j|epsilon_j|` throws away the observed cancellation;
a scalar asymptotic fitted from two heights is contradicted by the projected
rows.  No physical-height run is interpretable until this componentwise
theorem, or an interval-certified replacement, is available.

## Proof Boundary

This is a finite diagnostic on two saved error vectors and twenty-four saved
kernel rows.  It proves neither asymptotic mode bounds nor persistence of the
observed correlations.  It supplies no physical carrier value, interval
enclosure, retained observation, determinant sign or bound, current
inequality, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.

built Newman C1 Hardy joint correlated-error ledger gate: 2 residual vectors, 15 shifts, 24 physical rows, 192 restart records inherited, 21 constant-annihilating rows, 21 projected errors increase, 6 cutoff replays, 0 physical values and 1 open modewise-error theorem
