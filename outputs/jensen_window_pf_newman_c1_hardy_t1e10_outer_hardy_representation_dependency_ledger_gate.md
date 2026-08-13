# Outer-Hardy representation dependency ledger

Date: 2026-08-09

Status: diagnostic validated; not a proof; Hardy-Z and Xi certification remain blocked

The fifteen saved source outputs have been traced from their exact target
definitions through the hybrid Fortran path.  The audit distinguishes exact
mathematical identities, finite enclosed internal columns, and open outer
columns.  It contains 4 exact rows, 1
enclosed internal row, and 13 open rows.

The decisive boundary is upstream of the recurrence work.  The source follows
the paper's hybrid equations (126)--(127), and the paper explicitly says that
representation is not exact.  Its relative `O(epsilon_t)` and comparable
Gaussian-sum errors are not supplied here with constants that can be summed
into a Hardy-Z remainder.  Therefore the input `et=0.005` is a requested
accuracy parameter, not a proved error theorem.

The current bound

```text
1.18534368330532563659717863467183140923914676684614E-3
```

remains valid for the finite saved-height corrected internal blocks-20--35
column.  It is deliberately not relabelled as `|Z_source-Z_exact|`.

The missing outer ledger includes shifted-phase and scale Taylor remainders,
the `H(t)` normalization, lower-alpha and high-alpha truncations, block zero,
blocks 1--19, the finite lower Riemann-Siegel arithmetic, omitted higher
Riemann-Siegel corrections and final remainder, and complete source rounding.
Every one of the fifteen output rows is therefore marked
`hardy_z_certified=false` and `xi_certified=false`.

The next finite falsification gate is an independent Arb Hardy-Z calibration
at the same fifteen heights.  That finite calibration is not evidence for RH;
it tests whether the saved estimator is even inside its requested tolerance.
The proof route is separate: start from an exact classical Riemann-Siegel
identity (or an exact Hardy integral), partition it exactly, and attach an
explicit constant-bearing remainder to every transformation.
