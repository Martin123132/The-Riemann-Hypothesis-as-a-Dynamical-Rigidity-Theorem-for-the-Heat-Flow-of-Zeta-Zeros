# Hardy block-20 paired exceptional-cancellation gate

Date: 2026-08-07

Status: rigorous all-374 finite-roster cancellation audit; not a height-uniform proof and not a proof of RH

The exact correction is grouped according to the paper approximation it
replaces:

```text
Q_exact-Q_paper=(R_b-W2)+(R_c+zero_mode-(W3 or W4)).
```

Both grouped endpoint families sum to the exact correction on all 374 calls.
The maximum identity gap is
`8.71336313038452263401989483782017487101256847381592E-14`.  Pairing reduces the worst
component-triangle/correction ratio from
`1.38940364315770530666447663559917886409495020516720E+5` to
`2.10274136979634595059826824866236778629728459933153E+2`, a factor
of at least `6.60758219301250250144007846506532286901389920874051E+2`.  The
remaining worst witness is chain
`260`.

Thus W2--W4 must be bounded jointly with the exact endpoint remainder they
replace.  Even after that pairing, separate absolute values for the upper and
lower endpoint families can lose a factor above 210, so a competitive uniform
theorem must preserve their common Euler--Maclaurin structure or cancellation.
These finite ratios are diagnostics, not fitted constants, recurrence budgets,
or proofs of `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
