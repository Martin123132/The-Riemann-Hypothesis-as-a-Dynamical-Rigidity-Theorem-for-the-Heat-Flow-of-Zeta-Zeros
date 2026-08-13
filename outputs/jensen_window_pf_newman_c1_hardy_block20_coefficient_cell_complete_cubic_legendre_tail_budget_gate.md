# Hardy block-20 coefficient-cell complete cubic Legendre-tail budget

Date: 2026-08-09
Status: complete_cubic_legendre_tail_transported_over_all_selector_cells_and_combined_budget_closes; finite block-20 theorem only, not a proof of the complete evaluator or RH

## Stable Exact Remainder

For `s=sqrt(1+z)`, exact algebra factors the complete phase remainder as

```text
H-H3=(u^2/y) z^2 (s^2+4s+1)/[12(1+s)^4].
```

This form is nonnegative on every admitted cell and avoids subtracting two
nearly equal interval phases.  The common normalized child phase has unit
modulus and is removed before estimation.  Each term is therefore enclosed as

```text
|exp(i*tpm*(H-H3))/(1+z)^(1/4)
 -1/(1+z/4-3z^2/32)|.
```

The full-cell recurrence multiplier is applied only after the termwise child
sum, and the exact output weights are applied only after the local majorant.

## Result

```text
selector cells                           374
child indices                            753
precision overlaps                       374 / 374
maximum cell |z|                         1.29853160979109816253185272216796875000000000000000000000000000000000000000000000E-3
minimum cell 1+z                         9.98701468390208901837468147277832031250000000000000000000000000000000000000000000E-1
maximum cell phase tail                  8.43638063940943538909777998924255371093750000000000000000000000000000000000000000E-6
maximum cell/point tail inflation        4.41594821479731632618701505419083776887468746196395296037951968202526105645017549E+0
maximum transported cell tail            1.69448620739517091116194404388982545403400354434574724762406923871886315234412849E-4
maximum endpoint plus cell tail          2.36597199906583183576986751275047104430259514781929472476240692387188631523441285E-3
outputs below 0.005                      15 / 15
```

No signed output cancellation is used.

## Boundary

This closes only the coefficient-cell endpoint and complete cubic Legendre-tail
columns at finite block 20.  Parent/saddle recurrence arithmetic, source
rounding, other blocks, cross-block accumulation, the outer Hardy
representation, `Lambda<=0`, PF-infinity, RH, and a prize-level theorem remain
open.
