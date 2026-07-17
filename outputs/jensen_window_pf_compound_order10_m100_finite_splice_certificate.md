# Order-Ten Lambda=-100 Finite Splice Certificate

Date: 2026-07-16

Status: rigorous two-index endpoint splice. This is not a proof of
the analytic order-ten tail, delayed entry, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_m100_finite_splice_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_m100_finite_splice_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_finite_splice_certificate.py
```

## Added Coefficients

Retained-integral quadrature with `n_sum=70`, cutoff 7, 220 decimal
digits, and absolute tolerance `1e-150` gives rigorous balls for
`A_1259(-100)` and `A_1260(-100)`.

## Direct Signs

At 4096-bit precision, signed Hankel determinants, Jacobi-Trudi
determinants, Toda numerators, and relative `Q_9` margins all give
the same strict signs:

```text
n=1241: L_n=[0.004566462556691569913871292746695554317516566763851582361307629401337215646269214861788723 +/- 9.24E-91]; Q_(10,n)>0
n=1242: L_n=[0.004563235162986120518706369887288213379505210924489440582398810375415555830810838269783069 +/- 2.73E-91]; Q_(10,n)>0
```

## Splice

```text
Q_(10,n)(-100)>0 for n=1241,1242
Q_(10,n)(-100)>0 for every 4<=n<=1242
Q_(10,n)(-100)<0 for n=0,1,2,3
n>=1243, equivalently k=n+9>=1252
```

The finite theorem now reaches two rows beyond the old scan. The
remaining endpoint task is a full-kernel continuum curvature bound
starting at `n=1243`; no sign is inferred there by this artifact.
