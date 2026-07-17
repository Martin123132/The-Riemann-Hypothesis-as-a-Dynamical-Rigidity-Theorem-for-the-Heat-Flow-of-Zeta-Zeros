# Jensen-Window PF Compound Order-Five Lambda=-100 Tail Curvature Reduction

Date: 2026-07-13

Status: exact order-five `lambda=-100` tail reduction to one open
stable log-curvature ceiling. This is not a proof of the analytic tail,
all-shift order five, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order5_m100_tail_curvature_reduction`.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py
```

## Exact Scalar Target

Put `k=n+4` and define

```text
C_n=log(F_n*F_(n+2)/F_(n+1)^2)+log(d_k^2/(d_(k-1)*d_(k+1))), k=n+4
```

The stable factorization gives the exact equivalence

```text
J_n>0 iff C_n<-4*log(x_k)
```

and the finite-difference form is

```text
C_n=Delta^2 log(F_n)-Delta^2 log(d_(n+3))
```

Thus the raw nested determinant has been reduced to one difference of
stable log-curvatures.

## Coarse Sufficient Ceiling

The proved lambda=-100 defect anchor gives

```text
-4*log(x_k)>=4*d_k>=502/(125*(2*k+1)), k>=320
```

A very loose sufficient curvature theorem is

```text
C_n<=100/k^2 for every k=n+4>=321
```

because exact rational arithmetic gives

```text
100/k^2<502/(125*(2*k+1)), k>=321
after k=321+m: 502*m**2 + 297284*m + 43689082>0.
```

Every coefficient of the shifted polynomial is positive. Therefore the
curvature ceiling implies `C_n<-4log(x_k)`, hence `J_n>0`, on the
whole tail.

## Boundary Calibration

The Arb prefix gives at its final row

```text
at n=316, (n+4)^2*C_n=3.5869277550969014082... while the sufficient cap is 100
relative_316=[0.00626902754512557813924681177 +/- 7.83E-30].
```

The proposed constant `100` has a factor-above-27 reserve over the
observed scaled curvature at the splice; it is chosen for analytic
robustness, not numerical sharpness.

## Conditional Completion

```text
[C_n<=100/(n+4)^2 for every n>=317] => [J_n(-100)>0 and H_(5,n)(-100)>0 for every n>=317]
prefix 0<=n<=316 plus the curvature tail => H_(5,n)(-100)>0 for every n>=0
the completed uniform tail and cooperative flow then imply H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0
```

The sole new analytic target in this reduction is

```text
prove C_n<=100/(n+4)^2 for every n>=317 at lambda=-100
```

No finite data or eventual non-effective asymptotic is promoted into
that ceiling.

```text
outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
