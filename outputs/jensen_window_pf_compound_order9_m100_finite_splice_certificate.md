# Order-Nine Lambda=-100 Finite Splice Certificate

Date: 2026-07-13

Status: rigorous two-index endpoint splice. This is not a proof of
the analytic order-nine tail, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order9_m100_finite_splice_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order9_m100_finite_splice_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_finite_splice_certificate.py
```

## Added Coefficients

Retained-integral quadrature with `n_sum=70`, cutoff 7, and 220
decimal digits gives rigorous balls for `A_1257` and `A_1258`.
Together with the established sources this covers `A_0,...,A_1258`.

## Stable Signs

```text
Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)
M_n=Q_(8,n+1)^2/(Q_(8,n)*Q_(8,n+2))-1
Q_(9,n)(-100)>0 for n=1241,1242
Q_(9,n)(-100)>0 for every 0<=n<=1242
minimum relative margin at n=1242: [0.00406015342664670077197221967240111812209211818491693044255200 +/- 4.35E-63]
```

The finite prefix now meets the analytic bridge with no index gap.
The remaining handoff is the continuous sixth-nested theorem
`w_1''(t)<=4200/t^2` on `t>=1250`.

```text
outputs/jensen_window_pf_compound_order9_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order9_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md
```
