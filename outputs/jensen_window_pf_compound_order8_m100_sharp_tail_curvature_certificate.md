# Order-Eight Sharp Tail Curvature Certificate

Date: 2026-07-13

Status: rigorous order-eight sharp tail theorem at `lambda=-100`.
This is not a proof of order-nine entry, PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.py
```

## Recovered Reserve

The same aligned H2-H14 caches are adaptively repartitioned from
`t=1249` to saddle mode two. Every common-collar block proves all
five stable-log arguments positive before bounding the curvature.

```text
s_1''(t)<=2500/t^2 for every real 1249<=t<=V'(2)
all 95 sharp compact blocks pass,
largest scaled upper=2.49866277771572865399021660209073411363375179357860310142842E+3<2500,
t^2*s_1''(t)<=355.867<2500 for every saddle mode 2<=u<=20
t^2*s_1''(t)<200<2500 for u>=20
s_1''(t)<=2500/t^2 for every real t>=1249
```

## Complete Kernel

```text
W_k^(1)<=2500*[-log(1-1/k^2)]<2501/k^2, k>=1250
-log(1-x)<x/(1-x) and 2500/(k^2-1)<2501/k^2 because k^2>2501
|W_k-W_k^(1)|<190/k^2, k>=1250
W_k<2501/k^2+190/k^2=2691/k^2, k>=1250
```

The earlier 4300 ceiling was sufficient for order eight but too loose
for the sixth stable-log floor. The certified 2691 ceiling preserves
the reserve needed by the order-nine first/full bridge.

```text
outputs/jensen_window_pf_compound_order8_m100_entry_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md
```
