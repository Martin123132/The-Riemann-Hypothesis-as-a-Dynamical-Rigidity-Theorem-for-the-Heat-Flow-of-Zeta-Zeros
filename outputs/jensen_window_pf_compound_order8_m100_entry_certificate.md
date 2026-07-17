# Jensen-Window PF Compound Order-Eight Lambda=-100 Entry Certificate

Date: 2026-07-13

Status: all-shift signed contiguous order-eight entry theorem at
`lambda=-100`. This certificate is not a proof of forward order-eight
invariance, arbitrary-column order eight, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_m100_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_m100_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_entry_certificate.py
```

## Continuous First-Summand Theorem

The fifth-nested stable curvature is covered on four closed ranges:

```text
s_1''(t)<=2000/t^2 for every 699<=t<=999
s_1''(t)<=4000/t^2 for every real 999<=t<=V'(2)
s_1''(t)<=4000/t^2 for every saddle mode 2<=u<=20
t^2*s_1''(t)<200<4000 for u>=20
```

The finite saddle interval is covered by 17999 rational mode blocks.

The asymptotic scaled upper is `1.34489839907184243202209472656240124460170043090398255726938E+2`.
The certified mode-to-t coverage gives

```text
s_1''(t)<=4000/t^2 for every real t>=699
[699<=t<=999] union [999<=t<=V'(2)] union [2<=u<=20] union [u>=20] covers every real t>=699
```

## Discrete And Full-Kernel Transfer

```text
s_1''(t)<=4000/t^2 => W_k^(1)<=4000*[-log(1-1/k^2)]<4001/k^2, k>=1250
|W_k-W_k^(1)|<=C_(k-1)+2*C_k+C_(k+1)<190/k^2, k>=1250
W_k<W_k^(1)+|W_k-W_k^(1)|<4001/k^2+190/k^2=4191/k^2<4300/k^2, k>=1250
```

The last line is strictly below the scalar ceiling isolated by the
tail reduction, so

```text
W_k<4300/k^2 => Q_(8,n)(-100)>0, k=n+7, n>=1243
```

## Entry Theorem

The rigorous retained-integral prefix proves

```text
Q_(8,n)(-100)>0 for every 0<=n<=1242
```

Combining prefix and tail gives

```text
Q_(8,n)(-100)=H_(8,n)(-100)>0 for every integer n>=0
```

The sole endpoint input in the conditional order-eight heat-flow
reduction is therefore available for the forward composition.

```text
outputs/jensen_window_pf_compound_order8_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md
outputs/formal_core.md
```
