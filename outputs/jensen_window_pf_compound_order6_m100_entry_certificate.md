# Jensen-Window PF Compound Order-Six Lambda=-100 Entry Certificate

Date: 2026-07-13

Status: all-shift signed contiguous order-six entry theorem at
`lambda=-100`. This certificate is not a proof of forward order-six
invariance, arbitrary-column order six, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_m100_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_entry_certificate.py
```

## Continuous First-Summand Theorem

The third-nested stable curvature is covered on three overlapping ranges:

```text
p_1''(t)<=200/t^2 for 321<=t<=V'(2)
p_1''(t)<=200/t^2 for every saddle mode 2<=u<=20
p_1''(t)<=200/t^2 for every mode u>=20
```

The asymptotic scaled upper is `2.27683610696345567703247070312471784171958323101989922910682E+1`.
These ranges cover every real `t>=321`, hence

```text
p_1''(t)<=200/t^2 for every real t>=321
```

## Discrete And Full-Kernel Transfer

```text
P_k^(1)<=200*[-log(1-1/k^2)]<201/k^2, k>=322
|P_k-P_k^(1)|<=Y_(k-1)+2*Y_k+Y_(k+1)<100/k^2, k>=322
P_k<=P_k^(1)+|P_k-P_k^(1)|<201/k^2+100/k^2=301/k^2<320/k^2, k>=322
```

The last line is strictly below the scalar ceiling isolated by the
prior tail reduction, so

```text
P_k<=320/k^2 => Q_(6,n)(-100)>0, k=n+5, n>=317
```

## Entry Theorem

The 1024-bit prefix proves

```text
Q_(6,n)(-100)>0 for every 0<=n<=316
```

Combining prefix and tail gives

```text
Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every integer n>=0
```

Thus `target_compound_order6_m100_entry` is discharged. The next
composition is the already-derived uniform tail and cooperative flow
from `lambda=-100` through `lambda=0`.

```text
outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
