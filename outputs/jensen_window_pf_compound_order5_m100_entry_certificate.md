# Jensen-Window PF Compound Order-Five Lambda=-100 Entry Certificate

Date: 2026-07-13

Status: all-shift contiguous order-five entry theorem at `lambda=-100`.
This is not a proof of forward order-five invariance by itself, arbitrary-column order five, PF-infinity, RH, or `Lambda <= 0`.
This is not by itself a proof of the forward or arbitrary-column theorem.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_m100_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_entry_certificate.py
```

## Continuous First-Summand Theorem

The nested stable curvature is now covered on three overlapping ranges:

```text
q_1''(t)<=60/t^2 for 320<=t<=V'(2)
q_1''(t)<=60/t^2 for every mode 2<=u<=20
q_1''(t)<=60/t^2 for every mode u>=20
```

The asymptotic scaled upper is `9.15835270172636955976486206054574636687947832827158772903837E+0`. These ranges cover every real `t>=320`, hence

```text
q_1''(t)<=60/t^2 for every real t>=320
```

## Discrete And Full-Kernel Transfer

```text
C_n^(1)<=60*[-log(1-1/k^2)]<=60/(k^2-1)<63/k^2, k=n+4>=321
|C_n-C_n^(1)|<=E_(k-1)+2*E_k+E_(k+1)<=37/k^2, k=n+4>=321
C_n<=C_n^(1)+|C_n-C_n^(1)|<63/k^2+37/k^2=100/k^2, k>=321
```

The last line is exactly the scalar ceiling isolated by the prior tail reduction, so

```text
C_n<=100/(n+4)^2 => J_n(-100)>0 and H_(5,n)(-100)>0, n>=317
```

## Entry Theorem

The 1024-bit prefix proves

```text
H_(5,n)(-100)>0 for every 0<=n<=316
```

Combining prefix and tail gives

```text
H_(5,n)(-100)>0 for every integer n>=0
```

Thus `target_compound_order5_m100_entry` is discharged. The next composition is the already-derived uniform tail and cooperative flow from `lambda=-100` through `lambda=0`.

```text
outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
