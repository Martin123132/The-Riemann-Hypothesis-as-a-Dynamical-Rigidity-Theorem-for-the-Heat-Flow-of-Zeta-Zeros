# Jensen-Window PF Compound Order-Seven Lambda=-100 Entry Certificate

Date: 2026-07-13

Status: all-shift signed contiguous order-seven entry theorem at
`lambda=-100`. This certificate is not a proof of forward order-seven
invariance, arbitrary-column order seven, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order7_m100_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_entry_certificate.py
```

## Continuous First-Summand Theorem

The fourth-nested stable curvature is covered on four closed ranges:

```text
r_1''(t)<=600/t^2 for 320<=t<=1000
r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2)
r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20
t^2*r_1''(t)<100<600 for u>=20
```

The asymptotic scaled upper is `5.55401953430641442537307739257756068343861012881347792016427E+1`.
The certified mode-to-t coverage gives

```text
r_1''(t)<=600/t^2 for every real t>=320
[320<=t<=1000] union [1000<=t<=V'(2)] union [2<=u<=20] union [u>=20] covers every real t>=320
```

## Discrete And Full-Kernel Transfer

```text
R_k^(1)<600*[-log(1-1/k^2)]<601/k^2, k>=321
|R_k-R_k^(1)|<=N_(k-1)+2*N_k+N_(k+1)<262/k^2, k>=321
R_k<R_k^(1)+|R_k-R_k^(1)|<601/k^2+262/k^2=863/k^2<900/k^2, k>=321
```

The last line is strictly below the scalar ceiling isolated by the
tail reduction, so

```text
R_k<900/k^2 => Q_(7,n)(-100)>0, k=n+6, n>=315
```

## Entry Theorem

The 1024-bit prefix proves

```text
Q_(7,n)(-100)>0 for every 0<=n<=314
```

Combining prefix and tail gives

```text
Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every integer n>=0
```

The sole endpoint input in the conditional order-seven heat-flow
reduction is therefore available for the forward composition.

```text
outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md
outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
outputs/formal_core.md
```
