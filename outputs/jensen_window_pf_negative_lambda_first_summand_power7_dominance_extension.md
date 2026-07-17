# Jensen-Window PF First-Summand Power-Seven Dominance Extension

Date: 2026-07-13

Status: rigorous inverse-seventh-power first-summand dominance at
`lambda=-100`. This is not a proof of order six, PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py
```

## Strengthened Tail

The exact summand-ratio monotonicity and tilted-integrand geometry from
the power-six proof are unchanged. Re-running only the endpoint power
comparison at `k=316` gives

```text
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<=2/k^7 for every integer k>=316
0<=log(1+delta_k)<=2/k^7, k>=316
```

The two decisive outward-rounded endpoint balls are

```text
high log comparison=[-124.4151088715342004284394386786649561154590148631923285631932424956892 +/- 4.89E-68]<0,
low log comparison=[-0.09644772993769165999532627795735302274578455572009859714240907548187765 +/- 1.34E-72]<0.
```

All eleven derivative and endpoint gates are strictly positive. The
same ratio-log derivative checks used in the original certificate show
that both comparisons strengthen for every `k>=316`.

## Wall Error

For the adjacent logarithmic coordinate this yields

```text
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-7)+2*j^(-7)+(j+1)^(-7)), j>=317
```

That extra inverse power is what keeps the third nested stable-log
transfer below a fixed multiple of `k^-2`.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md
outputs/formal_core.md
```
