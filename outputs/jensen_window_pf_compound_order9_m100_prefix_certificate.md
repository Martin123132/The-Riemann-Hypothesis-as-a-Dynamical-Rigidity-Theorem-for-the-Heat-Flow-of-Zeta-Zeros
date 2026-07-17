# Jensen-Window PF Compound Order-Nine Lambda=-100 Prefix Certificate

Date: 2026-07-13

Status: rigorous signed order-nine endpoint prefix through `n=1240`
with one open analytic tail. This is not a proof of all-shift order
nine, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order9_m100_prefix_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order9_m100_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_prefix_certificate.py
```

## Stable Endpoint Coordinate

```text
Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)
M_n=Q_(8,n+1)^2/(Q_(8,n)*Q_(8,n+2))-1
Q_(9,n)=Q_(8,n)*Q_(8,n+2)*M_n/Q_(7,n+2)
```

Thus `Q_(9,n)>0` is exactly `M_n>0`. No raw nine-by-nine
near-cancelling determinant is evaluated.

## Coefficient Cover And Repair

The inherited sources and two rigorous retained-integral repairs give
2048-bit Arb coefficient balls through `A_1256`. The new repair covers
`A_153,...,A_178`; the inherited repair covers `A_179,...,A_190`.
All 38 repaired rows use `n_sum=70`, cutoff 7, and 220 decimal digits.

## Prefix Theorem

```text
Q_(8,n)(-100)>0 for every 0<=n<=1242,
M_n(-100)>0 for every 0<=n<=1240,
Q_(9,n)(-100)>0 for every 0<=n<=1240.
```

The weakest row is the final cached shift:

```text
n=1240,
M_1240=[0.00406590691191712528081394310994514114249520440684677115477570 +/- 3.65E-63],
M_1240 lower=4.06590691191712528081394310994514114249520440684677115477569E-3>1/250.
```

## Remaining Tail

```text
Q_(9,n)(-100)>0 for every n>=1241.
```

```text
outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order8_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
