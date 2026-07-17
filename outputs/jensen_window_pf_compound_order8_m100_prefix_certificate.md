# Jensen-Window PF Compound Order-Eight Lambda=-100 Prefix Certificate

Date: 2026-07-13

Status: rigorous signed order-eight endpoint prefix through `n=1242`
with one open analytic tail. This is not a proof of all-shift order
eight, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_m100_prefix_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_m100_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_prefix_certificate.py
```

## Stable Endpoint Coordinate

Inside the completed positive order-seven cone, signed condensation gives

```text
Q_(8,n)*Q_(6,n+2)=Q_(7,n+1)^2-Q_(7,n)*Q_(7,n+2)
M_n=Q_(7,n+1)^2/(Q_(7,n)*Q_(7,n+2))-1
Q_(8,n)=Q_(7,n)*Q_(7,n+2)*M_n/Q_(6,n+2)
```

Thus `Q_(8,n)>0` is exactly `M_n>0`. The calculation evaluates
`H_4`, `H_5`, `Q_6`, and `Q_7` through stable ratio coordinates;
it never subtracts a raw eight-by-eight determinant.

## Coefficient Cover

The inherited hashed sources, twelve-row precision repair, and two
retained-integral extensions totaling 930 rows give 2048-bit Arb evaluation
from outward-rounded coefficient balls for

```text
A_k(-100)>0 for every 0<=k<=1256.
```

The extension was generated directly by rigorous Arb quadrature.

## Prefix Theorem

Direct stable Arb evaluation proves

```text
M_n(-100)>0 for every 0<=n<=1242,
Q_(8,n)(-100)>0 for every 0<=n<=1242.
```

The weakest row is the final cached shift:

```text
n=1242,
M_1242=[0.00355609273479377717793692780206927653148465783678139703539852 +/- 1.36E-63],
M_1242 lower=3.55609273479377717793692780206927653148465783678139703539851E-3>1/300.
```

## Remaining Tail

The complete endpoint problem is reduced to

```text
Q_(8,n)(-100)>0 for every n>=1243.
```

The all-fixed-order eventual theorem supplies a finite but ineffective
tail threshold, so it cannot splice this prefix to every shift.

```text
outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
