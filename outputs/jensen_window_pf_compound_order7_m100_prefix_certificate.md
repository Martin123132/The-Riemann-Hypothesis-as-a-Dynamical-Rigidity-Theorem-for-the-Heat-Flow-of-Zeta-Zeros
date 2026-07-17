# Jensen-Window PF Compound Order-Seven Lambda=-100 Prefix Certificate

Date: 2026-07-13

Status: rigorous signed order-seven endpoint prefix through `n=314`
with one open analytic tail. This is not a proof of all-shift order
seven, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order7_m100_prefix_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_prefix_certificate.py
```

## Stable Endpoint Coordinate

Inside the completed positive order-six cone, signed condensation gives

```text
Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)
L_n=Q_(6,n+1)^2/(Q_(6,n)*Q_(6,n+2))-1
Q_(7,n)=Q_(6,n)*Q_(6,n+2)*L_n/H_(5,n+2)
```

Thus `Q_(7,n)>0` is exactly `L_n>0`. The calculation evaluates
`H_4`, `H_5`, and `Q_6` through their exact stable ratio
factorizations; it never subtracts a raw seven-by-seven determinant.

## Coefficient Cover And Repair

The inherited sources and the new twelve-row local repair, merged in
precedence order and hashed in the JSON artifact, give 2048-bit Arb
evaluation from outward-rounded coefficient balls for

```text
A_k(-100)>0 for every 0<=k<=326.
```

The dedicated retained-integral repair covers `A_179,...,A_190` at
220 decimal digits with the established analytic n-tail and cutoff-tail
bounds. It removes all eight interval-inconclusive rows from the broad
source; no midpoint sign is promoted.

## Prefix Theorem

Direct stable Arb evaluation proves

```text
L_n(-100)>0 for every 0<=n<=314,
Q_(7,n)(-100)>0 for every 0<=n<=314.
```

The weakest row is the final cached shift:

```text
n=314,
L_314=[0.00938326115496060603172 +/- 4.96E-24],
L_314 lower=9.38326115496060603172271048196819607046278167988820319015194E-3>9/1000.
```

## Remaining Tail

The complete endpoint problem is now reduced to

```text
Q_(7,n)(-100)>0 for every n>=315.
```

The all-fixed-order eventual theorem supplies a finite but ineffective
tail threshold, so it cannot splice this prefix to every shift.

```text
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
