# Jensen-Window PF Compound Order-Six Lambda=-100 Prefix Certificate

Date: 2026-07-13

Status: rigorous signed order-six endpoint prefix through `n=316` with
one open analytic tail. This is not a proof of all-shift order six,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_m100_prefix_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_prefix_certificate.py
```

## Stable Endpoint Coordinate

Inside the completed positive order-five cone, signed condensation gives

```text
Q_(6,n)*H_(4,n+2)=H_(5,n+1)^2-H_(5,n)*H_(5,n+2)
K_n=H_(5,n+1)^2/(H_(5,n)*H_(5,n+2))-1
Q_(6,n)=H_(5,n)*H_(5,n+2)*K_n/H_(4,n+2)
```

Thus `Q_(6,n)>0` is exactly `K_n>0`. Each `H_(5,n)` is evaluated
through the cancellation-preserving factorization `H_(5,n)=W_n*J_n`,
not by subtracting a raw six-by-six determinant.

## Coefficient Cover

The inherited sources, local repair, and two-row splice extension,
merged in precedence
order and hashed in the JSON artifact, give 1024-bit outward-rounded
balls for

```text
A_k(-100)>0 for every 0<=k<=326.
```

The dedicated repair covers `A_191,...,A_229` at 220 decimal digits.
It removes 24 interval-indeterminate rows caused by the older broad
source; no midpoint sign was promoted.

## Prefix Theorem

Direct stable Arb evaluation proves

```text
K_n(-100)>0 for every 0<=n<=316,
Q_(6,n)(-100)>0 for every 0<=n<=316.
```

The weakest row is

```text
n=316,
K_316=[0.00780960782881022083463485 +/- 3.80E-27],
K_316 lower=7.80960782881022083463484620439881399269827010621568534633766E-3>7/1000.
```

## Remaining Tail

The complete endpoint problem is now reduced to

```text
Q_(6,n)(-100)>0 for every n>=317.
```

The uniform eventual signed tail from the companion reduction is
non-effective, so it does not splice this finite prefix to every shift.

```text
outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
