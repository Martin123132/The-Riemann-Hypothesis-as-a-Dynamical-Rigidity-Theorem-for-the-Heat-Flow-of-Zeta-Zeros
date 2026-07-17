# Jensen-Window PF Compound Order-Eight Nested-Curvature Compact Certificate

Date: 2026-07-13

Status: rigorous first-summand curvature theorem on `999<=t<=V'(2)`.
This is not a proof of the finite/asymptotic saddle rays, order-eight
entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py
```

## Common-Collar Core

Four aligned 107,452-row caches and one right collar enclose H2-H14
on a strict t+-6 neighborhood. Five stable logarithms are evaluated
with outward-rounded common-collar arithmetic, and every accepted
block proves B,J,R,S,T,U>0 before bounding s_1''.

## Compact Theorem

```text
s_1''(t)<=4000/t^2 for every real 999<=t<=V'(2),
all 96 adaptive blocks pass,
largest scaled upper=3.99455119863063431363867738935234812098029095541283895880362E+3<4000,
weakest U lower=1.22933158930953976688262799664500547709540427726036454298503E-4>0.
```

The first block overlaps the shifted-jet theorem at t=999. The final
block reaches saddle mode u=2 with no mode gap.

## Remaining Ray

```text
prove s_1''(t)<=4000/t^2 for every saddle mode 2<=u<=20,
then prove the asymptotic ray u>=20.
```

```text
outputs/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.md
outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md
outputs/formal_core.md
```
