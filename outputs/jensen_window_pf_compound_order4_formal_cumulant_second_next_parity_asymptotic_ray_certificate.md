# Jensen-Window PF Order-Four Second-Next Asymptotic-Ray Certificate

Date: 2026-07-13

Status: exact analytic theorem for the global second-next coefficient layer.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.py
```

## Leading Bounds

```text
P_0=1; P_(r+1)=(u/2)*P_r'+2u*P_r
coefficient floor+1/(100u) < D_r^(infinity) < coefficient ceiling-1/(100u)
0<L_r^(infinity)<2-1/(100u), r=3,...,12
```

After `u=20+v`, both sides of all seven buffered coefficient bounds
have coefficient-positive numerators. Odd signs are fixed by
nonvanishing-square gates and rigorous Arb endpoint signs.

## Exact Jet Transfer

The eighteen sign gates through order ten are reused. Four new gates give

```text
|V^(r)-q*P_r(u)|<=100000*u^12, r=2,...,12
eta(u)=25000*u^10/q<=10^-17
|L_r-L_r^(infinity)|<=20*eta(u), r=3,...,12
|L_r|<2, r=3,...,12
|D_r-D_r^(infinity)|<=8000000000*eta(u)<1/(200u)
```

The transfer consumes less than half the leading `1/(100u)` buffer.
Together with the 3,600-block finite Taylor certificate, every signed
second-next coefficient bound holds for all `u>=2`.

## Remaining Boundary

The complete epsilon-ten subtraction is now globally bounded. What remains
is the exact-density theorem itself: intervalize only the residual after
this subtraction, and prove the left and right adaptive tails separately.
No exact cumulant corridor is promoted here.

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/formal_core.md
```
