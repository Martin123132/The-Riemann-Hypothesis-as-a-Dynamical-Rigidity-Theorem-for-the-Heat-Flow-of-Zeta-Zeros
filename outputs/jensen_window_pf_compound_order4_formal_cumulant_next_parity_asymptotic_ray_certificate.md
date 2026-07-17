# Jensen-Window PF Order-Four Next-Parity Asymptotic-Ray Certificate

Date: 2026-07-13

Status: exact analytic theorem for the global next-parity coefficient layer.
This is not a proof of the exact cumulant ray, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.py
```

## Leading Bounds

The q-leading jets satisfy

```text
P_0=1; P_(r+1)=(u/2)*P_r'+2u*P_r
coefficient floor+1/(100u) < C_r^(infinity) < coefficient ceiling-1/(100u)
0<L_r^(infinity)<8/5, r=3,...,10
```

After `u=20+v`, both buffered sides of all seven signed coefficient
bounds have coefficient-positive numerators. For the odd expressions,
nonvanishing-square gates plus rigorous Arb signs at `u=20` fix the sign
before squared inequalities are used.

## Exact Jet Transfer

The existing fourteen sign gates through order eight are reused. Four new
coefficient-positive gates prove the order-nine and order-ten extension:

```text
|V^(r)-q*P_r(u)|<=10000*u^10, r=2,...,10
eta(u)=2500*u^8/q<=10^-21
|L_r-L_r^(infinity)|<=16*eta(u), r=3,...,10
|L_r|<2, r=3,...,10
```

Exact coefficient norms then give

```text
|C_r-C_r^(infinity)|<=64000000*eta(u)<1/(200u)
```

This consumes less than half the leading `1/(100u)` buffer. Therefore the
signed coefficient bounds proved by the 1,800-block Taylor certificate on
`2<=u<=20` hold for every `u>=2`.

## Remaining Boundary

The first omitted formal layer is now global. The remaining analytic gap is
the actual exact-density error after epsilon order eight: a central
beyond-epsilon-eight remainder and two adaptive tails. No exact cumulant
corridor is promoted here.

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.md
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/formal_core.md
```
