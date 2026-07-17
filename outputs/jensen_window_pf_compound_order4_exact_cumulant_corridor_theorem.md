# Jensen-Window PF Order-Four Exact Cumulant Corridor Theorem

Date: 2026-07-13

Status: global exact cumulant corridor theorem with open localized-curvature composition.
This is not a proof of the curvature ray, order-four entry, PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.py
```

## Exact Theorem

The finite and asymptotic formal corridor theorems, the two globally
bounded correction layers, and the exact central-and-tail residual compose
to prove all seven candidate corridors for every `u>=2`:

```text
2/5<=q*(kappa_2-1)<=4/5
1<=(-1)^3*kappa_3*q^(3/2-1)/1<=6/5
1<=(-1)^4*kappa_4*q^(4/2-1)/2<=27/20
1<=(-1)^5*kappa_5*q^(5/2-1)/6<=3/2
1<=(-1)^6*kappa_6*q^(6/2-1)/24<=17/10
1<=(-1)^7*kappa_7*q^(7/2-1)/120<=2
1<=(-1)^8*kappa_8*q^(8/2-1)/720<=5/2
```

The finite common reserve is `79999/10000000`; the ray reserve
is `29/(1000u)`.

## Remaining Boundary

The exact-density theorem is closed, but the curvature ray is not yet
closed. The remaining task is a continuum interval/analytic composition
of these exact corridor boxes with the exact mode geometry in the localized
quantity `U(t)`. Seven sampled collars passed previously; samples are not
a substitute for the full `u>=2` theorem.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/formal_core.md
```
