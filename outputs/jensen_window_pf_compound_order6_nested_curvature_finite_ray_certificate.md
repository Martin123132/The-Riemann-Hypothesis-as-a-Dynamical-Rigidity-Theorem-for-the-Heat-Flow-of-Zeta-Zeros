# Jensen-Window PF Compound Order-Six Nested Curvature Finite-Ray Certificate

Date: 2026-07-13

Status: rigorous first-summand order-six curvature theorem on
`2<=u<=20`. This is not a proof of the asymptotic ray, order-six entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py
```

## Dimensionless Cover

The exact signed cumulant corridors through order eight and the coarse
absolute corridors for orders nine and ten produce `H^(2),...,H^(10)`
on a strict `t+-4` collar. The stable hierarchy is evaluated after
rescaling by `z=1/t`, preserving the cancellations explicitly.

A short aligned quadrature extension handles the endpoint `u=2`; the
remaining interval uses a rational mode grid of width `1/1000`.

```text
p_1''(t)<=200/t^2 for every saddle mode 2<=u<=20
ray blocks=17999,
largest t^2*p_1'' upper=2.13362219228344925466417117085597240579680675033580438106657E+1,
weakest S_1 lower=3.08397652081474752151148936849098251684064049884291231388636E+0.
```

## Remaining Ray

```text
Prove p_1''(t)<=200/t^2 for u>=20.
```

```text
outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md
outputs/formal_core.md
```
