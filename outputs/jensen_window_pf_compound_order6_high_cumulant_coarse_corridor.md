# Jensen-Window PF Compound Order-Six High-Cumulant Coarse Corridor

Date: 2026-07-13

Status: global coarse exact ninth- and tenth-cumulant corridor theorem.
This is not a proof of order-six curvature, PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py
```

## Formal Bounds

The exact tilted-Gaussian partition recurrence is extended only far
enough to extract the scaled epsilon-ten cumulants `r=9,10`. Rational
termwise substitution gives

```text
|scaled kappa_r^[10]|<1600 for 2<=u<=20,
|scaled kappa_r^[10]|<36000 for u>=20.
```

These constants are intentionally loose; the order-six stable hierarchy
suppresses the new derivative uncertainty by several powers of `q`.

## Exact-Density Transfer

The existing complex-disk theorem already controls one analytic
partition residual on `|z|<=1`. Cauchy's estimate is not limited to
order eight. For the two new derivatives,

```text
r!/(r-2)!=r*(r-1)<=90 for r=9,10
finite scaled residual < 19/1375
ray scaled residual < 1099/110000/u
```

so the complete exact corridor is

```text
|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10, u>=2.
```

No sign is asserted or needed for these two coarse corridors.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md
outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md
outputs/formal_core.md
```
