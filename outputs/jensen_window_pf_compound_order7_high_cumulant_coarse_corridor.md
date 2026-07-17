# Jensen-Window PF Compound Order-Seven High-Cumulant Coarse Corridor

Date: 2026-07-14

Status: global coarse exact eleventh- and twelfth-cumulant corridor
theorem. This is not a proof of order-seven curvature, PF-infinity,
RH, or `Lambda <= 0`.

## Formal Bounds

The exact epsilon-ten tilted-Gaussian recurrence needs no potential
orders beyond the already-certified `L_3,...,L_12` box. Rational
termwise substitution gives

```text
|scaled kappa_r^[10]|<14000 for 2<=u<=20,
|scaled kappa_r^[10]|<700000 for u>=20,
formal terms=72.
```

The exact Cauchy factor is `12*11=132`; both scaled residual budgets
are below one. Hence the finite and asymptotic source boxes give

```text
|kappa_r|*q^(r/2-1)/(r-2)!<14001, r=11,12, 2<=u<=20
|kappa_r|*q^(r/2-1)/(r-2)!<700001, r=11,12, u>=20
|kappa_r|*q^(r/2-1)/(r-2)!<1000000, r=11,12, u>=2
```

```text
work/rh_compute/results/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.json
outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md
outputs/formal_core.md
```
