# Jensen-Window PF Order-Four Exact Cumulant Remainder Budget

Date: 2026-07-13

Status: exact formal composition and sharpened exact-density remainder target.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_remainder_budget`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.py
```

## Finite Budget

Arb gives `q(2)>9000`. The epsilon-six formal corridor
margin is greater than `9/500` throughout
`2<=u<=20`, while the global next-parity coefficient theorem gives

```text
scaled |kappa_r^[8]-kappa_r^[6]|<1/1000, r=2,...,8.
```

Thus the epsilon-eight formal model remains at least `17/1000` inside
every corridor. It is sufficient to prove

```text
scaled |kappa_r-kappa_r^[8]|<1/100, 2<=u<=20.
```

This still reserves a final corridor margin of `7/1000`.

## Ray Budget

The epsilon-six formal ray lies `1/(20u)` inside every corridor. Since
the largest next-parity coefficient cap is `61/10`, exponential q growth
proves

```text
scaled |kappa_r^[8]-kappa_r^[6]|<1/(100u)
```

The epsilon-eight formal model is therefore at least `1/(25u)` inside
every corridor. It is sufficient to prove

```text
scaled |kappa_r-kappa_r^[8]|<1/(50u), u>=20.
```

## Epsilon-Ten Sharpening

The globally certified second-next coefficients give

```text
scaled |kappa_r^[10]-kappa_r^[8]|<1/10000000, 2<=u<=20,
scaled |kappa_r^[10]-kappa_r^[8]|<1/(1000u), u>=20.
```

It is therefore sufficient to prove the cancellation-preserving residual
bounds

```text
scaled |kappa_r-kappa_r^[10]|<9/1000, 2<=u<=20,
scaled |kappa_r-kappa_r^[10]|<1/(100u), u>=20.
```

These leave final corridor margins greater than `7/1000` on the finite
interval and `29/(1000u)` on the asymptotic ray.

## Remaining Boundary

These are sufficient budgets, not exact-density estimates. The remaining
work is a cancellation-preserving central residual after epsilon ten
plus rigorous left and right tails. Raw high-moment interval boxes are not
used because their dependency widths destroy the q-scaled cancellations.

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/formal_core.md
```
