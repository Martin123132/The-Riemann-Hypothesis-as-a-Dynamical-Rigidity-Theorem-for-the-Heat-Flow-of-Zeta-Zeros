# Jensen-Window PF State-Space Sign-Lift Obstruction Scout

Date: 2026-07-06

Status: state-space sign-lift obstruction diagnostic. This is not a proof
against every state-space doubled model, modified production matrix,
oscillatory theorem, Schur positivity, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

## Purpose

The modified signed-model target leaves a genuine state-space doubled model
open. This note rejects only the cheapest such idea:

```text
keep the raw ordinary Motzkin paths
split positive and negative paths into sign states
replace signed path weights by absolute values
count all lifted paths with nonnegative readout
```

That absolute-value sign-state cover already fails at the first nontrivial
coefficient.

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json
work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_state_space_sign_lift_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_state_space_sign_lift_obstruction_scout.py
```

Current result:

```text
validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues
```

## Obstruction

The raw Motzkin/J-fraction coefficient has:

```text
mu_2 = beta_0^2 + lambda_1
kappa_1 = -lambda_1 > 0
mu_2 = beta_0^2 - kappa_1
```

An absolute-value sign-state cover would instead contribute:

```text
beta_0^2 + kappa_1
```

The gap is:

```text
absolute_lift_mu2 - mu_2 = 2*kappa_1 = -2*lambda_1 > 0
```

So this lifted nonnegative cover cannot equal the actual reciprocal
coefficient.

## Finite Grid

This scout is derived from the validated Motzkin obstruction rows:

```text
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

It validates:

```text
735 / 735 derived mu_2 sign-lift rows have positive absolute-lift gap
```

with:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
degrees d=2..8
```

## Interpretation

This closes another cheap escape hatch. A state-space doubled model cannot be
only a sign-state absolute-value cover of the raw Motzkin path expansion.

A surviving state-space route must change more than signs. It must give a new
exact representation of:

```text
E(t)=1/H(-t)
```

whose final path, network, or coefficient formula is manifestly nonnegative
and equals every `mu_m`, not only a finite prefix.

## Integration

This sharpens:

```text
msm_04_state_space_doubled_model
rp_04_companion_or_production_matrix_total_positivity
rp_09_signed_or_modified_continued_fraction
outputs/jensen_window_pf_modified_signed_model_target.md
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
outputs/jensen_window_pf_signed_j_fraction_theorem_target.md
```

## Boundary

Passing this checker rejects only the absolute-value cover of raw Motzkin path
weights. It does not rule out a genuinely modified state-space doubled model,
modified production matrix, oscillatory sign-regular theorem, Xi/Phi positive
cancellation identity, all-order column recurrence theorem, Schur positivity,
Jensen-window PF-infinity, or `Lambda <= 0`.
