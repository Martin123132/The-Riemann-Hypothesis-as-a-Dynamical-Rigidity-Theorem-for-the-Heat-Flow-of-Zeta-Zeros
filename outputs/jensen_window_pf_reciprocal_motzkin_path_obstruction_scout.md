# Jensen-Window PF Reciprocal Motzkin Path Obstruction Scout

Date: 2026-07-06

Status: ordinary Motzkin-path obstruction diagnostic. This is not a proof against every modified signed continued-fraction
or production-matrix model, not an all-order column recurrence theorem, not
Schur positivity, not Jensen-window PF-infinity, not RH, and not
`Lambda <= 0`.

## Purpose

The signed J-fraction and signed Jacobi beta scouts make the finite ordinary
J-fraction parameters precise. This note asks whether the raw ordinary
Motzkin path model is already a manifest positivity proof.

For an ordinary J-fraction:

```text
M(t) = 1 / (1 - beta_0*t - lambda_1*t^2/(1 - beta_1*t - ...))
```

the coefficient sequence is a weighted Motzkin-excursion path sum with
horizontal weights `beta_h` and down-step weights `lambda_h`. A nonnegative
production-matrix proof would need those path weights, or an equivalent
transformation, to be nonnegative.

The first obstruction is already:

```text
mu_2 = beta_0^2 + lambda_1
beta_0^2 > 0
lambda_1 < 0
mu_2 > 0
```

So the raw ordinary Motzkin expansion has a negative path contribution even
when the resulting coefficient is positive.

There is also a diagonal obstruction:

```text
beta_1 < 0 for d >= 3
```

and diagonal sign conjugation cannot change a diagonal entry.

## Outputs

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py
```

Current result:

```text
validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues
```

## Finite Grid

The Arb grid checks:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
degrees d=2..8
dps=520
```

It validates:

```text
735 / 735 mu_2 rows have beta_0^2 > 0, lambda_1 < 0, and mu_2 > 0
630 / 630 d>=3 beta_1 diagonal rows are negative
```

This sharpens:

```text
rp_04_companion_or_production_matrix_total_positivity
rp_09_signed_or_modified_continued_fraction
```

and is finite evidence only.

## Interpretation

The raw ordinary Motzkin/J-fraction production matrix is not a positive path
model for the reciprocal coefficients. Positivity of `mu_2` already comes
from cancellation between a positive `beta_0^2` contribution and a negative
`lambda_1` contribution.

This does not kill the signed route. It says the next viable theorem must use
something stronger or different, such as a parity-doubled state space, a
modified production matrix, an oscillatory/sign-regular theorem with an
explicit positivity conclusion, or an Xi/Phi-specific representation whose
weights are manifestly nonnegative after the correct transformation.

## Boundary

Passing this checker rejects only the raw ordinary Motzkin path model as a
manifest positivity proof. It does not rule out every modified signed
continued fraction, every production matrix, every oscillatory matrix theorem,
the all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
or `Lambda <= 0`.
