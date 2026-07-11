# Jensen-Window PF Modified Signed-Model Target

Date: 2026-07-06

Status: modified signed-model theorem-search audit. This is not a proof of
the all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The signed reciprocal J-fraction diagnostics now give a sharper situation:

```text
signed Hankel/J-fraction signatures look coherent on finite grids
raw ordinary Motzkin path positivity is obstructed
diagonal sign conjugation and global length-parity signs are obstructed
direct reciprocal coefficients are positive on a larger finite grid
```

This note classifies which modified signed-model routes remain alive after
those obstructions, and which routes are now dead as proof mechanisms.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_modified_signed_model_target.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_modified_signed_model_target.py
```

Current result:

```text
validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates
```

## Target

For:

```text
h_j = binom(d,j) * A_{n+j}(0)
g_j = h_j / h_0
H(t) = 1 + g_1*t + ... + g_d*t^d
E(t) = 1 / H(-t)
mu_m = [t^m] E(t)
```

the modified signed-model target is to prove:

```text
mu_m >= 0
```

for every `m`, degree `d`, and shift `n`, without assuming endpoint PF,
Jensen hyperbolicity, Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Model Rows

```text
msm_01_raw_ordinary_motzkin_j_fraction:
  rejected_by_motzkin_obstruction
  The ordinary J-fraction Motzkin path expansion has negative path weights
  and cannot be used as manifest positivity.

msm_02_diagonal_sign_conjugation:
  rejected_by_sign_repair_obstruction
  Closed path signs are invariant under diagonal sign conjugation, and the
  same-length mixed-sign witnesses survive.

msm_03_global_length_parity_sign:
  rejected_by_sign_repair_obstruction
  Same-length paths H0^m and U D H0^(m-2) have opposite signs, so a global
  length-parity sign cannot make both nonnegative.

msm_04_state_space_doubled_model:
  live_if_constructed
  A two-state or parity-resolved path model remains possible only if it gives
  exact nonnegative lifted path sums for every mu_m. The cheap absolute-value
  sign-state cover of raw Motzkin paths is rejected at mu_2.

msm_05_modified_production_matrix_or_riordan:
  live_if_constructed
  A nonnegative production matrix, Riordan array, or planar network remains
  possible only if it generates the actual reciprocal coefficients and states
  whether it is column-only or all-Schur/Toeplitz.

msm_06_oscillatory_sign_regular_resolvent:
  live_if_theorem_proved
  Oscillatory or sign-regular resolvent language remains possible only if a
  theorem converts the signed Jacobi/resolvent data into coefficientwise
  nonnegativity of E(t). The oscillatory/resolvent fit matrix validates 8
  fit rows, 0 ready rows: only a noncircular positive spectral transform or
  an Xi/Phi-specific positive resolvent kernel remains live. The
  positive-readout theorem target turns those two options into exactness,
  positivity, noncircularity, and all-Schur/Toeplitz-scope obligations.
  The ordinary positive moment version is blocked by Delta_2=-g_2.

msm_07_cancellation_positive_identity:
  live_if_identity_proved
  The finite evidence may be explained by an Xi/Phi-specific cancellation
  identity, but the final formula must be manifestly nonnegative.

msm_08_finite_grid_extrapolation:
  finite_evidence_only
  The signed-Hankel, beta, Motzkin, parity-lift, and 72,600 reciprocal
  coefficient rows guide theorem search but cannot become the theorem.

msm_09_endpoint_real_rooted_or_pf_input:
  circular_if_used_as_input
  Endpoint real-rootedness or Jensen-window PF would imply reciprocal
  positivity, but using it as input assumes the target.
```

## Required Conditions

A surviving modified signed model must:

```text
1. recover exactly E(t)=1/H(-t) for all degrees and shifts;
2. prove mu_m>=0 for every m, not only a finite prefix;
3. avoid endpoint PF, Jensen hyperbolicity, Laguerre-Polya, RH, and Lambda <= 0 assumptions;
4. state whether it proves column shapes only or also an all-Schur/Toeplitz lift;
5. survive the raw Motzkin path obstruction, beta_1 diagonal obstruction, and same-length mixed-sign path obstruction;
6. end with nonnegative terms, not merely a difference of positive sums.
```

## Evidence Anchors

The raw Motzkin obstruction is:

```text
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues
```

The parity/sign-lift obstruction is:

```text
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues
```

The signed J-fraction and beta evidence is:

```text
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json
outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json
```

It validates:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues
```

The state-space sign-lift obstruction is:

```text
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_state_space_sign_lift_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues
```

The oscillatory/resolvent theorem-fit matrix is:

```text
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py
```

It validates:

```text
validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows
```

The positive-readout theorem target is:

```text
outputs/jensen_window_pf_positive_readout_theorem_target.md
work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json
python work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py
```

It validates:

```text
validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes
```

The ordinary positive spectral moment interpretation is blocked by:

```text
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json
python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py
```

It validates:

```text
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
```

The extended reciprocal coefficient stress is:

```text
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_coefficient_extended_stress.py
```

It validates:

```text
validated Jensen-window PF reciprocal coefficient extended stress: 72600 rows, 0 issues
```

## Kill Gates

Reject a proposed proof if it:

```text
1. uses the raw ordinary Motzkin/J-fraction path expansion as manifest positivity;
2. repairs that raw model only by diagonal sign conjugation;
3. repairs that raw model only by a global path-length parity sign;
4. uses only an absolute-value sign-state cover of raw Motzkin paths as a state-space doubled model;
5. promotes finite signed-Hankel, beta, Motzkin, parity-lift, or reciprocal grids into an all-order theorem;
6. assumes endpoint real-rootedness, Jensen-window PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0;
7. proves only formal signed-fraction existence without coefficientwise nonnegativity of E(t);
8. claims to close jwpf_06 from a column-only argument without an all-Schur/Toeplitz lift.
```

## Boundary

Passing this checker means the signed-model theorem search has a clean
classification of dead repairs and live modified candidates. It does not
construct a state-space model, production matrix, oscillatory theorem,
positive cancellation identity, all-order column recurrence theorem, or
`Lambda <= 0` proof.
