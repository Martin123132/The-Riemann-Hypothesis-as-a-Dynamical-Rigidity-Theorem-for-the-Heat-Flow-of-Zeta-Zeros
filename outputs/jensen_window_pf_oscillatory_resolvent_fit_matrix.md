# Jensen-Window PF Oscillatory Resolvent Fit Matrix

Date: 2026-07-06

Status: oscillatory/resolvent theorem-search matrix. This is not a proof of
the all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The signed J-fraction route naturally suggests an oscillatory or sign-regular
Jacobi/resolvent theorem. This note separates what that language can and
cannot currently do for:

```text
E(t) = 1/H(-t)
mu_m = [t^m]E(t)
```

The target is coefficientwise positivity:

```text
mu_m >= 0
```

for every `m`, degree `d`, and shift `n`.

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py
```

Current result:

```text
validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows
```

## Fit Rows

```text
or_01_nonnegative_jacobi_or_production_resolvent:
  rejected_by_entry_signs
  Entrywise nonnegative Jacobi powers would be a clean path-sum proof, but
  the actual signed Jacobi parameters have lambda_n<0 and beta_1<0.

or_02_diagonal_similarity_to_nonnegative_matrix:
  rejected_by_sign_invariance
  Diagonal sign similarity cannot change beta_1, and closed path signs are
  invariant.

or_03_absolute_value_resolvent_majorant:
  rejected_by_mu2_gap
  The absolute-value path/resolvent cover gives beta_0^2+kappa_1 at mu_2
  instead of beta_0^2-kappa_1.

or_04_classical_oscillatory_matrix_spectrum:
  route_mismatch
  Classical oscillatory or sign-regular spectral conclusions do not by
  themselves prove the scalar moments e_0^T J^m e_0 are nonnegative.

or_05_indefinite_or_krein_space_moment_model:
  language_only
  Indefinite moment language describes the signature but does not imply
  coefficient positivity without a positive readout theorem.

or_06_positive_spectral_measure_after_transform:
  live_if_transform_proved
  A nonordinary noncircular transform to a positive scalar readout would work
  if its Taylor coefficients are exactly mu_m and it is verified for the zeta
  windows. The ordinary positive moment interpretation is rejected by
  Delta_2=-g_2.

or_07_xi_phi_resolvent_kernel:
  live_if_theorem_proved
  An Xi/Phi-specific positive resolvent kernel remains possible if it gives
  exact coefficients and a true positivity theorem.

or_08_finite_signed_jacobi_pattern:
  finite_evidence_only
  The finite signed lambda/beta pattern guides theorem search but cannot
  become the theorem.
```

## Required Conditions

A usable oscillatory/resolvent theorem must:

```text
1. construct an exact matrix, operator, or kernel model whose scalar
   resolvent or moment generating function is E(t)=1/H(-t);
2. prove mu_m>=0 for every m, not merely spectral reality or formal
   continued-fraction existence;
3. verify all hypotheses for the actual zeta heat-flow Jensen windows;
4. avoid endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 assumptions;
5. survive lambda_n<0, beta_1<0, raw Motzkin cancellation, same-length mixed
   signs, and the absolute-value sign-lift obstruction;
6. state whether it is column-only or supplies an all-Schur/Toeplitz lift.
```

## Evidence Anchors

Signed J-fraction evidence:

```text
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md
```

These validate:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues
```

Path and sign-lift obstructions:

```text
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
```

These validate:

```text
validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues
validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues
validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues
```

The broader theorem ecosystem is tracked in:

```text
outputs/sign_regularity_theorem_map.md
outputs/sign_regularity_theorem_fit_matrix.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/jensen_window_pf_signed_j_fraction_theorem_target.md
outputs/jensen_window_pf_modified_signed_model_target.md
```

The two live rows are sharpened into the positive-readout theorem target:

```text
outputs/jensen_window_pf_positive_readout_theorem_target.md
work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json
python work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py
```

It validates:

```text
validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes
```

The ordinary positive moment interpretation is blocked by:

```text
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json
python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py
```

It validates:

```text
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
```

## Kill Gates

Reject a proposed proof if it:

```text
1. uses formal J-fraction existence without coefficientwise nonnegativity;
2. treats oscillatory or sign-regular spectral conclusions as coefficient
   positivity without a positive scalar readout;
3. diagonal-conjugates the raw Jacobi model and ignores beta_1<0 or closed
   path sign invariance;
4. replaces signed path weights by absolute values and claims equality with E(t);
5. assumes endpoint real-rootedness, Jensen-window PF-infinity, Laguerre-Polya
   membership, RH, or Lambda <= 0;
6. claims to close jwpf_06 from a column-only resolvent theorem without an
   all-Schur/Toeplitz lift.
```

## Boundary

Passing this checker means the oscillatory/resolvent route has explicit
fit/misfit rows, live conditional subroutes, and kill gates. It does not
construct a positive spectral measure, prove an Xi/Phi resolvent kernel,
prove reciprocal coefficient positivity, or prove `Lambda <= 0`.
