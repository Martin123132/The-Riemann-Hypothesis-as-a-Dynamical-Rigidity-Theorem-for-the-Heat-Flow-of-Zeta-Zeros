# Jensen-Window PF Signed J-Fraction Theorem Target

Date: 2026-07-06

Status: signed J-fraction theorem target specification. This is not a proof of
the all-order determinant signature, the all-order column recurrence, Schur
positivity, Jensen-window PF-infinity, Jensen hyperbolicity, RH, or
`Lambda <= 0`.

## Purpose

The reciprocal signed J-fraction scout gives a precise finite pattern:

```text
(-1)^(r(r-1)/2) Delta_r > 0
lambda_n = Delta_{n+1} Delta_{n-1} / Delta_n^2 < 0
kappa_n = -lambda_n > 0
```

The signed Jacobi beta scout adds the finite diagonal pattern:

```text
Delta_r^* = determinant with only the final Hankel column shifted
Q_r = Delta_r^*/Delta_r
beta_n = Q_{n+1} - Q_n
beta_0 > 0, beta_1 < 0, beta_n > 0 for n >= 2
```

This note states the missing theorem needed to turn that pattern into the
column recurrence:

```text
mu_m = [t^m] 1/H(-t) >= 0
```

for every degree and shift of the actual zeta heat-flow Jensen windows.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_signed_j_fraction_theorem_target.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_signed_j_fraction_theorem_target.py
```

Current result:

```text
validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows
```

## Objects

```text
h_j = binom(d,j) * A_{n+j}(0)
g_j = h_j / h_0
H(t) = 1 + g_1*t + ... + g_d*t^d
E(t) = 1 / H(-t)
mu_m = [t^m] E(t)
Delta_r = det(mu_{i+j})_{0<=i,j<r}
Delta_r^* = determinant with only the final Hankel column shifted
Q_r = Delta_r^*/Delta_r
beta_n = Q_{n+1} - Q_n
lambda_n = Delta_{n+1} Delta_{n-1} / Delta_n^2
kappa_n = -lambda_n
```

## Missing Theorem

Prove a theorem of the following kind.

Hypotheses:

```text
1. H(t) is the actual normalized zeta heat-flow Jensen-window polynomial for
   every degree d and shift n at lambda=0.
2. The reciprocal coefficients mu_m have nonzero Hankel determinants Delta_r.
3. The all-order signed signature (-1)^(r(r-1)/2) Delta_r > 0 holds, or
   equivalently the signed J-fraction parameters kappa_n=-lambda_n are
   positive for every n.
4. The ordinary diagonal Jacobi beta parameters are compatible with the
   finite beta signature, or the proof explains why that signature is not
   part of the required mechanism.
5. Any extra analytic hypotheses are stated in Xi/Phi language and do not
   assume Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya
   membership, RH, or Lambda <= 0.
```

Conclusion:

```text
mu_m >= 0 for every m, degree d, and shift n.
```

This would prove the column-shape recurrence subtarget only. It would not by
itself prove all Schur/Jacobi-Trudi shapes or close `jwpf_06`.

## Fit Rows

```text
sj_01_standard_positive_stieltjes_j_fraction:
  rejected_by_fraction_sign
  ordinary positive Stieltjes/Jacobi S-fraction or J-fraction machinery has
  the wrong first sign for E(t).

sj_02_signed_j_fraction_signature:
  open_missing_theorem
  the finite signed-Hankel, signed-lambda, and signed Jacobi beta evidence
  matches this route, but no all-order theorem is known.

sj_03_oscillatory_jacobi_matrix_after_sign_conjugation:
  live_if_conjugation_theorem_proved
  a sign-conjugated Jacobi or oscillatory matrix model could be relevant if it
  recovers E(t), explains the beta signature, and proves coefficient
  positivity. The oscillatory/resolvent fit matrix rejects ordinary entrywise
  nonnegative Jacobi powers, diagonal similarity, absolute-value majorants,
  and classical oscillatory spectral conclusions unless a positive scalar
  readout is proved.

sj_04_production_matrix_or_lattice_path_model:
  live_if_positive_model_proved
  a modified path or production-matrix model could prove positivity if the
  signed parameters become nonnegative weights after the correct transform;
  the raw ordinary Motzkin/J-fraction model, diagonal sign conjugation, and
  global length-parity sign repair are rejected, as is the absolute-value
  sign-state cover of raw Motzkin paths.

sj_05_indefinite_moment_problem:
  language_only
  indefinite moment language may interpret the signature but does not imply
  coefficient positivity by itself.

sj_06_finite_signature_extrapolation:
  finite_evidence_only
  finite signed-Hankel/J-fraction rows guide theorem search but cannot close
  the theorem.

sj_07_endpoint_real_rooted_model:
  circular_endpoint_only
  real-rooted endpoint language explains the signature but is circular as a
  proof input.
```

## Evidence

The finite signed J-fraction diagnostic is:

```text
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
```

The finite signed Jacobi beta diagnostic is:

```text
outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues
```

The raw ordinary Motzkin path obstruction is:

```text
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues
```

The global parity/sign-lift obstruction is:

```text
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues
```

The modified signed-model classification is:

```text
outputs/jensen_window_pf_modified_signed_model_target.md
work/rh_compute/results/jensen_window_pf_modified_signed_model_target.json
python work/rh_compute/scripts/check_jensen_window_pf_modified_signed_model_target.py
```

It validates:

```text
validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates
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

The oscillatory/resolvent fit matrix is:

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

The standard positive fraction rejection is:

```text
outputs/jensen_window_pf_reciprocal_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues
```

## Source Anchors

These are theorem-search anchors, not claims that their hypotheses are already
verified:

```text
https://epubs.siam.org/doi/10.1137/090781127
https://arxiv.org/pdf/2202.03793
outputs/sign_regularity_theorem_map.md
outputs/sign_regularity_theorem_fit_matrix.md
outputs/jensen_window_pf_modified_signed_model_target.md
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
outputs/jensen_window_pf_positive_readout_theorem_target.md
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
```

The first source anchors the classical ecosystem connecting structured
matrices, Hankel matrices, Jacobi-type continued fractions, total positivity,
and root localization. The second anchors the modern production-matrix and
Hankel-total-positivity ecosystem.

## Acceptance Conditions

A proposed signed J-fraction theorem can enter the proof route only if it
supplies:

```text
1. an exact all-degree/all-shift theorem statement;
2. verified hypotheses for the actual zeta heat-flow Jensen windows;
3. a proof that signed Hankel/J-fraction data implies mu_m>=0 for all m;
4. no endpoint real-rootedness, Jensen hyperbolicity, Laguerre-Polya, RH, or
   Lambda <= 0 assumption as input;
5. a clear statement whether it proves only column shapes or also supplies an
   all-Schur/Toeplitz lift.
```

## Kill Gates

Reject a proposed proof if it:

```text
1. uses only the finite d=2..8, n=0..20 signed J-fraction or signed Jacobi beta scout;
2. uses the raw ordinary Motzkin/J-fraction path expansion as a manifestly
   nonnegative production matrix;
3. repairs the raw ordinary Motzkin/J-fraction path expansion only by a global
   length-parity sign or diagonal sign conjugation;
4. uses only an absolute-value sign-state cover of raw Motzkin paths as a
   state-space doubled model;
5. cites ordinary positive Stieltjes/Jacobi fractions for E(t) despite the
   negative first parameter;
6. assumes endpoint real-rootedness, Jensen hyperbolicity, Laguerre-Polya
   membership, RH, or Lambda <= 0;
7. proves only formal existence of a signed fraction without coefficientwise
   nonnegativity of E(t);
8. claims to close jwpf_06 from column recurrence positivity without an
   all-Schur/Toeplitz lift.
```

## Boundary

Passing this checker means the signed J-fraction route has a precise missing
theorem, source anchors, fit/misfit rows, acceptance conditions, kill gates,
and a modified signed-model classification. It does not prove that theorem, the all-order column recurrence,
Jensen-window PF-infinity, or `Lambda <= 0`.
