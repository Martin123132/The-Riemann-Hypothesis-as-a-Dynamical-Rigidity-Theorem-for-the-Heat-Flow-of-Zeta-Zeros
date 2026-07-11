# Jensen-Window PF Positive Readout Theorem Target

Date: 2026-07-06

Status: positive-readout theorem target. This is not a proof of the all-order
column recurrence, Schur positivity, Jensen-window PF-infinity, Jensen
hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The oscillatory/resolvent fit matrix leaves only two live subroutes:

```text
or_06_positive_spectral_measure_after_transform
or_07_xi_phi_resolvent_kernel
```

This note turns that narrowing into a sharper theorem target. For:

```text
H(t) = 1 + g_1*t + ... + g_d*t^d
E(t) = 1 / H(-t)
mu_m = [t^m] E(t)
```

the missing positive-readout theorem must prove:

```text
mu_m >= 0
```

for every `m`, degree `d`, and shift `n` by an exact positive scalar readout.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py
```

Current result:

```text
validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes
```

## Theorem Obligation

A usable positive readout must:

```text
1. be exactly E(t)=1/H(-t), not a majorant, approximation, interpolation, or finite fit;
2. prove mu_m>=0 for every m by a positive measure, positive kernel,
   nonnegative path model, or equivalent positive functional;
3. verify all hypotheses from Xi/Phi or zeta heat-flow data without endpoint
   PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0;
4. survive lambda_n<0, beta_1<0, raw Motzkin cancellation, same-length mixed
   signs, and the absolute-value sign-lift gap;
5. state whether it is column-only or supplies an all-Schur/Toeplitz lift;
6. state the all-degree, all-shift, and all-m coverage.
```

## Candidate Rows

```text
pr_01_positive_spectral_transform:
  live_if_transform_constructed
  A nonordinary noncircular transform from signed Jacobi data to a positive
  scalar readout would work if its Taylor coefficients are exactly mu_m and
  all hypotheses are verified. The ordinary positive moment interpretation is
  rejected by Delta_2=-g_2, and the nonordinary ansatz matrix narrows the
  surviving forms.

pr_02_xi_phi_positive_resolvent_kernel:
  live_if_kernel_identity_proved
  An Xi/Phi-specific positive kernel remains live if its scalar resolvent is
  exactly E(t), it proves coefficient positivity, and it passes the
  nonordinary ansatz acceptance tests.

pr_03_abstract_stieltjes_pick_or_hamburger_wrapper:
  conditional_wrapper_only
  Classical moment or analytic-function machinery can package a proof only
  after the exact positive readout exists.

pr_04_indefinite_signed_spectral_readout:
  language_only
  Indefinite spectral language describes the signature but does not imply
  coefficient positivity without a positive scalar functional.

pr_05_endpoint_root_factorization_readout:
  circular_if_endpoint_used
  Product factorization from real nonpositive roots would prove positivity,
  but as input it assumes the endpoint target.

pr_06_finite_quadrature_or_interpolation_fit:
  finite_evidence_only
  Finite fitted measures or quadrature patterns can guide theorem search but
  do not prove all m,d,n.

pr_07_raw_signed_jacobi_scalar_resolvent:
  rejected_by_signed_readout
  The raw signed resolvent is not a positive readout because lambda_n<0,
  beta_1<0, mu_2 cancellation, and same-length mixed signs remain.

pr_08_absolute_value_or_majorant_readout:
  rejected_by_exactness_gap
  Absolute-value path weights are positive but not exact; at mu_2 the gap is
  2*kappa_1>0.
```

## Evidence Anchors

The parent oscillatory/resolvent fit matrix is:

```text
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py
```

It validates:

```text
validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows
```

The signed-fraction theorem target and route matrix are:

```text
outputs/jensen_window_pf_signed_j_fraction_theorem_target.md
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
outputs/jensen_window_pf_modified_signed_model_target.md
```

They validate:

```text
validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates
```

The signed-data and obstruction anchors are:

```text
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md
```

They validate:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues
validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues
validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues
validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows
```

Finite positive rows remain evidence only:

```text
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
outputs/arb_jensen_window_column_recurrence_stress.md
```

## Kill Gates

Reject a proposed readout if it:

```text
0. claims mu_m are ordinary power moments of a positive measure;
1. is a positive majorant, interpolation, fitted quadrature, or finite prefix
   rather than exactly E(t)=1/H(-t);
2. infers positive spectral measure from endpoint real-rootedness,
   Jensen-window PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0;
3. gives only signed or indefinite spectral language without a positive scalar
   functional proving mu_m>=0;
4. ignores lambda_n<0, beta_1<0, same-length mixed signs, or the
   absolute-value sign-lift gap;
5. cites a theorem wrapper without verifying its hypotheses for the actual
   zeta windows;
6. claims to close jwpf_06 from a column-only readout without an
   all-Schur/Toeplitz lift.
```

## Boundary

Passing this checker means the positive-readout search has explicit live
subroutes, exactness requirements, and rejection gates. It does not construct
a positive spectral transform, prove an Xi/Phi positive kernel, prove
reciprocal coefficient positivity, or prove `Lambda <= 0`.
