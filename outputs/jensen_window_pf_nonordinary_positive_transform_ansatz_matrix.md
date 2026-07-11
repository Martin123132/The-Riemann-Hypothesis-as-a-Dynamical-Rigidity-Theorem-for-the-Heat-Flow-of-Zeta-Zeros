# Jensen-Window PF Nonordinary Positive Transform Ansatz Matrix

Date: 2026-07-06

Status: nonordinary positive-transform ansatz matrix. This is not a proof of
the all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The positive spectral moment obstruction rejects the ordinary interpretation:

```text
mu_m = int x^m dnu(x), with nu a positive measure
```

for the raw reciprocal coefficients of:

```text
H(t) = 1 + g_1*t + ... + g_d*t^d
E(t) = 1/H(-t)
mu_m = [t^m]E(t)
```

because:

```text
Delta_2(mu) = det [[mu_0, mu_1], [mu_1, mu_2]] = -g_2
```

and the finite signed-J-fraction grid has `Delta_2<0`.

This note asks what a surviving nonordinary positive transform would have to
be. It does not construct one.

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.py
```

Current result:

```text
validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows
```

## Acceptance Tests

A usable nonordinary positive transform must:

```text
1. prove the proposed readout has Taylor coefficients exactly
   mu_m=[t^m]1/H(-t) for every m,d,n;
2. explain why Delta_2(mu)=-g_2 does not apply, without claiming the raw
   mu_m are ordinary power moments of a positive measure;
3. identify the positive object: positive kernel, positive functional on a
   non-power basis, nonnegative transfer model, or Xi/Phi integral identity;
4. verify all hypotheses from Xi/Phi or zeta heat-flow data without endpoint
   PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0;
5. survive lambda_n<0, beta_1<0, raw Motzkin cancellation, same-length mixed
   signs, and the absolute-value sign-lift gap;
6. state all-degree, all-shift, and all-m coverage;
7. state whether the result is column-only or supplies an all-Schur/Toeplitz
   lift.
```

## Ansatz Rows

```text
npt_01_raw_power_moment_measure:
  rejected_by_delta2_obstruction
  Ordinary positive power moments for the raw mu_m are rejected by
  Delta_2(mu)=-g_2 and the 735 finite Delta_2<0 rows.

npt_02_positive_majorant_or_absolute_value_basis:
  rejected_by_exactness_gap
  Positive majorants and absolute-value covers fail because the known cover
  overshoots mu_2 by 2*kappa_1>0.

npt_03_fixed_triangular_basis_change:
  conditional_wrapper_only
  A window-independent triangular basis change is only a wrapper unless its
  inverse extraction proves the original mu_m>=0.

npt_04_nonpower_positive_functional:
  live_if_functional_constructed
  A positive functional L on a non-power basis K_m with L(K_m)=mu_m remains
  live, but the low-degree scout shows it must absorb the signed mu_2 and
  mu_3 cancellation identities; no cone C, K_m, or L is currently
  constructed.

npt_05_xi_phi_positive_kernel_identity:
  live_if_kernel_identity_proved
  An Xi/Phi-specific positive kernel identity remains live if it is exactly
  E(t)=1/H(-t) and is noncircular.

npt_06_modified_state_space_transfer_model:
  live_if_exact_model_constructed
  A genuinely modified nonnegative transfer model remains live, but the cheap
  absolute-value and parity-lift repairs are already rejected.

npt_07_finite_quadrature_or_prefix_fit:
  finite_evidence_only
  Finite positive fits can suggest kernels but do not prove all m,d,n.

npt_08_endpoint_or_lp_factorization_basis:
  circular_if_endpoint_used
  Endpoint real-rootedness or Laguerre-Polya factorization cannot be used as
  input to prove the missing readout.
```

## Evidence Anchors

Parent target and obstruction:

```text
outputs/jensen_window_pf_positive_readout_theorem_target.md
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json
python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py
```

They validate:

```text
validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
```

Signed and modified-model obstructions:

```text
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
outputs/jensen_window_pf_modified_signed_model_target.md
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md
outputs/signed_hankel_jensen_dependency_graph.md
```

The nonpower-functional low-degree scout validates:

```text
validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows
```

## Kill Gates

Reject a proposed transform if it:

```text
1. claims mu_m are ordinary power moments of a positive measure;
2. uses a positive majorant, absolute-value cover, fitted quadrature, or
   finite prefix without exact equality to E(t);
3. proves positivity only after a transform and never proves positive inverse
   extraction for the original mu_m;
4. assumes endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH,
   or Lambda <= 0;
5. cites generic spectral, Pick, Hamburger, or kernel machinery without
   constructing the positive object and exact identity;
6. advertises a column-only theorem as an all-Schur/Toeplitz lift.
```

## Boundary

Passing this checker means the nonordinary positive-transform search has a
sharper ansatz matrix and explicit rejection gates. It does not construct a
positive functional, prove an Xi/Phi positive kernel, prove an exact modified
state-space transfer model, prove reciprocal coefficient positivity, or prove
`Lambda <= 0`.
