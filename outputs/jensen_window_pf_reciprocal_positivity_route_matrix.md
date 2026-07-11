# Jensen-Window PF Reciprocal Positivity Route Matrix

Date: 2026-07-06

Status: reciprocal-positivity theorem-search matrix. This is not a proof of Schur positivity,
Jensen-window PF-infinity, Jensen hyperbolicity,
Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Purpose

The column recurrence contract reduces the hard column-shape frontier to
coefficient positivity of a reciprocal polynomial:

```text
h_j = binom(d,j) A_{n+j}
g_j = h_j / h_0
H(t) = 1 + g_1*t + ... + g_d*t^d
E(t) = 1 / H(-t)
```

The column target is:

```text
[t^m] E(t) >= 0
for every m, degree d, and shift n.
```

This note asks which theorem mechanisms could prove that target without
assuming Jensen-window PF-infinity or using finite recurrence grids as an
all-order theorem.

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
```

Current result:

```text
validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows
```

## Route Rows

```text
rp_01_real_negative_zero_or_pf_window_endpoint:
  endpoint_equivalence_only
  If H(t) already has only real nonpositive zeros, then 1/H(-t) has
  nonnegative coefficients.  This is endpoint language, not an independent
  proof input for jwpf_06.

rp_02_positive_renewal_or_convolution_kernel:
  live_if_representation_proved
  A nonnegative renewal, convolution, or integral representation for
  [t^m]1/H(-t) would prove the column recurrence.  No Xi/Phi-specific
  representation is currently known.

rp_03_positive_stieltjes_or_j_fraction:
  route_mismatch
  The standard positive Stieltjes/Jacobi S-fraction or J-fraction has the
  wrong first nontrivial sign for E(t)=1/H(-t): a_2 = -g_2/g_1 and
  lambda_1 = -g_2.

rp_04_companion_or_production_matrix_total_positivity:
  live_if_representation_proved
  A nonnegative or totally nonnegative production matrix would turn the
  recurrence coefficients into path sums.  The obvious companion recurrence
  has alternating signs, the raw ordinary Motzkin/J-fraction path model has
  negative weights, global length-parity signs cannot repair same-length
  mixed signs, the absolute-value sign-state cover overshoots mu_2 by
  2*kappa_1, and no positive modified conjugate/network is known.

rp_05_kaluza_reciprocal_sign_theorem:
  route_mismatch
  Kaluza-type reciprocal sign theorems do not currently match the
  sign-alternating denominator H(-t) and the available Jensen-window
  hypotheses.

rp_06_generic_ratio_or_log_concavity_conditions:
  rejected_by_countermodel
  Generic log-concavity, ratio monotonicity, and finite low-degree
  positivity have already been rejected as standalone bridge conditions.

rp_07_finite_recurrence_stress_grid:
  finite_evidence_only
  The 12,600 positive Arb recurrence stress rows and 72,600 extended
  normalized reciprocal coefficient rows are useful theorem-search evidence
  but cannot become an all-order recurrence theorem by wording.

rp_08_full_schur_or_toeplitz_lift:
  blocked_until_input_pf
  Column recurrence positivity is necessary but not sufficient for all
  Schur/Jacobi-Trudi shapes.  A positive specialization, planar network, or
  full Toeplitz total-positivity lift is still required.

rp_09_signed_or_modified_continued_fraction:
  live_if_representation_proved
  A signed or modified continued fraction remains possible only if it comes
  with a separate oscillatory, sign-regular, or production-matrix theorem
  compatible with the negative ordinary fraction parameter. The modified
  signed-model audit now separates dead raw sign repairs from the four
  conditionally live modified classes, and the oscillatory/resolvent fit
  matrix leaves only positive spectral-transform or Xi/Phi positive-kernel
  variants live inside the oscillatory subroute. The positive-readout theorem
  target states the exact scalar-readout obligations for those variants.
  The ordinary positive moment-measure version is blocked by Delta_2=-g_2.
```

## Fraction Sign Scout

The continued-fraction row is sharpened by:

```text
outputs/jensen_window_pf_reciprocal_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_fraction_sign_lamgrid_n0_n20_d2_d8_dps520.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues
```

The signed J-fraction row is sharpened by:

```text
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_lamgrid_n0_n20_d2_d8_dps520.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
```

The signed Jacobi beta row is sharpened by:

```text
outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_lamgrid_n0_n20_d2_d8_dps520.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues
```

The raw ordinary Motzkin path model is blocked by:

```text
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues
```

The global parity/sign-lift repair is blocked by:

```text
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_lamgrid_n0_n20_d2_d8_m2_m8_dps520.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py
```

It validates:

```text
validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues
```

The precise missing theorem target is:

```text
outputs/jensen_window_pf_signed_j_fraction_theorem_target.md
work/rh_compute/results/jensen_window_pf_signed_j_fraction_theorem_target.json
python work/rh_compute/scripts/check_jensen_window_pf_signed_j_fraction_theorem_target.py
```

It validates:

```text
validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows
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

The oscillatory/resolvent route is sharpened by:

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

The extended direct reciprocal coefficient stress is:

```text
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json
work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_lamgrid_n0_n32_d2_d12_m1_m40_dps520.jsonl
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_coefficient_extended_stress.py
```

It validates:

```text
validated Jensen-window PF reciprocal coefficient extended stress: 72600 rows, 0 issues
```

## Integration

This route matrix sharpens:

```text
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
outputs/arb_jensen_window_column_recurrence_stress.md
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_contraction_log_concavity_scout.md
outputs/jensen_window_pf_schur_shape_contract.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
```

The main conclusion is deliberately modest: No row is `ready_to_apply`, and
no row closes `jwpf_06`.  The live routes are now narrowed to a positive
renewal/convolution representation, a signed or modified continued fraction,
or a positive production-matrix/network representation for the actual zeta
heat-flow Jensen windows.  The standard positive Stieltjes/Jacobi continued
fraction route is not compatible with the first symbolic fraction signs. The
signed J-fraction and signed Jacobi beta signatures have strong finite
support, but the raw ordinary Motzkin path model is already blocked by
cancellation, and a global length-parity sign cannot repair the same-length
mixed path signs.  The absolute-value sign-state cover is also blocked at
mu_2 by the gap 2*kappa_1.  A surviving signed or production-matrix route must
therefore be genuinely modified: a state-space doubled model, a modified
production/Riordan/network model, an oscillatory sign-regular resolvent
theorem, or an Xi/Phi-specific positive cancellation identity.  All four still
lack the all-order theorem or construction that would make them proof routes.

The finite stress evidence is now stronger:

```text
12,600 positive unnormalized column recurrence rows
72,600 positive normalized reciprocal coefficient rows
```

but finite grids still cannot close the theorem.

## Boundary

Passing this checker means the reciprocal-positivity routes have been
classified without promoting endpoint equivalences, generic ratio conditions,
or finite stress grids into a proof.  It does not prove the all-order column
recurrence, Schur positivity, Jensen-window PF-infinity, or `Lambda <= 0`.
