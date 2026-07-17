# Countermodel Library

Date: 2026-07-09

Status: proof-safety artifact. This is not evidence against RH. It records small models, finite-prefix traps, and finite-grid traps that block invalid bridge lemmas in the RH dynamical-rigidity programme.

Executable gate:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

Current result:

```text
validated 11 countermodel gate examples
```

## Purpose

The programme is trying to prove the missing Newman direction:

```text
Lambda <= 0.
```

The dangerous failure mode is to turn one of these into a proof:

```text
local zero repulsion
finite signed-Hankel evidence
finite Jensen hyperbolicity
finite Jensen-window PF/Sturm rectangles
finite Toeplitz/PF certificates
finite Schur/Toeplitz shape-prefix evidence
finite Edrei moment-recurrence evidence
Stieltjes/Hankel moment positivity without Toeplitz total positivity
```

Each item is useful, but none is logically enough by itself. This file gives concrete gates that proposed proof steps must pass before they can enter a manuscript as more than evidence or a conditional claim.

## Gate 1: Local Heat Birth

The exact model:

```text
F_tau(t) = t^2 - 2 tau
```

satisfies the same Newman heat sign:

```text
partial_tau F = -2
partial_tt F = 2
partial_tau F = - partial_tt F
```

Its zeros are:

```text
tau < 0:  t = +/- i sqrt(-2 tau)
tau = 0:  double real zero at t = 0
tau > 0:  t = +/- sqrt(2 tau)
```

On the real-zero side, the gap is:

```text
g(tau) = 2 sqrt(2 tau)
```

and:

```text
g'(tau) = 4/g(tau).
```

Blocked proof step:

```text
The local law g' = 4/g excludes positive Newman birth.
```

Correct use:

```text
The local law describes the already-real side after a square-root birth.
```

Therefore any proof of `Lambda <= 0` must use global Xi structure, a sign-regularity theorem, or another nonlocal invariant. Local repulsion alone cannot do the job.

## Gate 2: Finite Prefix Is Not All-Order PF

Let:

```text
c_k(lambda) = A_k(lambda) / k!
```

be the ordinary coefficient sequence used in the coefficient-PF route. Suppose a finite prefix through `c_K` has passed every certified Toeplitz/PF test we have run.

That still cannot prove PF-infinity. Preserve the whole prefix through `K`, and define one positive next coefficient by:

```text
c_{K+1} = 2 c_K^2 / c_{K-1}.
```

Then the order-2 Toeplitz minor:

```text
det [[c_K,   c_{K+1}],
     [c_{K-1}, c_K]]
= c_K^2 - c_{K-1} c_{K+1}
= -c_K^2
< 0.
```

This does not say the actual zeta coefficient sequence fails. It says a proof step using only a finite prefix is invalid unless it adds a structural all-order theorem.

Blocked proof step:

```text
The certified finite Toeplitz/PF ledger makes c_k PF-infinity.
```

Correct use:

```text
The finite Toeplitz/PF ledger is falsification pressure and theorem-search evidence.
```

## Gate 2A: Finite Schur Prefix Is Not A Positive Specialization

The positive Schur-specialization target studies:

```text
h_k -> d_k(0)
```

and tries to prove all skew-Schur evaluations nonnegative. A finite set of
Schur or Toeplitz checks still cannot prove such a specialization exists.

The exact gate starts with:

```text
h_k = 1/k!
```

whose generating function is `exp(z)`, a clean restricted PF-infinity model.
It preserves:

```text
h_0, h_1, ..., h_6
```

and exactly checks a finite Toeplitz/Schur grid:

```text
N = 7
orders <= 4
2,940 finite tests
1,204 structurally nonzero positive minors
1,736 structural zeros
```

Then it chooses the next positive complete-homogeneous coordinate:

```text
h_7 = 1/2160
```

At the first untested Jacobi-Trudi shape:

```text
lambda = (6,6)
mu = (0,0)
```

the determinant becomes:

```text
s_(6,6) = det [[h_6, h_7],
               [h_5, h_6]]
        = -1/518400 < 0.
```

Blocked proof step:

```text
A finite Schur/Toeplitz shape ledger proves h_k -> d_k is a positive
specialization.
```

Correct use:

```text
Finite Schur/Toeplitz checks test proposed formulas. A proof needs an
all-order positive specialization, planar network, production matrix,
continued fraction, positive determinant integral, or equivalent theorem.
```

## Gate 3: Finite Signed-Hankel Is Not All-Order Signed Regularity

For:

```text
D_{m,s} = det(A_{i+j+s})_{i,j=0}^m
sigma(m) = (-1)^(m(m+1)/2)
```

the observed signed-Hankel condition is:

```text
sigma(m) D_{m,s} > 0.
```

For `m = 1`, this says:

```text
-(A_s A_{s+2} - A_{s+1}^2) > 0.
```

Given any positive prefix through `A_K`, preserve it and define:

```text
A_{K+1} = 2 A_K^2 / A_{K-1}.
```

At shift `s = K-1`:

```text
D_{1,K-1} = A_{K-1} A_{K+1} - A_K^2
          = A_K^2
          > 0
```

so:

```text
sigma(1) D_{1,K-1} = -A_K^2 < 0.
```

Blocked proof step:

```text
The finite signed-Hankel certificate grid proves all-order signed regularity.
```

Correct use:

```text
The finite signed-Hankel grid supports a conjectural all-order signed-regularity target.
```

## Gate 4: Finite Jensen Hyperbolicity Is Not Jensen Criterion

For degree 2 and shift `K-1`, the Jensen polynomial has the form:

```text
P_{2,K-1}(x)
= A_{K-1} + 2 A_K x + A_{K+1} x^2.
```

With the same positive one-term extension:

```text
A_{K+1} = 2 A_K^2 / A_{K-1}
```

the discriminant becomes:

```text
4(A_K^2 - A_{K-1} A_{K+1})
= -4 A_K^2
< 0.
```

So the next degree-2 Jensen polynomial is not hyperbolic, despite every earlier finite check being preserved.

Blocked proof step:

```text
Many finite Jensen/Sturm passes imply all Jensen polynomials are hyperbolic.
```

Correct use:

```text
Finite Jensen/Sturm passes guide theorem search and catch numerical failures.
```

## Gate 4A: Finite Consecutive Signed-Hankel Grid Is Not All Shifts

The executable gate also includes an exact finite-grid trap independent of
the zeta coefficients. Start with:

```text
a_k = 1/k!
```

For the shifted-principal signed-Hankel determinants:

```text
sigma(m) det(a_{i+j+s})_{i,j=0}^m
sigma(m) = (-1)^(m(m+1)/2)
```

the exact rational check validates the whole grid:

```text
m = 0..4
s = 0..8
45/45 signed determinants > 0
coefficients used: a_0..a_16
minimum signed grid value about 4.048066611965355946e-53
```

Now preserve all those coefficients and choose the next positive coefficient:

```text
a_17 = 1/167382319104000.
```

At the next untested shift:

```text
s = 15
```

the `m = 1` signed-Hankel value becomes negative:

```text
-2.284340357080483672e-27
```

and the degree-2 Jensen discriminant at the same shift is negative:

```text
-9.137361428321934688e-27.
```

Blocked proof step:

```text
The finite shifted-principal signed-Hankel grid proves all shifts, all-order
sign-regularity, or Jensen hyperbolicity.
```

Correct use:

```text
The finite grid is a certified finite diagnostic. A promotion needs a known
Hankel sign-consistency reduction plus an all-order proof, or a new bridge
theorem matching the zeta coefficient sequence.
```

## Gate 4B: Current Jensen-Window Rectangle Is Not All Shifts

The current Arb Jensen-window diagnostics are deliberately finite:

```text
PF obligation checks:
  degrees d = 3,4
  shifts n = 0..20
  coefficients used through A_24

Sturm/root-count checks:
  degrees d = 3,4,5
  shifts n = 0..20
  coefficients used through A_25
```

The executable gate preserves every coefficient those finite Jensen-window
checks can see:

```text
A_0, A_1, ..., A_25
```

for all five lambda rows, then chooses a positive next coefficient:

```text
A_26 = 2 A_25^2 / A_24.
```

At the next untested degree-2 Jensen window, shift `24`, the discriminant is:

```text
4(A_25^2 - A_24 A_26) = -4 A_25^2 < 0.
```

So all existing finite Jensen-window PF/Sturm inputs are preserved, but the
next degree-2 Jensen polynomial is not hyperbolic.

Blocked proof step:

```text
The finite Jensen-window PF/Sturm rectangle proves all-shift Jensen
hyperbolicity.
```

Correct use:

```text
The finite Jensen-window rectangle is a strong stress test and normalization
check. A proof still needs an all-degree/all-shift theorem for the actual
coefficient sequence.
```

## Gate 5: Finite Moment Recurrence Is Not An Edrei Representation

The Edrei reconstruction scout checks finite recurrence data for the shifted
moment sequence:

```text
a_n = p_{n+1}.
```

Positive recurrence data through any fixed order still does not prove an
all-order positive measure or Edrei zero-parameter representation.

The executable gate uses exact rational arithmetic. Start with the factorial
moments:

```text
m_n = n!
```

These are genuine Stieltjes moments of `exp(-x) dx` on `[0, infinity)`, so the
finite recurrence/Hankel prefix is honestly positive. For the current default
order `12`, preserve:

```text
m_0, m_1, ..., m_23
```

Then choose the next even moment `m_24` as a positive number below the exact
Schur-complement threshold for the next Hankel matrix. The gate reports:

```text
all preserved leading Hankel determinants > 0
adversarial m_24 > 0
next Hankel determinant < 0
```

Blocked proof step:

```text
The finite Arb recurrence scout through order 12 proves the all-order Edrei
moment representation.
```

Correct use:

```text
The recurrence scout is a constructive finite diagnostic and precision
frontier. It becomes a proof only after an all-order moment theorem,
positive parameter construction, or analytic recurrence formula is supplied.
```

## Gate 6: Stieltjes Moment Positivity Is Not Coefficient PF

The coefficient-PF route uses:

```text
c_k = mu_k / (2k)!.
```

A tempting but invalid bridge is:

```text
mu_k is a Stieltjes moment sequence
and 1/(2k)! has a restricted Laguerre-Polya generating function
therefore c_k is PF-infinity.
```

The executable gate blocks this with an exact positive measure:

```text
10 delta_0 + delta_1 + delta_2 on [0, infinity).
```

Its moments begin:

```text
mu_0..mu_6 = 12, 3, 5, 9, 17, 33, 65
```

and the leading Hankel determinants are:

```text
size 1..4 = 12, 51, 40, 0
```

so the moment sequence is Stieltjes/Hankel-nonnegative. But after the
coefficient-route normalization:

```text
c_0 = 12
c_1 = 3/2
c_2 = 5/24
```

the first order-2 Toeplitz/PF minor is:

```text
c_1^2 - c_0 c_2 = -1/4 < 0.
```

Blocked proof step:

```text
Stieltjes/Hankel positivity of the moments mu_k, together with the
factor 1/(2k)!, proves coefficient PF-infinity.
```

Correct use:

```text
Moment positivity can motivate the route, but a valid proof needs a
Toeplitz-total-positivity theorem, a positive determinant integral formula,
or an explicit restricted Laguerre-Polya factorization.
```

## Current Cache-Based Gate Run

The executable gate reads:

```text
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
```

and applies the finite-prefix traps to the five lambda rows:

```text
0
1e-6
1e-4
1e-2
1e-1
```

With the current cache, the prefix is preserved through:

```text
K = 32
```

For each lambda, a positive one-term extension at `K+1` breaks:

```text
order-2 Toeplitz/PF
m = 1 signed-Hankel
degree-2 Jensen hyperbolicity
```

The same run also validates the exact rational moment-recurrence trap:

```text
recurrence order <= 12
moments 0..23 preserved
positive edited moment 24
next Hankel/moment gate breaks
```

It validates the exact finite shifted-principal signed-Hankel grid trap:

```text
base a_k = 1/k!
m <= 4, shifts <= 8 preserved
positive a_17 breaks the next shifted m = 1 signed-Hankel/Jensen gate
```

It validates the exact finite Schur-prefix trap:

```text
base h_k = 1/k!
h_0..h_6 preserved
2,940 finite Toeplitz/Schur tests preserved
positive h_7 breaks s_(6,6)
```

It validates the finite Jensen-window rectangle trap:

```text
current Jensen-window PF/Sturm coefficient inputs A_0..A_25 preserved
positive A_26 breaks degree-2, shift-24 Jensen hyperbolicity
```

It also validates the Stieltjes multiplier trap:

```text
positive measure = 10 delta_0 + delta_1 + delta_2
leading Hankel determinants = 12, 51, 40, 0
c_1^2 - c_0 c_2 = -1/4
```

Again, these are proof-safety models, not claims about the actual zeta coefficients. Their role is logical: they show that finite-prefix evidence, finite shifted-principal grids, finite recurrence evidence, and generic Stieltjes/Hankel positivity cannot be promoted into an all-order theorem by wording alone.

## Manuscript Rule

Additional signed-model gate:

```text
The raw ordinary Motzkin/J-fraction path model, diagonal sign conjugation,
global length-parity sign repairs, and absolute-value sign-state covers cannot
be used as manifest positivity proofs for E(t)=1/H(-t). A surviving signed
route must be a genuinely modified state-space, production/Riordan/network,
oscillatory resolvent, or Xi/Phi positive-cancellation construction, and it
must still prove all-order coefficientwise nonnegativity.
```

Before accepting any proposed bridge lemma, ask:

```text
Does the lemma fail on the local heat-birth model?
Does the lemma rely only on a finite coefficient prefix?
Does the lemma rely only on a finite Schur/Toeplitz shape prefix?
Does the lemma rely only on a finite shifted-principal signed-Hankel grid?
Does the lemma rely only on a finite moment or recurrence prefix?
Does the lemma prove only Stieltjes/Hankel positivity when Toeplitz PF is required?
Does the lemma silently assume PF-infinity, Jensen hyperbolicity, or Laguerre-Polya membership?
Does the lemma assume RH at lambda = 0?
```

If yes, it cannot be used to prove `Lambda <= 0`.

The only acceptable promotions are:

```text
finite certificate
conditional theorem with explicit hypotheses
all-order theorem with noncircular structural proof
```

Executable result-language scan:

```text
python work/rh_compute/scripts/check_output_reference_integrity.py
python work/rh_compute/scripts/check_output_status_manifest.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py
python work/rh_compute/scripts/check_arb_jensen_window_column_recurrence_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_signed_j_fraction_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_modified_signed_model_target.py
python work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_cone_candidate_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_cone_frontier_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_column_extension_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree6_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_all_m_counterexample.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_monotone_closure_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_boundary_threshold_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_cone_entry_asymptotic_target.py
python work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_log_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_uniform_remainder_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_state_space_sign_lift_obstruction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py
python work/rh_compute/scripts/check_result_language_boundaries.py
```

Current result:

```text
validated output references: scanned 191 markdown files, 2523 path references, 0 missing required paths, 3 planned missing deliverables
validated output artifact statuses: scanned 191 markdown files, 0 status issues
validated proof-claim ledger: 149 claims, 0 issues, 14 open theorem targets
validated signed-Hankel/Jensen dependency graph with 0 issues
validated Jensen-window PF bridge obligations: 11 obligations, 0 issues, 3 open obligations
validated Jensen-window PF theorem machinery fit matrix: 7 rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF sign-regular transfer gap matrix: 9 transfer rows, 2 countermodel gates, 3 open requirements, 3 rejected shortcuts, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF factorial multiplier split audit: 5 exact rows, 315 raw degree-2 anti-hyperbolic rows, 315 normalized degree-2 positive rows, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows
validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows
validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues
validated 12600 Arb Jensen-window column recurrence stress rows with 0 issues
validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates
validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows
validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes
validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues
validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows
validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows
validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows
validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows
validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues
validated Jensen-window PF monotone-contraction column extension scout: 25 column rows, 3329 Bernstein coefficients, 3 beyond-frontier rows, 0 negative Bernstein rows, 0 issues
validated Jensen-window PF monotone-contraction sparse degree-6 scout: 10 degree-6 rows, 63347 Bernstein coefficients, m<=10, 0 negative Bernstein rows, 0 zero Bernstein rows, 0 issues
validated Jensen-window PF monotone-contraction sparse degree-7 frontier scout: 9 positive rows, 1 certificate-obstruction row, 932691 Bernstein coefficients, first obstruction m=10, 126 negative Bernstein coefficients, 0 zero Bernstein coefficients, 0 issues
validated Jensen-window PF monotone-contraction sparse degree-7 subdivision scout: 3 dyadic slabs, 785400 slab Bernstein coefficients, 0 negative slab coefficients, 0 zero slab coefficients, repaired m=10 obstruction, 0 issues
validated Jensen-window PF monotone-contraction all-m counterexample: degree 7, m=11, exact monotone witness, negative normalized value, 0 issues
validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
validated Jensen-window PF kernel Mellin upper-wall certificate: 8 rows, 0 issues, 200 positive compact intervals, 1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows
validated Jensen-window PF log-concave Mellin monotone-wall countermodel: 6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation
validated Jensen-window PF T=1156 monotone-wall counterexample certificate: 7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation
validated Jensen-window PF negative-lambda kernel summand-shift lemma: 8 rows, 0 issues, 6 exact rows, 1 compact interval row, 1 open far-tail row, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand dominance certificate: 10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: 6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand saddle-wall target: 9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, 1 open requirement, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand cumulant bridge: 8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, 1 open requirement, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: 8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 0 ready-to-apply rows, 1 live routes
validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 69 coefficient rows, 63 lower-wall rows, 63 upper-wall rows, 60 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 183 coefficient rows, 177 lower-wall rows, 177 upper-wall rows, 174 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 243 coefficient rows, 237 lower-wall rows, 237 upper-wall rows, 234 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 303 coefficient rows, 297 lower-wall rows, 297 upper-wall rows, 294 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 453 coefficient rows, 447 lower-wall rows, 447 upper-wall rows, 444 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda cone-entry prefix scout: 603 coefficient rows, 597 lower-wall rows, 597 upper-wall rows, 594 monotone-gap rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=19, 57 active lower rows, 57 active upper rows, 57 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=57, 171 active lower rows, 171 active upper rows, 171 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=77, 231 active lower rows, 231 active upper rows, 231 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=97, 291 active lower rows, 291 active upper rows, 291 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=147, 441 active lower rows, 441 active upper rows, 441 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=197, 591 active lower rows, 591 active upper rows, 591 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 63 cone-buffer rows, 60 defect-monotone rows, 63 one-third-buffer rows, 60 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 177 cone-buffer rows, 174 defect-monotone rows, 139 one-third-buffer rows, 174 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 177 scaled rows, 177 cone rows, 177 half-width rows, 139 one-third rows, 38 one-third failures, 174 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 237 cone-buffer rows, 234 defect-monotone rows, 159 one-third-buffer rows, 234 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 237 scaled rows, 237 cone rows, 237 half-width rows, 159 one-third rows, 78 one-third failures, 234 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 297 cone-buffer rows, 294 defect-monotone rows, 179 one-third-buffer rows, 294 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 297 scaled rows, 297 cone rows, 297 half-width rows, 179 one-third rows, 118 one-third failures, 294 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 447 cone-buffer rows, 444 defect-monotone rows, 179 one-third-buffer rows, 444 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda tail-barrier scout: 597 cone-buffer rows, 594 defect-monotone rows, 179 one-third-buffer rows, 594 scaled-defect increase rows, 1 rejected candidate, 0 issues
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 597 scaled rows, 597 cone rows, 521 half-width rows, 179 one-third rows, 418 one-third failures, 594 scaled-increase rows, 0 issues
validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues
validated Jensen-window PF negative-lambda log-curvature bridge: 63 simple log-buffer rows, 63 exact defect-buffer rows, 60 curvature-monotone rows, 5 bridge rows, 0 issues
validated Jensen-window PF negative-lambda bounded log-curvature target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes, 63 raw-threshold rows
validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda Gaussian curvature matrix: 7 matrix rows, 63 positive-deficit rows, 63 bounded-deficit rows, 63 raw-threshold rows, 0 issues
validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: 8 matrix rows, 2 certified Taylor signs, 1 fixed-k activation estimates, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda uniform remainder target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 open requirements, 3 leading-scale rows
validated Jensen-window PF negative-lambda Taylor moment budget: 9 budget rows, 7 tail-start samples, 4 invalid truncation rows, 2 bounded truncation rows, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda high-order Taylor scout: 8 coefficient rows, 35 truncation rows, 9 invalid normalizers, 2 upper-wall violations, 3 overbound rows, 0 ready-to-apply rows, 0 issues
validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
validated Jensen-window PF negative-lambda half-width tail target: 9 rows, 0 issues, 0 ready-to-apply rows, 0 live routes, 430 half-width rows, 17 half-width failures
validated Jensen-window PF negative-lambda adaptive scaled-defect target: 8 rows, 0 issues, 2 live routes, 597 exact-cone rows, 76 half-width failures
validated Jensen-window PF negative-lambda adaptive envelope matrix: 7 matrix rows, 0 issues, 594 k-increase rows, 398 lambda-order rows, 76 half-width failures
validated Jensen-window PF negative-lambda adaptive envelope obligations: 9 obligation rows, 0 issues, 3 exact rows, 3 open requirements, 1 rejected routes
validated Jensen-window PF negative-lambda raw-moment bridge matrix: 8 matrix rows, 0 issues, 597 raw-cone rows, 594 corridor rows, 76 half-width failures
validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: 9 rows, 0 issues, 594 decrement-corridor rows, 591 theta-k-monotone rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda k300 precision-repair audit: 7 rows, 0 issues, 894 repaired decrement-corridor rows, 891 repaired theta-k-monotone rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda raw-log decrement bridge: 8 rows, 0 issues, 894 log-corridor rows, 894 log-decrease rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: 9 rows, 0 issues, 894 curvature-corridor rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda linear curvature-barrier scout: 8 rows, 0 issues, 894 linear-barrier rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: 10 rows, 0 issues, 2 live routes, 894 scaled-curvature increase rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: 8 rows, 0 issues, 894 scaled-ceiling rows, 894 scaled-log-corridor rows, 894 ceiling-dominance rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: 8 rows, 0 issues, 897 B-positive rows, 894 B-decrease rows, 894 C-increase rows, 598 C-lambda-order rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: 8 rows, 0 issues, 3 certified leading-sign rows, 35 truncation rows, 4 all-positive stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: 9 rows, 0 issues, 4 positive baseline rows, 31 blocked baseline rows, 4 exact stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: 9 rows, 0 issues, 4 positive baseline rows, 4 budget rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: 8 rows, 0 issues, 2 tested next-increment rows, 2 pointwise budget failures, 2 stencil-sign-preserving rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: 7 rows, 0 issues, 4 tested continuation rows, 3 stencil-sign-preserving rows, 1 stencil-sign-failure rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: 7 rows, 0 issues, 1301 scan rows, 1045 continuation-positive rows, 718 half-safety rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified surrogate stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified stencil rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: 8 rows, 0 issues, 13 degree levels, max degree 40, 0 failed Bernstein rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: 8 rows, 0 issues, 5 budget inequalities, 4 finite tail profile rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: 8 rows, 0 issues, 4 profile rows, 4 formal-tail turnaround rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: 6 rows, 0 issues, 4 first-omitted rows, 4 optimized-window rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout: 6 rows, 0 issues, 4 endpoint rows, 5 quadrature orders, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: 6 rows, 0 issues, 20 grid rows, 5 T values, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: 6 rows, 0 issues, 8 obligations, 5 open requirements, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: 6 rows, 0 issues, 3 tail bounds below 1e-1000, 2 conditional requirements, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: 5 rows, 0 issues, 16 Laguerre bound rows, 2 certified side conditions, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: 6 rows, 0 issues, 3 certified tail sources, 2 certified side conditions, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: 5 rows, 0 issues, 7 ladder rows, 320 reference order, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: 7 rows, 0 issues, derivative order 640, 2 derivative-sup caps, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: 7 rows, 0 issues, split y=200, 2 tail ratios below quadrature cap, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout: 7 rows, 0 issues, 6 panels, plain interval Riemann rejected, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout: 7 rows, 0 issues, 5 degrees, 4 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout: 7 rows, 0 issues, 4 degrees, 3 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: 6 rows, 0 issues, 20 denominator rows, 2 ratio-cap rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: 6 rows, 0 issues, 22 coefficient rows, 20 propagation rows, 2 intervalization rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: 6 rows, 0 issues, 320 root brackets, 30 zero floating weights, 2 intervalization rows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout: 6 rows, 0 issues, 320 midpoint weights, 30 repaired floating underflows, 320 direct interval obstructions, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate: 6 rows, 0 issues, 320 interval weights, 0 Taylor denominator obstructions, 30 repaired floating underflows, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate: 6 rows, 0 issues, 320 refined nodes, 320 interval weights, 2 below-one ratios, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate: 6 rows, 0 issues, 2 composed ratios, 3 tail sources, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda raw-moment obstruction matrix: 7 matrix rows, 0 issues, 3 exact counterexamples, 0 ready-to-apply rows
validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: 9 rows, 0 issues, 2 live routes, 2 rejected shortcuts, 0 ready-to-apply rows
validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues
validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues
validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found
validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues
validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction
validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows
validated result-language boundaries: scanned 191 markdown files, 0 overclaims
```
