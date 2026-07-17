# Signed-Hankel/Jensen Dependency Graph

Date: 2026-07-17

Status: dependency hygiene gate. This is not a proof of PF-infinity,
Laguerre-Polya membership, RH, or `Lambda <= 0`; it records how finite
evidence, countermodel gates, and open theorem targets are allowed to depend on
one another.

## Purpose

The signed-Hankel and Jensen-window route now has several exact algebra gates,
finite Arb diagnostics, and countermodel guards. The dependency graph makes
the proof boundary executable: finite evidence may support open theorem
targets, and countermodels may block invalid promotions, but finite evidence
has no direct proving edge to `lambda_le_0_goal`.

Machine-readable graph:

```text
work/rh_compute/results/signed_hankel_jensen_dependency_graph.json
```

Checker:

```text
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
```

Current result:

```text
validated signed-Hankel/Jensen dependency graph with 0 issues
```

## Core Invariant

The graph has an explicit conclusion node:

```text
lambda_le_0_goal:
  status `not_proved`
```

The only edges into that node are `would_imply_if_proved` edges from ledger
nodes whose status is still `open_target`, plus documentation edges from this
hygiene gate. The graph also has `blocked_by` edges back from
`lambda_le_0_goal` to the missing theorem targets, so the direction of the
remaining burden stays visible.

## Active Bridge Region

```text
signed_hankel_finite_certificate
hankel_sign_consistency_reduction_finite_certificate
shifted_hankel_sign_consistency_finite_certificate
arb_jensen_window_pf_obligation_diagnostic
arb_jensen_window_sturm_hyperbolicity_diagnostic
arb_jensen_window_sturm_d5_hyperbolicity_diagnostic
arb_jensen_window_sturm_d6_d12_hyperbolicity_diagnostic
jensen_window_sturm_pf_finite_consequence
```

These nodes support:

```text
target_signed_hankel_jensen_bridge
target_jensen_window_pf_bridge
```

Both targets remain open. The graph checker compares their statuses with
`work/rh_compute/results/proof_claim_ledger.json`.

## Obligation-Decomposition Region

```text
jensen_window_pf_bridge_obligation_ledger
jensen_window_pf_theorem_machinery_fit_matrix
jensen_window_pf_sign_regular_transfer_gap_matrix
jensen_window_pf_factorial_multiplier_split_audit
jensen_window_pf_structural_ansatz_matrix
jensen_window_pf_cauchy_binet_low_degree_scout
jensen_window_pf_log_concavity_frontier_scout
jensen_window_pf_ratio_condition_scout
jensen_window_pf_contraction_log_concavity_scout
jensen_window_pf_schur_shape_contract
jensen_window_pf_column_recurrence_contract
jensen_window_pf_column_recurrence_finite_coverage
arb_jensen_window_column_recurrence_stress
jensen_window_pf_reciprocal_positivity_route_matrix
jensen_window_pf_reciprocal_fraction_scout
jensen_window_pf_reciprocal_signed_j_fraction_scout
jensen_window_pf_signed_j_fraction_theorem_target
jensen_window_pf_modified_signed_model_target
jensen_window_pf_state_space_sign_lift_obstruction_scout
jensen_window_pf_oscillatory_resolvent_fit_matrix
jensen_window_pf_positive_readout_theorem_target
jensen_window_pf_positive_spectral_moment_obstruction
jensen_window_pf_nonordinary_positive_transform_ansatz_matrix
jensen_window_pf_nonpower_functional_low_degree_scout
jensen_window_pf_nonpower_functional_cone_candidate_matrix
jensen_window_pf_cauchy_binet_cone_frontier_matrix
jensen_window_pf_monotone_contraction_frontier_scout
jensen_window_pf_monotone_contraction_column_extension_scout
jensen_window_pf_monotone_contraction_sparse_degree6_scout
jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout
jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout
jensen_window_pf_monotone_contraction_all_m_counterexample
jensen_window_pf_monotone_contraction_theorem_target
jensen_window_pf_heat_flow_monotone_closure_scout
jensen_window_pf_heat_flow_boundary_threshold_lemma
jensen_window_pf_heat_flow_ratio_cone_invariance_lemma
jensen_window_pf_heat_flow_cone_entry_asymptotic_target
jensen_window_pf_heat_flow_infinite_cone_invariance_certificate
jensen_window_pf_defect_complete_monotonicity_scout
jensen_window_pf_multiplier_complete_monotonicity_frontier_scout
jensen_window_pf_multiplier_hausdorff_uniqueness_bridge
jensen_window_pf_multiplier_leading_atom_bound_certificate
jensen_window_pf_multiplier_unit_atomic_obstruction_certificate
jensen_window_pf_heat_flow_jensen_hierarchy_lemma
jensen_window_pf_rank_two_boundary_family_lemma
jensen_window_pf_cubic_reciprocal_defect_invariance_lemma
jensen_window_pf_cubic_m100_tail_entry_certificate
jensen_window_pf_cubic_forward_uniform_tail_certificate
jensen_window_pf_quartic_boundary_flow_obstruction
jensen_window_pf_quartic_double_root_threshold_lemma
jensen_window_pf_quartic_quintic_polar_contact_lemma
jensen_window_pf_cofinal_degree_polar_closure_lemma
jensen_window_pf_cofinal_scaling_limit_equivalence_gate
jensen_window_pf_multiplier_counting_measure_target
jensen_window_pf_mellin_multiplier_power_sum_obstruction
jensen_window_pf_phi_taylor_cone_entry_sign_scout
jensen_window_pf_negative_lambda_cone_entry_prefix_scout
jensen_window_pf_negative_lambda_finite_collar_contract
jensen_window_pf_negative_lambda_tail_barrier_scout
jensen_window_pf_negative_lambda_scaled_defect_frontier_scout
jensen_window_pf_negative_lambda_defect_recurrence_scout
jensen_window_pf_negative_lambda_log_curvature_bridge
jensen_window_pf_negative_lambda_bounded_log_curvature_target
jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction
jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target
jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge
jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge
jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout
jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations
jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget
jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress
jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation
jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan
jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout
jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate
jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress
jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget
jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout
jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target
jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout
jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout
jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target
jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout
jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate
jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate
jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix
jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix
jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate
jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate
jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate
jensen_window_pf_negative_lambda_gaussian_curvature_matrix
jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix
jensen_window_pf_negative_lambda_defect_tail_theorem_target
jensen_window_pf_negative_lambda_half_width_tail_target
jensen_window_pf_monotone_contraction_stress
```

This node sharpens:

```text
target_jensen_window_pf_bridge
```

It records exact, finite, open, conditional, rejected, and route-separated
obligations in `outputs/jensen_window_pf_bridge_obligations.md`. It is a
hygiene node only; it has no edge to `lambda_le_0_goal`.

The theorem-machinery fit matrix sharpens that obligation ledger by auditing
source-anchored total-positivity, PF, zero-preserver, sign-regularity, and
Laguerre-Polya theorem families. It records `0` ready-to-apply rows and also
has no edge to `lambda_le_0_goal`.

The sign-regular transfer gap matrix sharpens the same theorem-machinery audit
with exact degree-2 contact, degree-3/4 countermodel gates, and explicit
all-order/binomial/shift acceptance tests:

```text
outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md
work/rh_compute/results/jensen_window_pf_sign_regular_transfer_gap_matrix.json
```

The factorial multiplier split audit sharpens the route boundary by separating
the conditional Polya-Schur preservation step for `gamma_k=k!/(2*k)!` from the
false raw moment-window input theorem:

```text
outputs/jensen_window_pf_factorial_multiplier_split_audit.md
work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json
```

The reciprocal-gamma mixture gate now resolves the fixed-scale half of that
split and identifies the actual mixture obstruction:

```text
outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md
work/rh_compute/results/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.json
jensen_window_pf_reciprocal_gamma_mixture_sign_gate
```

Karlin's reciprocal-gamma theorem gives the required strict signature in every
order at one common positive scale. Determinant multilinearity produces an
exact independent-row scale integral for the Xi mixture, but its order-two
symmetrization changes sign and a positive three-atom measure reverses the
required minor. For the actual Xi measure, the completed ratio cone is exactly
the tilted concentration law `CV_n^2<=2/(2n+1)`. The graph therefore hands
higher compound concentration or a new positive coupling theorem to the
signed-Hankel and determinant-integral targets. This node has no edge to
`lambda_le_0_goal`.

The first higher compound layer now has an exact reciprocal-defect coordinate:

```text
outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md
work/rh_compute/results/jensen_window_pf_reciprocal_defect_compound_order3_gate.json
jensen_window_pf_reciprocal_defect_compound_order3_gate
```

The contiguous `3` by `3` signed-Hankel sign is equivalent to

```text
C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0.
```

This is unit-buffered log-concavity of reciprocal defects. An exact rational
countermodel lies strictly inside the full ratio cone, decreasing-defect and
increasing-scaled-defect cones, and both neighboring strict cubic Jensen
cones, yet has `C_n=-181/100` and the wrong positive Hankel sign. The graph
therefore blocks promotion from the unanchored scalar/cubic cones.

For the actual Xi sequence, the lambda=-100 entry layer is now closed:

```text
outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.json
jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate
```

A repaired Arb prefix through `n=317`, the certified anchor
`s_319>251/500`, and exact all-index scaled-defect growth prove

```text
C_n(-100)>57613471/66107054971>0, n>=318.
```

Thus every shifted contiguous order-three minor has the required sign at
`lambda=-100`. Forward propagation is now closed by a second node:

```text
outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order3_forward_invariance_certificate.json
jensen_window_pf_compound_order3_forward_invariance_certificate
```

Its exact cooperative equation and weighted infinite maximum principle prove

```text
C_n(lambda)>0 for every n>=0 and finite lambda>=-100,
D_(3,n)(0)<0 for every n>=0.
```

Arbitrary columns are closed by planar secant transfer:

```text
outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md
work/rh_compute/results/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.json
jensen_window_pf_order3_noncontiguous_secant_transfer_lemma
```

Contiguous negativity makes the planar edge slopes strictly decrease, and
positive weighted secant averaging transfers that orientation to every
increasing column triple. Thus reshaped-Hankel orders two and three are
complete at `lambda=0`. The graph now records compound order four and higher
as the live boundary. None of the order-three nodes has an edge to
`lambda_le_0_goal`.

The first order-four node makes that boundary explicit:

```text
outputs/jensen_window_pf_compound_order4_condensation_gate.md
work/rh_compute/results/jensen_window_pf_compound_order4_condensation_gate.json
jensen_window_pf_compound_order4_condensation_gate
```

Desnanot-Jacobi condensation reduces the contiguous sign to

```text
G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2).
```

The node contains 317 positive repaired lambda=-100 prefix margins and the
stronger finite cap `P_n<2/(n+3)^2`; its exact sufficient tail target is
`P_n<=4/(n+3)^2` for `n>=317`. That tail is now discharged downstream. The
strict rational lower-order countermodel with the wrong negative `H_(4,0)`
remains an active guard: entry uses the new curvature theorem, not promotion
from orders two and three.

The first-summand curvature bridge now sharpens that tail further:

```text
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_compound_order4_first_summand_curvature_bridge.json
outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
work/rh_compute/results/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json
outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.json
outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md
outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.md
outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md
outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md
outputs/arb_xi_lambda0_order4_prefix_certificate.md
outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md
jensen_window_pf_compound_order4_first_summand_curvature_bridge
jensen_window_pf_compound_order4_localized_curvature_compact_certificate
jensen_window_pf_compound_order4_gaussian_cumulant_ray_target
jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate
jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate
jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate
jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate
jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate
jensen_window_pf_compound_order4_exact_cumulant_remainder_budget
jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate
jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate
jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate
jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract
jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate
jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate
jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate
jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate
jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem
jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate
jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate
target_compound_order4_first_summand_curvature_ceiling
jensen_window_pf_compound_order4_m100_entry_certificate
jensen_window_pf_compound_order4_forward_flow_reduction
arb_xi_lambda0_order4_prefix_certificate
jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate
target_compound_order4_forward_invariance
```

It proves the stable factorization `g_1=d_1^2*(1-exp(-J_1))`, the real
floor `J_1(t)>=1/(7t)` for `t>=319`, and the full-kernel perturbation budget
`|P_n-P_n^(1)|<=2/(5k^2)`. A centered-difference Taylor bound localizes the
remaining curvature to same-point derivatives and three explicit envelopes.
The localized certificate node uses `107,452` outward-rounded Arb tiles and
`1,073` positive central blocks to prove

```text
(log g_1(t))''<=7/(2t^2), 319<=t<=V'(2).
```

The exact-corridor finite theorem preserves the common mode geometry on
`20,700` rational blocks and proves the localized ceiling on `2<=u<=20`.
The analytic ray theorem uses midpoint Hurwitz-zeta bounds, normalized
`H`-jet boxes, and a convergent logarithmic-defect majorant to prove

```text
t^2*U(t)<3011223637/866377000<7/2, u>=20.
```

The Gaussian-cumulant ray target rewrites that ray as explicit two-sided
standardized cumulant corridors through order eight. Exact epsilon-six algebra
proves the alternating factorial formal signature, and seven conditional
full-collar boxes clear the localized ceiling. A 1,800,000-block Arb theorem
proves every formal corridor on `2<=u<=20`; coefficient-positive leading and
jet-remainder gates prove the formal corridors for `u>=20`. Thus the formal
model is complete on `u>=2`. The exact epsilon-eight extension is audited and
its first omitted coefficient layer is also proved globally: a 1,800-block
centered Taylor theorem covers `2<=u<=20`, while coefficient-positive leading
and order-nine/ten jet gates cover `u>=20`. Composition reduces the remaining
exact-density theorem to explicit epsilon-ten residual budgets. The second-next
coefficient layer is now proved globally. A cancellation-preserving complex-
disk contract, two formal-tail theorems, two exact-tail theorems, a 5,400-block
partition extension, and a Bell-15 central theorem discharge those budgets.
Consequently all seven exact cumulant corridors hold for every `u>=2`, their
continuum corridor-to-`U` implication is complete, and the former node
`target_compound_order4_first_summand_curvature_ceiling` is now a discharged
exact theorem:

```text
K_1(t)<=7/(2t^2), every real t>=319.
```

The exact tent transfer and full-kernel perturbation then give
`P_n<=4/(n+3)^2` for every `n>=317`. The strict scaled-defect buffer converts
that estimate into the order-four sign, while the repaired prefix supplies
`0<=n<=316`. Thus the graph now contains the entry theorem

```text
H_(4,n)(-100)>0 for every n>=0.
```

The node `jensen_window_pf_compound_order4_forward_flow_reduction` proves the
exact one-sided system

```text
F_n'=alpha_n*F_(n+1)+beta_n*F_n, alpha_n>0.
```

Its original weighted-infimum reduction to a spatial coefficient ceiling
remains valid. The following compact-heat route now avoids needing that global
ceiling:

```text
outputs/jensen_window_pf_lambda0_first_summand_dominance_transfer.md
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md
outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md
jensen_window_pf_lambda0_first_summand_dominance_transfer
jensen_window_pf_uniform_superpolynomial_first_summand_dominance
jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem
```

Writing `T=-lambda`, covariance monotonicity transfers the first-summand
dominance estimate from `T=100` to every `0<=T<=100`. The same explicit
negative exponent proves superpolynomial suppression of all higher theta
summands, uniformly on that compact heat interval; fixed local logarithmic
differences inherit the suppression.

For the first theta summand, O'Sullivan's published all-order saddle expansion
applies uniformly to the suitable multiplier
`f_T(t)=exp(-T*(log t)^2/16)`. Exact differentiation of
`w=W(2k/pi)` and the compact-uniform remainder give

```text
Delta_k^m log R_T^(1)(k)=O(log(k)/k^m),  m=2,...,7,
```

uniformly for `0<=T<=100`. Exact seven-node Newton interpolation transfers
these differences to the normalized ratio coefficients. In the resulting
order-four determinant, all terms through `h^5` vanish and the universal first
term is

```text
768*G_2(lambda,M)^6*h(lambda,M)^6,
```

with `G_2->1` uniformly. Consequently there is one, possibly ineffective,
integer `N` such that

```text
H_(4,n)(lambda)>0 for n>=N and -100<=lambda<=0.
```

This uniform eventual tail confines the cooperative flow to finitely many
indices. Variation of constants and backward induction from index `N` combine
it with the all-shift entry theorem at `lambda=-100` to prove

```text
H_(4,n)(lambda)>0 for every n>=0 and every -100<=lambda<=0.
```

Thus `target_compound_order4_forward_invariance` is a discharged exact node;
in particular the previously non-effective lambda-zero splice is no longer
needed. This conclusion is only for contiguous order-four Hankel minors. It
does not establish arbitrary-column order four, any compound order at least
five, PF-infinity, the all-degree Jensen bridge, RH, or `Lambda<=0`.

The downstream finite-block transfer now closes the arbitrary-column part:

```text
outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md
work/rh_compute/results/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json
jensen_window_pf_order4_noncontiguous_total_positivity_transfer
```

For `B_(i,q)^(n,N)=A_(n+i+N-q)(0)`, every solid minor is the corresponding
signed contiguous Hankel minor. The completed layers through order four make
all initial minors positive, and the rectangular Gasca-Pena criterion makes
every reversed finite block strictly totally positive. Reversing an arbitrary
selected column set proves

```text
R_(4,n)(j_1,j_2,j_3,j_4)>0
for every n>=0 and 0<=j_1<j_2<j_3<j_4.
```

The same argument is an exact fixed-order theorem: signed contiguous layers
through `m` imply every arbitrary-column layer through `m`. It supplied the
arbitrary-column transfer used after the downstream order-five closure.
PF-infinity, the all-degree Jensen bridge, RH, and `Lambda<=0` remain open.

The order-five chain begins with the exact uniform-tail and flow reduction:

```text
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
work/rh_compute/results/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json
jensen_window_pf_compound_order5_uniform_tail_flow_reduction
```

The all-order Xi ratio expansion, heat-tilt extension through order eleven,
and higher-theta suppression give the exact determinant signature

```text
[h^0,...,h^10] det K=[0,0,0,0,0,0,0,0,0,0,294912*G_2^10]
```

and hence one eventual positive `H_5` tail uniformly on `-100<=lambda<=0`.
The completed `H_4>0` layer makes the exact order-five flow cooperative.
At the time of that reduction, finite confinement left the endpoint theorem

```text
H_(5,n)(-100)>0 for every n>=0,
equivalently H_(4,n+1)(-100)^2>H_(4,n)(-100)*H_(4,n+2)(-100).
```

That endpoint theorem was recorded under the stable identifier
`target_compound_order5_m100_entry`. It is now discharged by the chain below;
the identifier remains in the graph as an exact theorem node, not an open
target.

The stable-coordinate prefix certificate sharpens that endpoint target:

```text
outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order5_m100_prefix_certificate.json
jensen_window_pf_compound_order5_m100_prefix_certificate
```

Exact algebra factors `H_(5,n)=W_n*J_n` with `W_n>0` in the completed lower
cone. Seven source files provide 325 outward-rounded Arb coefficient balls,
and direct evaluation proves

```text
J_n(-100)>0 and H_(5,n)(-100)>0 for every 0<=n<=316.
```

The weakest relative lower bound is above `0.006269` at `n=316`. This supplies
the finite half of the endpoint composition.

That tail has a scalar curvature reduction:

```text
outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md
work/rh_compute/results/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.json
jensen_window_pf_compound_order5_m100_tail_curvature_reduction
```

With `k=n+4`, exact logarithmic algebra gives

```text
J_n>0 iff C_n<-4log(x_k),
C_n=Delta^2 log(F_n)-Delta^2 log(d_(n+3)).
```

The inherited defect anchor and coefficient-positive rational arithmetic show
that the single ceiling

```text
C_n<=100/(n+4)^2 for every n>=317
```

is sufficient for the complete tail. The constant has a factor-above-27
reserve over the measured scaled curvature at the splice.

The exact first/full bridge decomposes that budget:

```text
outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_compound_order5_first_summand_curvature_bridge.json
jensen_window_pf_compound_order5_first_summand_curvature_bridge
```

Positive floors at both stable logarithmic layers prevent division by a raw
near-zero determinant. Exact coefficient-positive perturbation arithmetic
proves

```text
|C_n-C_n^(1)|<=37/(n+4)^2, n>=317,
```

and reduces the remaining `63` part to

```text
q_1''(t)<=60/t^2 for every real t>=320.
```

Three rigorous interval nodes close that continuous theorem:

```text
outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.json
jensen_window_pf_compound_order5_nested_curvature_compact_certificate

outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.json
jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate

outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.json
jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate
```

The compact certificate reuses the hashed `107452`-tile H-derivative cache
and proves the ceiling on `320<=t<=V'(2)`. One hundred collar-extension tiles
and 1850 exact-cumulant blocks cover `2<=u<=20`. The normalized asymptotic
box proves `t^2q_1''(t)<9.159` for every `u>=20`. Their union covers every
real `t>=320`.

The endpoint composition is:

```text
outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order5_m100_entry_certificate.json
jensen_window_pf_compound_order5_m100_entry_certificate
target_compound_order5_m100_entry
```

The tent bound contributes less than `63/(n+4)^2`; the full-kernel transfer
contributes at most `37/(n+4)^2`. The resulting scalar ceiling proves the
complete analytic tail, and the prefix supplies the remaining shifts. Hence

```text
H_(5,n)(-100)>0 for every integer n>=0.
```

Finally,

```text
outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.json
jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate
```

composes endpoint entry, the uniform eventual tail, cooperative variation of
constants, and the fixed-order initial-minor transfer. It proves

```text
H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0,
R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0
```

for every increasing column quintuple on the same interval. Thus fixed order
five is complete.

The downstream order-six chain begins with

```text
outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md
work/rh_compute/results/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json
jensen_window_pf_compound_order6_uniform_tail_flow_reduction
```

The heat-tilt expansion through order sixteen, higher-theta suppression, and
exact 720-permutation determinant audit give the signed first term

```text
1132462080*G_2^15*h^15
```

and one uniform eventual positive `Q_6=-H_6` tail. The completed `H_5>0`
layer makes the exact signed flow cooperative.

The finite endpoint half is

```text
outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order6_m100_prefix_certificate.json
jensen_window_pf_compound_order6_m100_prefix_certificate
```

Stable signed condensation and 327 outward-rounded coefficient balls prove

```text
Q_(6,n)(-100)>0 for every 0<=n<=316.
```

The canonical tail reduction and first/full bridge are

```text
outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md
jensen_window_pf_compound_order6_m100_tail_curvature_reduction

outputs/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.md
jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension

outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md
jensen_window_pf_compound_order6_first_summand_curvature_bridge
```

Exact factorization gives `H_(5,n)=A_(n+4)^5exp(p(n+4))` and reduces entry to
`P_k<=320/k^2` for `k=n+5>=322`. The inverse-seventh-power defect theorem and
a degree-75 coefficient-positive transfer prove

```text
|P_k-P_k^(1)|<100/k^2
```

and reduce the remaining `201` part to `p_1''(t)<=200/t^2` on `t>=321`.

The required derivative corridors and three continuous ranges are graph
nodes:

```text
outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md
jensen_window_pf_compound_order6_high_cumulant_coarse_corridor

outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md
jensen_window_pf_compound_order6_nested_curvature_compact_certificate

outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md
jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate

outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md
jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate
```

The exact high-cumulant node supplies orders nine and ten. Thirty-eight
adaptive compact blocks cover `321<=t<=V'(2)`; a mode-two collar and 17,999
rational blocks cover `2<=u<=20`; the normalized asymptotic box proves
`t^2p_1''(t)<22.769` for `u>=20`. Their union proves the continuous theorem.

The endpoint composition

```text
outputs/jensen_window_pf_compound_order6_m100_entry_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order6_m100_entry_certificate.json
jensen_window_pf_compound_order6_m100_entry_certificate
target_compound_order6_m100_entry
```

combines the strict `301<320` analytic budget with the prefix and proves

```text
Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every n>=0.
```

Finally,

```text
outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.json
jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate
```

composes endpoint entry, the uniform tail, cooperative variation of constants,
and the all-fixed-order initial-minor transfer. It proves signed contiguous and
arbitrary-column order six for every shift on `-100<=lambda<=0`. The graph has
only `sharpens` edges from this theorem into the two all-order bridge targets.

The new all-fixed-order tail node is

```text
outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md
work/rh_compute/results/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json
jensen_window_pf_graded_kernel_vandermonde_all_order_lemma
```

Mixed-kernel valuation and Cauchy-Binet force determinant degree
`D=binom(m,2)` at every fixed order. The unique equality block is the
Vandermonde block and gives

```text
[h^D]Q_m=2^D*prod_(j=1)^(m-1)j!*G_2^D>0.
```

Arbitrary fixed-order suitable-multiplier truncation, direct finite-difference
control of its remainder, and higher-theta suppression prove

```text
for every fixed m, exists N_m such that
Q_(m,n)(lambda)>0 for n>=N_m and -100<=lambda<=0.
```

The threshold may depend on `m` and is non-effective. The graph therefore
records only `supports` edges from the three analytic inputs and does not
interpret this as a common-tail or all-shift theorem.

Its order-seven specialization and exact flow reduction are

```text
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
work/rh_compute/results/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.json
jensen_window_pf_compound_order7_uniform_tail_flow_reduction
```

Here `D=21`, `epsilon_7=-1`, and the positive signed leading coefficient is
`52183852646400*G_2^21`. Exact condensation gives

```text
Q_(7,n)H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)Q_(6,n+2),
```

while the completed order-six cone makes the order-seven flow cooperative.
Thus all-shift endpoint entry at `lambda=-100` would imply contiguous and
arbitrary-column order seven on the complete heat interval. The exact
lower-cone countermodel in the node prevents promotion from orders at most six.

The rigorous finite half of endpoint entry is

```text
outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order7_m100_prefix_certificate.json
jensen_window_pf_compound_order7_m100_prefix_certificate
```

Stable Arb evaluation, including a twelve-coefficient retained-integral
repair, proves

```text
Q_(7,n)(-100)>0 for every 0<=n<=314,
L_314>9/1000.
```

The remaining endpoint tail has the canonical reduction

```text
outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md
work/rh_compute/results/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.json
jensen_window_pf_compound_order7_m100_tail_curvature_reduction
```

One additional stable gap gives `Q_(6,n)=A_(n+5)^6exp(r(n+5))`. With
`k=n+6`, order-seven positivity is equivalent to

```text
6log(x_k)+R_k<0,
R_k=r(k-1)-2r(k)+r(k+1).
```

The completed defect anchor and an exact coefficient-positive comparison
reduce every missing `n>=315` sign to the single estimate later discharged
by the continuous certificates:

```text
R_k<=900/k^2 for every k>=321.
```

The first-summand input is now strengthened by

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.json
jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension
```

The rebalanced split `a(k)=log(k)/10` uses the high-region exponential
reserve to prove

```text
0<=delta_k<2/k^8 for every k>=300.
```

Fourteen endpoint and half-line derivative gates are strict. This is the
power needed to cross one additional stable logarithm without exhausting the
inverse-square endpoint budget.

The resulting exact transfer node is

```text
outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_compound_order7_first_summand_curvature_bridge.json
jensen_window_pf_compound_order7_first_summand_curvature_bridge
```

The completed order-six theorems and two finite relative-H5 margins give

```text
min(T_j,T_j^(1))>=3/(2j), j>=320.
```

The full four-layer perturbation chain has a degree-102 shifted numerator
with all 103 coefficients positive and proves

```text
|R_k-R_k^(1)|<262/k^2, k>=321.
```

Tent integration therefore reduces the endpoint tail to exactly

```text
r_1''(t)<=600/t^2 for every real t>=320
 =>R_k<601/k^2+262/k^2=863/k^2<900/k^2.
```

The compact part of this continuous target is now discharged by

```text
outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.json
jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate
```

Exact-potential shifted jets and dimensionless `H2..H21` collars prove

```text
r_1''(t)<=600/t^2 for every real 320<=t<=1000.
```

All 186 rational blocks pass. The worst scaled upper is below `50.911`, and
the weakest dimensionless `T` floor exceeds `2.069`.

The next compact interval is discharged by

```text
outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.json
jensen_window_pf_compound_order7_nested_curvature_compact_certificate
```

Aligned `H2..H12` quadrature caches and strict `t+-5` common collars prove

```text
r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2).
```

All 82 adaptive rational mode blocks pass. The worst scaled upper is below
`358.733`, and the weakest `T` floor exceeds `0.00010247`.

The two new cumulant orders needed beyond mode two are certified by

```text
outputs/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.md
work/rh_compute/results/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.json
jensen_window_pf_compound_order7_high_cumulant_coarse_corridor
```

The exact epsilon-ten recurrence gives normalized `kappa_11,kappa_12` caps
below `14001` on `2<=u<=20` and below `700001` on `u>=20`. These feed the two
closed ray nodes:

```text
outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.md
jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate
jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate
```

The mode-two collar and 17,999 rational blocks prove the target on
`2<=u<=20`; one normalized stable-log interval proves
`t^2r_1''(t)<55.541<100` for every `u>=20`. Consequently the four closed
ranges prove

```text
r_1''(t)<=600/t^2 for every real t>=320.
```

The exact endpoint composition is

```text
outputs/jensen_window_pf_compound_order7_m100_entry_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order7_m100_entry_certificate.json
jensen_window_pf_compound_order7_m100_entry_certificate
```

Tent integration and the complete-kernel transfer give

```text
R_k<601/k^2+262/k^2=863/k^2<900/k^2, k>=321,
Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every n>=0.
```

The formerly conditional flow theorem then composes to

```text
outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md
work/rh_compute/results/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.json
jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate
```

It proves signed contiguous and arbitrary-column order seven for every shift
and every `-100<=lambda<=0`. Only `sharpens` edges leave this fixed-order
theorem for the two all-order bridge targets. The next layer cannot be
promoted from this lower cone alone.

The completed fixed-order-eight chain begins with

```text
jensen_window_pf_compound_order8_uniform_tail_flow_reduction
jensen_window_pf_compound_order8_m100_prefix_certificate
jensen_window_pf_compound_order8_m100_tail_curvature_reduction
jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension
jensen_window_pf_compound_order8_first_summand_curvature_bridge
```

At `m=8`, `D=28`, `epsilon_8=1`, and the universal positive signed term is

```text
33664847019245568000*G_2^28*h^28.
```

The all-fixed-order theorem supplies a compact-uniform eventual `Q_8` tail.
Signed condensation and the affine-Hankel/Plucker identity give

```text
Q_(8,n)Q_(6,n+2)=Q_(7,n+1)^2-Q_(7,n)Q_(7,n+2),
Q_n'=a_n Q_(n+1)+b_n Q_n,
a_n=(4n+58)Q_(7,n)/Q_(7,n+1)>0.
```

Thus all-shift endpoint entry at `lambda=-100` is the sole new zeta-specific
input for fixed order eight. It cannot be inferred from the completed lower
cone: the
exact sequence `(1,1,1/2!,...,1/13!,1/87120000000)` has all available signed
contiguous minors through order seven positive but

```text
Q_(8,0)=-463/69210459277322870541707704182767616000000000000000<0.
```

The finite endpoint component uses 1,257 outward-rounded coefficient balls
at 2048 bits. The stable relative coordinate

```text
M_n=Q_(7,n+1)^2/(Q_(7,n)Q_(7,n+2))-1
```

is rigorously positive for every `0<=n<=1242`. Signed condensation proves
`Q_(8,n)(-100)>0` throughout the same 1,243-shift prefix, and the weakest
margin, at `n=1242`, exceeds `1/300`.

For the tail, exact factorization reduces `n>=1243`, or `k>=1250`, to

```text
W_k<=4300/k^2.
```

The inverse-ninth-power dominance theorem and fifth-gap bridge prove

```text
|W_k-W_k^(1)|<190/k^2,
s_1''(t)<=4000/t^2 on t>=999
 => W_k<4001/k^2+190/k^2=4191/k^2<4300/k^2.
```

The continuous premise is supplied by the following exact and interval
nodes:

```text
jensen_window_pf_compound_order8_high_cumulant_coarse_corridor
jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate
jensen_window_pf_compound_order8_nested_curvature_compact_certificate
jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate
jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate
```

The shifted cover proves the stronger `2000/t^2` ceiling on `699<=t<=999`.
The compact cover handles `999<=t<=V'(2)`, 17,999 rational blocks handle
`2<=u<=20`, and normalized `H2-H14` boxes handle every `u>=20`, where the
scaled curvature upper is below `134.49`. The exact high-cumulant corridor
supplies the orders thirteen and fourteen caps on both saddle rays.

These components compose in

```text
jensen_window_pf_compound_order8_m100_entry_certificate
jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate
```

to prove

```text
Q_(8,n)(-100)=H_(8,n)(-100)>0 for every n>=0,
Q_(8,n)(lambda)=H_(8,n)(lambda)>0
  for every n>=0 and -100<=lambda<=0.
```

The fixed-order initial-minor theorem then transfers the conclusion to every
arbitrary-column signed reshaped-Hankel minor through order eight. The graph
uses only `supports` edges inside this fixed-order composition and `sharpens`
edges from its final theorem to the all-order targets.

The completed fixed-order-nine chain is now recorded by the nodes

```text
jensen_window_pf_compound_order9_uniform_tail_flow_reduction
jensen_window_pf_compound_order9_m100_prefix_certificate
jensen_window_pf_compound_order9_m100_tail_curvature_reduction
jensen_window_pf_compound_order9_first_summand_curvature_bridge
jensen_window_pf_compound_order9_high_cumulant_coarse_corridor
jensen_window_pf_compound_order9_shifted_point_h0_h8_cache
jensen_window_pf_compound_order9_localized_lower_bridge_certificate
jensen_window_pf_compound_order9_nested_curvature_compact_certificate
jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate
jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate
jensen_window_pf_compound_order9_first_summand_curvature_certificate
jensen_window_pf_compound_order9_m100_finite_splice_certificate
jensen_window_pf_compound_order9_m100_entry_certificate
jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate
```

The eventual tail and endpoint coordinate are

```text
det K=347485857744891213250560000*G_2^36*h^36+o(h^36),
Q_(9,n)Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)Q_(8,n+2).
```

The finite prefix proves `Q_(9,n)(-100)>0` through `n=1240`. The tail is
reduced to `Y_k<=4900/k^2`, while the six-layer complete-to-first transfer
contributes less than `550/k^2`. The hash-bound localized cover, compact
cover, finite saddle ray, and asymptotic ray compose to

```text
w_1''(t)<=4200/t^2 for every real t>=1250.
```

Tent integration contributes less than `4201/k^2`; hence the full endpoint
tail is below `4751/k^2<4900/k^2`. Two retained-integral splice signs at
`n=1241,1242` join that tail to the prefix. The endpoint and cooperative-flow
nodes therefore prove

```text
Q_(9,n)(lambda)=H_(9,n)(lambda)>0
for every n>=0 and -100<=lambda<=0,
```

and transfer the result to arbitrary columns through fixed order nine.

The order-ten graph now records a delayed recovery rather than an all-shift
endpoint entry.  Three hash-bound compact caches, the localized cover, compact
cover, finite saddle ray, and asymptotic ray compose to

```text
z_1''(t)<=4200/t^2 for every real t>=1251.
```

The exact seven-layer transfer, finite endpoint block, and two-index splice
then prove

```text
Q_(10,n)(-100)<0 for n=0,1,2,3,
Q_(10,n)(-100)>0 for every integer n>=4.
```

The delayed cooperative heat lemma propagates only the positive `n>=4` ray;
an independent Arb certificate supplies the four missing signs at lambda
zero.  Consequently

```text
Q_(10,n)(0)>0 for every integer n>=0,
epsilon_k R_(k,n)(j_1,...,j_k)(0)>0
  for 1<=k<=10 and every admissible shift and column set.
```

Only `supports` edges enter this fixed-order composition, and only `sharpens`
edges leave its final node for the open bridge targets.  The graph therefore
records a genuine order-ten lambda-zero theorem without promoting it to order
eleven, PF-infinity, RH, or `Lambda<=0`.

The next fixed-order node is deliberately subordinate and remains open:

```text
jensen_window_pf_compound_order11_curvature_bridge_target
```

Its finite endpoint prefix through `n=1242`, four-row lambda-zero complement,
power-thirteen full-kernel transfer, and delayed heat theorem are all recorded
as `supports` inputs.  Exactly one premise remains: a half-line proof of
`y_1''(t)<=6000/t^2` for `t>=1252`.  The node has only `sharpens` edges to the
two global bridge targets, so neither selected-point evidence nor the tested
localized interval core can be promoted to order eleven or `Lambda<=0`.

The new node

```text
jensen_window_pf_all_order_endpoint_heat_reduction
```

then records the arbitrary-order determinant identities

```text
Q_(m,n)'=(4n+8m-6)delta(Q_(m,n)),
Q_(m-1,n+1)delta(Q_(m,n))
 =Q_(m-1,n)Q_(m,n+1)+delta(Q_(m-1,n+1))Q_(m,n),
```

and the exact quantifier composition

```text
[Q_(m,n)(-100)>0 for every m>=10,n>=0]
iff
[Q_(m,n)(lambda)>0 for every m>=10,n>=0
 and -100<=lambda<=0].
```

The eventual-tail threshold may depend on `m`: each fixed order is completed
before induction advances. The graph uses `supports` edges from the fixed-order
tail, completed order-nine base, and initial-minor theorem into this reduction.
Only `sharpens` edges leave it for the two open bridge targets. Its conditional
equivalence remains exact, but the static `m>=10` endpoint antecedent is now
rejected by the order-ten node below. PF-infinity, RH, and `Lambda<=0` are not
proved or disproved by that route failure.

The endpoint hierarchy now has an exact normalized Schur coordinate:

```text
outputs/jensen_window_pf_endpoint_deep_schur_coordinate.md
work/rh_compute/results/jensen_window_pf_endpoint_deep_schur_coordinate.json
jensen_window_pf_endpoint_deep_schur_coordinate
jensen_window_pf_endpoint_order10_counterexample
jensen_window_pf_endpoint_pf3_boundary_counterexample
jensen_window_pf_endpoint_deep_schur_literature_fit
```

For `h_k=A_k(-100)/A_0(-100)`, column reversal and Jacobi-Trudi give

```text
Q_(m,n)(-100)=A_0(-100)^m s_((n+m-1)^m)(h).
```

Arbitrary increasing Hankel columns map bijectively to partitions
`lambda_1>=...>=lambda_m>=m-1`. The graph therefore uses `supports` edges
from the all-order endpoint reduction and fixed-order initial-minor theorem
into the deep-Schur node, followed only by `sharpens` edges to the open bridge
targets. The deep-Schur and Toda coordinates also support the direct
order-ten counterexample node:

```text
outputs/jensen_window_pf_endpoint_order10_counterexample.md
work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json
python work/rh_compute/scripts/check_jensen_window_pf_endpoint_order10_counterexample.py
```

Independent `4096`-bit Arb Hankel, Jacobi-Trudi, and Toda checks prove

```text
Q_(10,n)(-100)<0 for n=0,1,2,3,
s_((N^10))(h)<0 for N=9,10,11,12,
Q_(10,4)(-100)>0.
```

These are four negative required order-ten deep rectangles, so the proposed
all-order endpoint hierarchy is false. The graph sends a `blocks_promotion`
edge from this node to `target_signed_hankel_jensen_bridge` and a `sharpens`
edge to the direct Jensen target. It has no edge to `lambda_le_0_goal`.

The rigorous endpoint coefficient balls also prove

```text
s_(1,1,1)(h)<0.
```

Since `(1,1,1)` lies outside the required deep order-three cone, the graph
uses a `blocks_promotion` edge to prevent strengthening the endpoint target to
full PF-infinity. The primary-literature fit node records that initial-minor
and orientation theorems support the coordinate, while full Edrei positivity
is too strong and matrix-power eventual positivity is a different notion.
None of these coordinate or guard nodes has an edge to `lambda_le_0_goal`.

The rectangular coordinate now has an exact Toda refinement and two explicit
route guards:

```text
outputs/jensen_window_pf_deep_schur_toda_boundary_gate.md
jensen_window_pf_deep_schur_rectangular_toda_coordinate
jensen_window_pf_deep_schur_moving_tail_boundary_counterexample
jensen_window_pf_strict_schur_jensen_counterexample
```

Writing `tau_(m,N)=s_((N^m))(h)`, Desnanot-Jacobi gives

```text
tau_(m+1,N) tau_(m-1,N)
 = tau_(m,N)^2-tau_(m,N-1) tau_(m,N+1).
```

Once the lower row is positive, the next rectangle is positive exactly when
the current width row is strictly log-concave. This is an exact coordinate,
not an induction theorem: the Toda numerator is subtractive. The order-ten
gate certifies that it is negative at the first four admissible widths.

The first guard concerns the ordinary tail `b_k^(s)=h_(s+k)/h_s`, with
negative indices reset to zero. Its Jacobi-Trudi determinants agree with
translated endpoint shapes only for deep tail shapes. At
`r=2,s=1,mu=(0,0)`, the exact mismatch is `h_0 h_2/h_1^2>0`. Thus deep
positivity does not imply moving-tail PF by this zero-reset translation.

The second guard uses the strictly Schur-positive Edrei specialization

```text
H(z)=exp(z/100)/((1-z)(1-2z)).
```

Its shift-zero cubic Jensen polynomial has exact discriminant
`-222484532394597/2000000000000<0`, so it is not hyperbolic. Therefore even
strict full unweighted Schur/PF positivity is not, by itself, the missing
Jensen implication. The graph sends this counterexample only to the open
bridge targets. Together with the direct endpoint failure, any successful
bridge must use a weaker Xi/Phi-specific structure that does not assume the
all-shift signed-Hankel cone. None of these nodes has an edge to
`lambda_le_0_goal`.

The current primary-source theorem audit adds a non-promotion edge:

```text
outputs/jensen_window_pf_order_moment_transport_fit_gate.md
work/rh_compute/results/jensen_window_pf_order_moment_transport_fit_gate.json
jensen_window_pf_order_moment_transport_fit_gate
```

The positive Gamma-normalized order-moment transport theorem would require a
completely monotone transformed first-summand kernel. The exact origin
derivative is positive, so that hypothesis fails, and the theorem has ordinary
positive Hankel orientation rather than reciprocal sign regularity. The node
records why that shortcut could not discharge the historical curvature
target; the successful proof instead uses the exact-corridor route.

The structural ansatz matrix sharpens the theorem-machinery audit into
candidate, blocked, and rejected proof-search rows:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
outputs/jensen_window_pf_schur_shape_contract.md
work/rh_compute/results/jensen_window_pf_schur_shape_contract.json
outputs/jensen_window_pf_column_recurrence_contract.md
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json
outputs/arb_jensen_window_column_recurrence_stress.md
work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json
outputs/jensen_window_pf_reciprocal_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json
outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json
outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json
outputs/jensen_window_pf_signed_j_fraction_theorem_target.md
work/rh_compute/results/jensen_window_pf_signed_j_fraction_theorem_target.json
outputs/jensen_window_pf_modified_signed_model_target.md
work/rh_compute/results/jensen_window_pf_modified_signed_model_target.json
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json
outputs/jensen_window_pf_positive_readout_theorem_target.md
work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json
outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.json
outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md
work/rh_compute/results/jensen_window_pf_nonpower_functional_low_degree_scout.json
outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md
work/rh_compute/results/jensen_window_pf_nonpower_functional_cone_candidate_matrix.json
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json
outputs/jensen_window_pf_monotone_contraction_frontier_scout.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_frontier_scout.json
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_theorem_target.json
outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md
work/rh_compute/results/jensen_window_pf_heat_flow_monotone_closure_scout.json
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
work/rh_compute/results/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.json
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.json
outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md
work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_contract.json
outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_scout.json
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.json
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k80_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k80_scout.json
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_defect_recurrence_scout.json
outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_log_curvature_bridge.json
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_target.json
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md
work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.json
outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.json
outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json
outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.json
outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.json
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_defect_tail_theorem_target.json
outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_half_width_tail_target.json
outputs/jensen_window_pf_monotone_contraction_stress.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json
```

It checks positive Cauchy-Binet, planar-network/production-matrix,
determinant-integral, preserver, direct-Hankel, and finite-grid ansatz rows
against the exact low-degree algebra and finite countermodel kill gate. It
records `0` ready-to-apply rows and has no edge to `lambda_le_0_goal`. The
Schur shape contract sharpens these structural rows by recording the bounded
finite-support Jacobi-Trudi shapes that a positive Schur, network, or
determinant-integral route must cover. The column recurrence contract sharpens
the Schur shape contract by isolating the elementary-symmetric recurrence
conditions for the hard frontier column shapes. The finite coverage map
supports that recurrence target on checked zeta-window grids, while preserving
that finite evidence has no proving edge to `lambda_le_0_goal`. The Arb
column-recurrence stress grid supports the finite coverage map across a wider
degree and recurrence-size range, but remains a diagnostic node only.
The extended reciprocal coefficient stress grid supports the reciprocal
target directly with 72,600 normalized positive rows, but also remains a
diagnostic node only.

The reciprocal-positivity route matrix sharpens the column recurrence contract
by classifying the exact theorem mechanisms that could prove
`[t^m]1/H(-t)>=0`: real-root/PF endpoint language, positive renewal or
convolution kernels, positive continued fractions, production matrices,
Kaluza-type sign theorems, generic ratio conditions, finite recurrence stress,
and the remaining all-Schur lift. It records `0` ready-to-apply rows and has
no edge to `lambda_le_0_goal`.

The reciprocal fraction scout sharpens that matrix by rejecting the standard
positive Stieltjes/Jacobi S-fraction or J-fraction route at the first
nontrivial parameter:

```text
a_2 = -g_2/g_1
lambda_1 = -g_2
```

It validates 3 symbolic rows and 735 finite Arb zeta-window sign rows. It is a
diagnostic node only and has no edge to `lambda_le_0_goal`.

The reciprocal signed J-fraction scout supports the surviving signed-fraction
subroute by validating:

```text
3,675 signed reciprocal-Hankel determinant rows
2,940 ordinary Jacobi lambda rows with negative sign
2,940 signed kappa_n=-lambda_n rows with positive sign
```

It is finite evidence only, remains a diagnostic node, and has no edge to
`lambda_le_0_goal`.

The reciprocal signed Jacobi beta scout sharpens the same surviving signed
subroute by validating:

```text
3,675 signed Jacobi beta rows
2,940 positive beta rows
630 negative beta_1 rows
105 terminal degree-2 zero-containing rows
```

It is finite evidence only, remains a diagnostic node, and has no edge to
`lambda_le_0_goal`.

The reciprocal Motzkin path obstruction scout rejects the raw ordinary
J-fraction path model as manifest positivity by validating:

```text
735 mu_2 cancellation rows
630 beta_1 diagonal obstruction rows
```

It is finite evidence only, remains a diagnostic/rejection node for the raw
ordinary path model, leaves modified signed models open, and has no edge to
`lambda_le_0_goal`.

The state-space sign-lift obstruction scout is derived from those mu_2 rows:

```text
outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json
```

It rejects the absolute-value sign-state cover of raw Motzkin paths by
recording:

```text
absolute_lift_mu2 - mu_2 = 2*kappa_1 > 0
```

across 735 derived rows. It is finite derived evidence only, leaves genuinely
modified state-space doubled models open, and has no edge to
`lambda_le_0_goal`.

The reciprocal Motzkin parity-lift obstruction scout rejects global
length-parity signs and diagonal sign conjugation as repairs of that raw
ordinary path model by validating:

```text
5,145 same-length mixed-sign witness rows
```

It is finite evidence only, remains a diagnostic/rejection node for cheap
sign repairs, leaves state-space modified models open, and has no edge to
`lambda_le_0_goal`.

The oscillatory/resolvent fit matrix sharpens the signed J-fraction theorem
target, the modified signed-model target, and the reciprocal-positivity route
matrix:

```text
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json
```

It validates 8 fit/misfit rows and 0 ready-to-apply rows. It rejects ordinary
entrywise Jacobi powers, diagonal similarity, absolute-value majorants,
classical oscillatory spectral conclusions, indefinite moment language, and
finite signed patterns as standalone coefficient-positivity proofs. It leaves
only positive spectral-transform and Xi/Phi positive-kernel variants live
conditionally. It is a hygiene node only and has no edge to `lambda_le_0_goal`.

The positive-readout theorem target sharpens those two surviving variants:

```text
outputs/jensen_window_pf_positive_readout_theorem_target.md
work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json
```

It validates 8 candidate rows, 0 ready-to-apply rows, and 2 live foundational
routes. It requires an exact positive scalar readout for E(t)=1/H(-t), either
by a noncircular positive spectral transform or by an Xi/Phi-specific positive
resolvent kernel. It rejects abstract wrappers, endpoint factorization, finite
quadrature, raw signed readouts, and absolute-value majorants as standalone
proofs. It is a hygiene node only and has no edge to `lambda_le_0_goal`.

The positive spectral moment obstruction rejects the ordinary moment-measure
interpretation of that spectral-transform row:

```text
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
work/rh_compute/results/jensen_window_pf_positive_spectral_moment_obstruction.json
```

It records the exact identity `Delta_2=-g_2` and validates 735 finite
`Delta_2<0` rows from the signed J-fraction scout. It leaves nonordinary
positive transforms, positive kernels, and Xi/Phi-specific readouts open. It
is a diagnostic node only and has no edge to `lambda_le_0_goal`.

The nonordinary positive-transform ansatz matrix then records what remains
after that obstruction:

```text
outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.json
```

It validates 8 ansatz rows with 0 ready-to-apply rows and 3 live ansatz rows:
non-power positive functionals, Xi/Phi positive kernels, and genuinely
modified exact state-space transfer models. It is a diagnostic theorem-search
node only and has no edge to `lambda_le_0_goal`.

The nonpower functional low-degree scout sharpens the live non-power
functional row:

```text
outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md
work/rh_compute/results/jensen_window_pf_nonpower_functional_low_degree_scout.json
```

It validates exact reciprocal coefficient formulas through `mu_6` and signed
composition counts through `m=8`, recording that npt_04 needs a genuine
positive cone, basis, and functional to absorb the signed `mu_2` and `mu_3`
cancellation identities. It is a diagnostic node only and has no edge to
`lambda_le_0_goal`.

The nonpower functional cone candidate matrix classifies candidate positive
cones for that contract:

```text
outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md
work/rh_compute/results/jensen_window_pf_nonpower_functional_cone_candidate_matrix.json
```

It validates 8 cone rows with 0 ready-to-apply rows and 2 live cone rows,
rejecting raw g-coordinate, standalone ratio/log-concavity, tautological, and
endpoint PF/LP cones while leaving Xi/Phi kernel and
Cauchy-Binet/determinant-integral cone routes live. It is a diagnostic node
only and has no edge to `lambda_le_0_goal`.

The Cauchy-Binet cone frontier matrix sharpens the live determinant-integral
cone:

```text
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json
```

It validates 8 frontier rows with 0 ready-to-apply rows and 2 live frontier
rows. It rejects treating selected low-degree Bernstein certificates or
adjacent log-concavity as a Cauchy-Binet cone proof, and leaves only a hard
column-frontier determinant integral or an all-shape Cauchy-Binet/Andreief
kernel live. It is a diagnostic node only and has no edge to
`lambda_le_0_goal`.

The monotone-contraction frontier scout sharpens the hard column-frontier
subroute:

```text
outputs/jensen_window_pf_monotone_contraction_frontier_scout.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_frontier_scout.json
outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_column_extension_scout.json
outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree6_scout.json
outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.json
outputs/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.json
outputs/jensen_window_pf_monotone_contraction_all_m_counterexample.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_all_m_counterexample.json
```

It validates 2 exact hard-row certificates with 88 positive Bernstein
coefficients and 210 finite zeta diagnostic rows. It shows that the first hard
column-frontier polynomials are positive under the sharper sufficient condition
`x1 <= x2 <= x3`, while the rational log-concavity countermodel violates that
condition. It remains a diagnostic node only: it does not prove that all zeta
windows satisfy monotone contractions, does not construct a determinant
integral, and has no edge to `lambda_le_0_goal`.

The column-extension scout extends that exact certificate in a bounded
column-family range: degree 3 through `m=10`, degree 4 through `m=7`, and
degree 5 through `m=8`, with 25 column rows, 3,329 positive Bernstein
coefficients, and 3 rows beyond the original first hard frontier. It is still
bounded theorem-search algebra, not an all-`m` or all-shape theorem.

The sparse degree-6 scout extends the same bounded certificate family through
degree 6 and `m=10`, using sparse exact arithmetic to validate 63,347 strictly
positive Bernstein coefficients with 0 zero or negative entries. It remains a
finite diagnostic and has no edge to `lambda_le_0_goal`.

The sparse degree-7 frontier scout extends the same one-shot global Bernstein
certificate to degree 7 through `m=9`, validating 670,891 strictly positive
coefficients, then records a certificate frontier at `m=10` with 126 negative
Bernstein coefficients. This is a certificate obstruction only, not a
polynomial-negativity counterexample, and has no edge to `lambda_le_0_goal`.

The sparse degree-7 subdivision scout repairs that `m=10` certificate frontier
by splitting only `s0` into the three dyadic slabs `[0,1/2]`, `[1/2,3/4]`,
and `[3/4,1]`, validating 785,400 strictly positive slab Bernstein
coefficients. It remains finite theorem-search algebra and has no edge to
`lambda_le_0_goal`.

The all-m counterexample blocks a tempting overpromotion: the exact shift-0
witness `x=(19/20,19/20,1,1,1,1)`, extended by `x_k=1` for `k>=7`, satisfies
every pointwise lower/upper wall and monotone contraction but makes the
normalized degree 7, `m=11` column recurrence negative. Thus even the full
propagated static ratio cone does not imply all column rows. This is an
abstract cone countermodel only, not a zeta-window row or a heat-flow orbit,
and has no edge to `lambda_le_0_goal`.

The former monotone-contraction theorem target is now the validated
parameter-specific ratio-cone theorem:

```text
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_theorem_target.json
```

It records `Delta^3 log A_{k-1}(lambda)>=0`, equivalently
`A_{k+2}A_k^3>=A_{k+1}^3A_{k-1}`, for every index and every finite
`lambda>=-100`. Full cone entry at `lambda=-100` and the infinite heat-flow
maximum principle close the previously live differential-inequality route.
This is an `interval_validated` degree-2 ratio theorem; it still does not imply
the all-shape PF/Jensen bridge or `lambda_le_0_goal` by itself.

The heat-flow monotone-closure scout sharpens the differential-inequality
route inside that target:

```text
outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md
work/rh_compute/results/jensen_window_pf_heat_flow_monotone_closure_scout.json
```

It validates 4 exact lambda-flow algebra rows, 315 finite Arb threshold rows,
and 305 finite Arb flow-bracket rows. It isolates the boundary subtarget
`q >= (2*k-1)/(2*k+5)` for a possible monotone-contraction flow proof. It is a
diagnostic node only: it does not prove a closed differential inequality, does
not prove the monotone-contraction theorem, and has no edge to
`lambda_le_0_goal`.

The heat-flow boundary-threshold lemma discharges that threshold subtarget:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json
```

It proves the stronger exact lower bound `x_k >= (2*k-1)/(2*k+1)` from
Phi positivity and Cauchy-Schwarz raw-moment log-convexity, hence the heat-flow
threshold `x_k >= (2*k-1)/(2*k+5)`. It is an exact lemma node only: it does not
prove adjacent log-concavity, global flow invariance, monotone contractions for
all zeta windows, or `lambda_le_0_goal`.

The heat-flow ratio-cone invariance lemma packages the remaining boundary
algebra into a conditional exact statement:

```text
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json
```

It proves inward-pointing lower-wall, upper-wall, and monotone-wall algebra for
the cone `(2*k-1)/(2*k+1)<=x_k<=1`, `x_{k+1}>=x_k`, assuming a full infinite or
collared finite cone. It is exact but conditional: it does not prove that the
actual zeta coefficients enter the cone at a suitable lambda, does not prove
the infinite flow legitimacy, and has no edge to `lambda_le_0_goal`.

The former heat-flow cone-entry target now records the closed entry-and-flow
theorem:

```text
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
work/rh_compute/results/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.json
```

The repaired finite collar, pointwise Mellin walls, and analytic adjacent tail
prove the full ratio cone at `lambda=-100`; the infinite-sequence maximum
principle then propagates it to every finite `lambda>=-100`. The historical
fixed-k Taylor rows remain diagnostics, while this node is now
`interval_validated`. It has no direct proving edge to `lambda_le_0_goal`.

The Phi Taylor cone-entry sign scout supports only the fixed-k local sign part:

```text
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
work/rh_compute/results/jensen_window_pf_phi_taylor_cone_entry_sign_scout.json
```

It certifies `2*b-a^2<0` and `2*(a^3-3*a*b+3*c)>0` with explicit Phi-tail
bounds. It remains a finite certificate node: it does not supply the
uniform-in-k or collared finite cone-entry theorem and has no edge to
`lambda_le_0_goal`.

The negative-lambda cone-entry prefix scout supports the same target with
actual ACB-enclosed coefficient ratios:

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.json
```

It certifies the finite prefix for `lambda=-25,-50,-100`: lower and upper
walls through `k<=21` and monotone gaps through `k<=20`. It remains finite
evidence only: it does not supply the all-`k` tail theorem or finite-collar
flow theorem and has no edge to `lambda_le_0_goal`.

The negative-lambda finite-collar contract sharpens that prefix:

```text
outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md
work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_contract.json
```

It extracts the collar-ready active depth `K=19`, with first collar `x_20` and
second collar `x_21`. It is a diagnostic accounting node only: it does not
provide the all-`k` tail theorem, does not extend the finite-collar flow beyond
that stated depth, and has no edge to `lambda_le_0_goal`.

The negative-lambda tail-barrier scout sharpens the remaining tail question:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_scout.json
```

It rewrites the tail as defect barriers `0 <= d_k <= 2/(2*k+1)` and
`d_{k+1} <= d_k`, certifies the exact finite cone and defect-monotone rows on
the checked prefix, records that the stronger one-third-width buffer has a
k150 frontier, and rejects the scaled-defect nonincreasing shortcut. It is a
theorem-search diagnostic only: it does not supply the all-`k` tail theorem,
does not prove finite-collar flow beyond the stated depth, and has no edge to
`lambda_le_0_goal`.

The scaled-defect frontier scout sharpens the tail-barrier diagnostic:

```text
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.json
```

It validates the exact scaled cone on all 597 checked k200 rows, while the
fixed half-width buffer holds on only 521/597 rows and first fails at
lambda=-50, k=191 and lambda=-25, k=133. The one-third buffer holds on only
179/597 rows and first fails at lambda=-25, k=31. This blocks promotion of
both fixed buffers as global tail routes, but it does not prove an all-`k`
exact-cone theorem and has no edge to `lambda_le_0_goal`.

The negative-lambda defect-recurrence scout sharpens that diagnostic:

```text
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_defect_recurrence_scout.json
```

It validates a finite-compatible buffered sufficient condition on the older
canonical prefix and rejects the direct width-preserving recurrence on all
checked adjacent rows. The k200 scaled-defect frontier now blocks promotion of
the one-third and fixed half-width buffers as global tail routes, so this remains only a
theorem-search diagnostic with no edge to `lambda_le_0_goal`.

The negative-lambda log-curvature bridge sharpens the buffered route:

```text
outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_log_curvature_bridge.json
```

It rewrites the sufficient condition using `B_k=-Delta^2 log A_k=-log(x_k)`.
The bridge validates 63 simple log-buffer rows and 60 curvature-monotone rows
on the finite prefix, then separates the all-`k` burden into a bounded
log-curvature estimate `B_k<=2/(3*(2*k+1))` and the already-open
`Delta^3 log A` monotone-contraction target. It is a diagnostic/reduction node
only and has no edge to `lambda_le_0_goal`.

The negative-lambda bounded log-curvature target is now historical:

```text
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_target.json
```

It records the formerly proposed quantitative estimate
`B_k=-Delta^2 log A_k<=2/(3*(2*k+1))` and its exact raw moment-curvature
equivalent on the old k<=22 prefix. The repaired k300 obstruction below
finite-rejects that fixed wall, so this node is retained as a diagnostic
record and no longer has `would_imply_if_proved` flow into the negative-lambda
defect-tail target.

The bounded log-curvature k300 obstruction retires that fixed wall:

```text
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md
work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.json
```

It rewrites the old target as `C_k=(2*k+1)*B_k<=2/3`, then classifies the
repaired k300 data: 179/897 rows satisfy the wall, 718/897 rows fail it, and
0 rows are inconclusive. The first checked failure is at `lambda=-25.0,
k=31`; the worst checked slack is at `lambda=-25.0, k=299`. This node has
`blocks_promotion` edges into the retired bounded-curvature target, the
defect-tail route that used it, and the signed-perturbation route unless that
route is revised around the newer raw-corridor or linear curvature-barrier
target.

The negative-lambda Gaussian curvature matrix sharpens that target:

```text
outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.json
```

It compares the actual prefix to the Gaussian raw-moment baseline, where
`x_k=1` and `B_k=0`. The actual finite prefix has a positive but bounded
Gaussian deficit on all checked rows. The matrix also records that positive
Gaussian scale-mixture arguments point the wrong way for the upper wall,
because they push `x_k>=1`. It is a diagnostic node only and has no edge to
`lambda_le_0_goal`.

The negative-lambda signed Gaussian perturbation matrix packages the surviving
fixed-k route:

```text
outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.json
```

It uses the certified Phi Taylor signs to record a positive leading Gaussian
deficit `D2=a^2-2*b` and a positive monotone correction
`2*(a^3-3*a*b+3*c)`. It also records a finite-depth activation scale while
marking it as a scale diagnostic only, since the uniform asymptotic remainder
needed for an all-`k` tail theorem remains open.

The negative-lambda uniform remainder target records that historical gap:

```text
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_uniform_remainder_target.json
```

It records that the fixed-`k` leading deficit has the right sign, but the
former fixed bounded-curvature wall shrinks like `1/k`; the leading-only
cutoff floors for `T=25,50,100` are `0,3,16`, all below the `k>=22` tail
start. The later lambda=-100 saddle, scaled-curvature, raw-corridor, and
adaptive-defect theorems close the destination by an independent route. This
node is therefore historical rather than an active blocker; no simultaneous
uniform remainder theorem is claimed.

The negative-lambda Taylor moment budget sharpens the local/mesoscopic side:

```text
outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md
work/rh_compute/results/jensen_window_pf_negative_lambda_taylor_moment_budget.json
```

It factors the Gaussian moments and makes `q/T`, with `q=k+1/2`, the explicit
local expansion parameter. At the `k=22` tail start, the low-order Taylor
multiplier is an invalid positive-moment model for `T=25,50,100,200`, and is
bounded only in two larger diagnostic samples. This is a diagnostic node only:
it does not replace the actual ACB-certified finite prefix, and it does not
prove the Taylor-tail positivity or log-curvature remainder theorem.

The high-order Taylor scout tests the tempting repair "use more Taylor terms":

```text
outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_high_order_taylor_scout.json
```

It automatically generates the local `P_j(q)` polynomials, certifies even
coefficients through `c14`, and classifies 35 truncation rows for degrees
6,8,10,12,14 at the `k=22` tail start. The matrix still contains invalid
normalizers, upper-wall violations, and overbound rows. This remains a
diagnostic node only; it sharpens the open uniform-remainder target without
proving any infinite Taylor-tail theorem.

The negative-lambda defect-tail theorem target is now a historical
simultaneous-target diagnostic:

```text
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_defect_tail_theorem_target.json
```

Its original all-three-lambda statement requires
`0 <= d_k <= 2/(2*k+1)` for `k>=150` and `d_(k+1)<=d_k` for `k>=149`.
The lambda=-100 theorem below proves both from `k=1`, so the
one-entry-parameter route is discharged. No corresponding all-k theorem at
lambda=-25 or -50 is claimed.

The negative-lambda half-width tail record documents a rejected scaled-defect
shortcut from the k150 frontier:

```text
outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_half_width_tail_target.json
```

The fixed threshold `0 <= s_k <= 1/2` fails on the checked k150 prefix:
only 430/447 rows pass, with first failure at lambda=-25, k=133. This record
blocks promotion of the half-width shortcut; a replacement exact-cone or
adaptive scaled-defect theorem would still need the separate monotone bridge,
and this node has no proving edge to `lambda_le_0_goal`.

The adaptive scaled-defect target records the historical simultaneous
replacement route:

```text
outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.json
```

The original target asks for the exact cone `0<=s_k<=1` or an explicit
adaptive envelope from `k>=200`, plus the separate monotone defect bridge
from `k>=199`, simultaneously on the exploratory lambda grid. The lambda=-100
composition below proves the exact cone and bridge from `k=1`; the stronger
simultaneous formulation is retained as diagnostic history and no longer has
a `would_imply_if_proved` edge.

The adaptive envelope matrix sharpens that replacement route:

```text
outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.json
```

It validates the finite k200 monotone-envelope pattern: 594/594 adjacent
`k`-increase rows and 398/398 cross-lambda order rows, with maximum checked
scaled defect 0.5376643171065356005 at lambda=-25, k=199. This isolates a
possible proof strategy but remains a finite diagnostic; it has no proving
edge to `lambda_le_0_goal`.

The adaptive envelope obligations artifact turns that pattern into exact
ratio inequalities:

```text
outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md
work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.json
```

It records that the lower threshold `s_k<=1` is already the exact
boundary-threshold lemma. The later kernel Mellin upper-wall certificate also
discharges `x_k<=1` for every real lambda and every k. Full cone entry and the
raw corridor now discharge the adjacent wall and scaled-envelope obligations
at lambda=-100. The all-three-lambda version remains historical, and this node
has no proving edge to `lambda_le_0_goal`.

The raw-moment bridge matrix rewrites those obligations in the moment-ratio
language:

```text
outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.json
```

It records the exact identities `x_k=((2*k-1)/(2*k+1))*R_k` and
`0<=s_k<=1 iff 1<=R_k<=(2*k+1)/(2*k-1)`, then turns the monotone bridge plus
scaled `k`-increase into an adaptive corridor for `R_(k+1)`. The k200 prefix
checks 597/597 raw-cone rows and 594/594 corridor rows, while preserving the
fixed-buffer failures. This is a diagnostic edge into the adaptive target and
has no proving edge to `lambda_le_0_goal`.

The raw-ratio decrement-corridor scout sharpens the direct recurrence route:

```text
outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.json
```

It rewrites the adaptive raw corridor exactly as
`2*(R_k-1)/(2*k+1)<=R_k-R_(k+1)<=4*R_k/(2*k+1)^2`. The k200 prefix validates
594/594 decrement-corridor rows, 591/591 theta-k monotonicity rows, and
396/396 theta lambda-order rows. Two exact raw-cone monotone counterexamples
block using raw decrease alone as a recurrence proof, so this remains a
zeta-specific theorem-search diagnostic with no proving edge to
`lambda_le_0_goal`.

The k300 precision-repair audit stress-tests that recurrence route at greater
depth:

```text
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
work/rh_compute/results/jensen_window_pf_negative_lambda_k300_precision_repair_audit.json
```

It records that the broad dps160/cutoff6 k300 run reports lambda=-100 high-k
precision alarms, but local dps220/cutoff7 repairs over k220..250 and
k245..320 restore 897/897 raw wall rows, 894/894 decrement-corridor rows,
891/891 theta-k monotonicity rows, and 596/596 theta lambda-order rows. This
supports the raw-ratio and zeta-specific raw-corridor nodes as finite stress
evidence only; it has no proving edge to `lambda_le_0_goal`.

The raw-log decrement bridge turns the same route into log-ratio language:

```text
outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.json
```

It sets `p_k=log(R_k)` and
`delta_k=p_(k+1)-p_k=log(R_(k+1)/R_k)`, then proves that the raw decrement
corridor is exactly equivalent to
`log(1-4/(2*k+1)^2)<=delta_k<=log(1-2*(1-exp(-p_k))/(2*k+1))`.
The repaired k300 table validates 897/897 raw-log wall rows and 894/894
log-corridor rows, while two exact raw-cone witnesses block raw-log decrease
alone. This sharpens the zeta-specific recurrence target but has no proving
edge to `lambda_le_0_goal`.

The coefficient-curvature corridor bridge refactors the same condition one
more step:

```text
outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.json
```

It sets `B_k=-log(((2*k-1)/(2*k+1))*R_k)` and proves that the raw corridor is
exactly equivalent to
`log((2*k+3)/(2+(2*k+1)*exp(-B_k)))<=B_(k+1)<=B_k`. The upper side is the
existing monotone-contraction/log-curvature target; the lower side is a new
zeta-specific lower curvature barrier. The repaired k300 table validates
897/897 B-wall rows and 894/894 curvature-corridor rows, but monotone
curvature alone is blocked by an exact cone witness, so this remains a
diagnostic node with no proving edge to `lambda_le_0_goal`.

The linear curvature-barrier scout gives a simpler sufficient target for that
lower side:

```text
outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.json
```

It proves the exact calculus lemma
`L_k(B)<=((2*k+1)/(2*k+3))*B` for `B>=0`, so the stronger linear theorem
`B_(k+1)>=((2*k+1)/(2*k+3))*B_k` would imply the nonlinear lower curvature
barrier. The repaired k300 table validates 894/894 linear-barrier rows. The
defect-width recurrence remains rejected, so this is specifically a
log-curvature route, still with no proving edge to `lambda_le_0_goal`.

The scaled-curvature monotonicity target names the open theorem supplied by
that scout:

```text
outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.json
```

It sets `C_k=(2*k+1)*B_k` and asks for `C_(k+1)>=C_k`, exactly equivalent to
`B_(k+1)>=((2*k+1)/(2*k+3))*B_k`. The repaired k300 table validates 894/894
scaled-curvature increase rows, while the retired fixed wall `C_k<=2/3` fails
on 718/897 rows. The target is now discharged by

```text
outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json
```

A `16,074`-block compact interval-Simpson theorem, an exact `u>=5` ray,
the repaired `k<=318` prefix, and the complete-kernel perturbation transfer
prove `C_(k+1)>=C_k` for every integer `k>=1` at `lambda=-100`. The node is
therefore `interval_validated`, not an open blocker.

The scaled-curvature log-ceiling bridge turns that target into raw-log
recurrence form:

```text
outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.json
```

It sets `p_k=log(R_k)` and `delta_k=p_(k+1)-p_k`, then proves that
`C_(k+1)>=C_k` is exactly the affine ceiling
`delta_k<=h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)`. The repaired k300
table validates 894/894 scaled-ceiling rows and 894/894 scaled-log-corridor
rows. The bridge also records that this affine ceiling is stronger than the
nonlinear raw-log upper wall, with a tiny high-k slack, so any promoted proof
needs sharp zeta-specific estimates.

The relative-Gaussian curvature bridge rewrites the same lower-side target as
a weighted four-point inequality for the log moment sequence after subtracting
the Gaussian baseline:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.json
```

It records `B_k=2*f_k-f_(k-1)-f_(k+1)` and rewrites `C_(k+1)>=C_k` as
`(2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)>=0`. The
repaired k300 table validates 897/897 B-positive rows, 894/894 B-decrease
rows, 894/894 C-increase rows, and 598/598 C-lambda-order rows, while keeping
the weighted four-point inequality as an open all-k zeta-specific target.

The relative-Gaussian Taylor stencil scout tests the local Taylor route into
that weighted target:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.json
```

It certifies the fixed-k leading signs for `B_k`, `B_k-B_(k+1)`, and
`C_(k+1)-C_k`, but its 35-row finite truncation matrix has only 4
all-positive stencil rows. This blocks promotion of finite Taylor truncations
and sharpens the uniform remainder target into explicit `B`, companion, and
weighted-gap stencil remainder obligations.

The relative-Gaussian stencil remainder obligations diagnostic decomposes the
Taylor-tail error into exact epsilon stencils:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.json
```

It records exact `E_B`, companion, and weighted-gap error stencils for the
35-row truncation matrix, with 4 positive baseline rows, 31 blocked baseline
rows, 4 exact stencil rows, and 0 ready-to-apply rows. The weakest positive
finite row leaves a half-margin budget of `1.166490564421582442E-8`, so this
is a sharpened uniform-remainder obligation, not a proof of the all-k theorem.

The relative-Gaussian pointwise tail budget converts those stencil margins
into per-index log-tail and multiplicative relative-tail tolerances:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.json
```

It records the exact absolute coefficient sums `4`, `8`, and `16*k+16`, so at
`k=22` the weighted-gap sum is `368`. Among the 4 positive finite baselines,
the companion stencil is limiting in 3 rows and the weighted-gap stencil in 1
row. The current weakest half-safety log-tail envelope is
`1.458113205526978052E-9`, equivalent to the multiplicative relative-tail
ratio bound `1.458113204463930993E-9`. These are required tolerances, not
proved analytic tail estimates.

The relative-Gaussian next-increment stencil stress tests that crude pointwise
budget against the known degree-to-degree Taylor increment:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.json
```

It validates 2 tested next-increment rows and 2 missing degree-16 frontier
rows. The tested increments exceed the pointwise half-safety budget in both
rows, with worst over-budget factor `3.010798908654295615E+3`, but the
structured `B`, companion, and weighted-gap stencil signs are preserved in
both tested rows. This blocks a crude pointwise triangle shortcut while
sharpening the live route toward direct signed stencil-tail estimates.

The relative-Gaussian degree-16 stencil continuation resolves the previous
degree-16 frontier for the degree-14 positive baselines:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.json
```

It validates the negative degree-16 Taylor ratio, tests all 4 current positive
baselines one Taylor step further, and records 3 stencil-sign-preserving rows
and 1 stencil-sign-failure row. The `T=1000` degree-14 baseline fails through
the companion stencil, while the `T=2000` degree-14 baseline survives. This is
a finite collar signal for the direct signed stencil-tail route, not a proof
of an infinite tail theorem.

The relative-Gaussian degree-16 collar scan makes that finite large-T signal
explicit on the integer `T=900..2200` grid:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.json
```

It validates 1301 scan rows. The M=7 baseline first has positive stencils at
`T=918`, the M=8 degree-16 continuation first preserves signs at `T=1156`
with `q/T=1.946366782006920415E-2`, and the stricter half-safety condition
first holds at `T=1483` with `q/T=1.517194875252865813E-2`. The pointwise
budget fails on all 1283 baseline-positive rows, so this node has a
`blocks_promotion` edge back to the pointwise shortcut and a `sharpens` edge
to the open uniform-remainder/direct signed stencil-tail target. It remains a
finite `k=22`, degree-16, integer-grid diagnostic only.

The relative-Gaussian degree-16 real-`T` collar scout upgrades the integer
sign threshold to a continuous finite-surrogate statement:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.json
```

For the rationalized degree-16 fixed-`k=22` surrogate, it certifies positive
normalizers `F_21` through `F_24` and all three structured stencil signs on
the real half-line `T>=1156`. The `B` and companion stencils use stripped
product-polynomial root counts, while the weighted-gap stencil uses positivity
of its derivative numerator. This sharpens the open uniform-remainder target
by isolating the next required upgrade: interval coefficient control and a
signed infinite-tail stencil theorem on the same real collar.

The relative-Gaussian degree-16 Arb real-`T` collar certificate upgrades that
finite surrogate from midpoint-rationalized ratios to coefficient-ratio balls:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.json
```

It uses Bernstein coefficients on `0<=u<=1/1156` to certify the Arb
finite-degree normalizers `F_21` through `F_24`, the stripped `B` product, the
stripped companion product, and the stripped weighted-gap derivative numerator.
The relative-Gaussian degree-40 Arb collar ladder stress then repeats that
same real-`T` collar test for every even finite surrogate degree `16` through
`40`:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.json
```

All tested ladder levels certify positive normalizers and all three structured
stencils with zero Bernstein failures. The residual-tail budget extracts the
next sufficient analytic target from those margins:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json
```

It records that a fixed-`k=22` proof on `0<=u<=1/1156` would be enough if it
proved `|R_i(u)| <= 5.382819486765314521E-01*u^3` and
`|R_i'(u)| <= 9.315354075509573936E-03*u` for `i=21..24`. The finite
degree-42..80 tail profile consumes less than `0.1%` of those budgets, but
this remains only a target: no analytic majorant for the infinite residual
tail is proved. The formal-tail obstruction scout then rejects the naive
infinite formal-tail summation shortcut:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json
```

Through degree 240 the formal terms decrease to a least-term region near
`j=103`, then all four value and derivative rows grow again from `j=103` to
`j=104`. It also records that a fixed-radius Cauchy coefficient bound cannot
be summed term-by-term after multiplication by the moment rising factors. This
blocks promotion of the finite formal-tail profile while sharpening the live
route: prove an actual contour, saddle, integral-remainder, or optimally
truncated asymptotic theorem. The asymptotic-remainder target turns that into
explicit constants:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.json
```

It records that a common `1000x` first-omitted-term theorem would fit inside
the fixed-`k=22` half-safety budgets, with common multiplier limit about
`1419.939`, and that an optimized finite window through `j=120` would leave
much larger least-term slack. This is still a target calibration, not a
remainder theorem. The actual endpoint remainder scout checks the real
relative-Gaussian multiplier at the collar endpoint:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.json
```

At `T=1156`, selected generalized Gauss-Laguerre order `N=320` gives value
and derivative residuals for `F_21` through `F_24` below one first omitted
term and below `0.1%` of the degree-40 half-safety budgets. The graph records
this as floating endpoint support for the asymptotic target only; it blocks
promotion to an infinite formal-tail proof and does not supply an
interval-certified or uniform collar remainder theorem. The
cancellation-reduced grid scout then subtracts the degree-40 polynomial inside
the Gamma expectation and tests a wider finite `T` grid:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
```

On `T=1156,1500,2000,5000,10000` and `F_21` through `F_24`, all sampled
value and derivative residuals are below one first omitted formal term across
four quadrature orders. The worst sampled value ratio is about
`0.9707100590`, and the worst sampled derivative ratio is about
`0.9693567775`, both at `T=10000`, `F_21`. This repairs the far-`T`
cancellation failure of separate floating subtraction as evidence, but it is
still floating-node finite-grid evidence rather than an interval-certified or
uniform collar theorem. The intervalization target converts that finite-grid
slack into a certification checklist:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
```

It records that the common slack to one first omitted term is about
`2.928994097967e-2`, so a future interval certificate with total ratio error
below `1.0e-2` would keep the finite grid below the first-omitted threshold.
The open requirements are Laguerre node/weight intervals, `Phi` and `Phi'`
n-tail bounds, quadrature-remainder bounds, coefficient-ratio propagation,
rounding aggregation, and a separate grid-to-collar bridge. These finite collar
nodes still have no edge to
`lambda_le_0_goal`; the residual-tail theorem and the all-`k` upgrade remain
open.

The Phi tail bound scout addresses one of those intervalization sources:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json
```

On the padded range `0<=x<=1`, the `n>30` tails in `Phi`, `Phi'`, and
`Phi(0)` are bounded far below the `2.0e-3` per-source intervalization cap
after normalizing by the diagnostic `Phi(0)>=0.44` proxy. This is conditional:
the graph keeps the node as a diagnostic that still requires interval
Laguerre-node range and `Phi(0)` lower-bound certificates, and it has no edge
to `lambda_le_0_goal`.

The node-c0 range certificate supplies those two finite-grid side conditions:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json
```

Gershgorin applied to the Laguerre Jacobi matrix, with AM-GM row bounds,
certifies every recorded Laguerre node below `809<T_min=1156`, so the padded
range condition `x<=1` holds. An Arb calculation certifies the `n=1` term of
`Phi(0)` above `0.44`. This supports the Phi-tail scout only; individual
Laguerre weights, quadrature error, coefficient propagation, rounding, and
grid-to-collar coverage remain open.

The Phi-tail grid certificate composes those two artifacts:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
```

It certifies the omitted `n>30` Phi, Phi', and Phi(0) tail source for the
recorded finite grid below the `2.0e-3` per-source intervalization cap. The
graph records this as source-level support only; finite `n<=30` node
evaluation, Laguerre weights, quadrature error, coefficient propagation,
rounding, and grid-to-collar coverage remain open.

The quadrature ladder scout targets the next finite-grid numerical weakness:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json
```

It recomputes the worst cancellation-reduced row `T=10000`, `F_21` at
floating generalized Gauss-Laguerre orders `96` through `320`. The largest
value ratio remains about `0.970710059020335`, the largest derivative ratio
about `0.969356777475809`, and the order spread is below `1e-14`. The graph
records this only as calibration for a future rigorous quadrature-radius
target, not as a quadrature-remainder theorem.

The quadrature-remainder route matrix turns that calibration into explicit
proof obligations:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json
```

It records the classical `N=320`, `alpha=41/2` Gauss-Laguerre remainder
prefactor, about `2.79e-159`, and converts the `1e-6` quadrature
ratio-radius target into concrete 640th-derivative supremum caps or,
alternatively, an interval adaptive integration target with a separate far
Gamma-tail bound. It also checks that adding the `1e-6` quadrature cap to the
worst-row finite-plus-tail budget still leaves both ratios below one first
omitted term. The graph records this as a route matrix only; no derivative
bound, interval integration certificate, or quadrature-remainder theorem is
proved, and it has no edge to `lambda_le_0_goal`.

The worst-row far-tail split certificate retires the first interval-integration
subpiece:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json
```

It splits the worst-row normalized Gamma expectation at `y=200`, proves with
Arb margins that the finite `n<=30` Phi and derivative-core majorants are
decreasing beyond the split, and bounds the polynomial part by upper
incomplete-Gamma moments. The resulting value and derivative tail ratios are
about `8.66e-26` and `4.50e-27`, far below the `1e-6` quadrature-radius
target. The graph records this as far-tail support only; compact integration
on `0<=y<=200`, aggregation, all-row coverage, and grid-to-collar coverage
remain open, and it has no edge to `lambda_le_0_goal`.

The compact-interval integration scout tests the next subpiece:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.json
```

It imports the far-tail split and evaluates six raw Arb panel hulls on
`0<=y<=200`. The raw value and derivative width-to-cap ratios are about
`1.38e39` and `4.28e37`, so plain interval-Riemann hulls are explicitly
rejected. The graph records this as a route-sharpening diagnostic only:
Taylor/Chebyshev panel models with exact Gamma-weighted moments remain open,
and it has no edge to `lambda_le_0_goal`.

The Chebyshev panel-moment scout calibrates that local-model route:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.json
```

It integrates high-precision floating Chebyshev interpolants against
incomplete-Gamma panel moments on the same six compact panels. Consecutive
degree deltas from `16->20` onward are below the unscaled quadrature caps in
both channels, with reference degree-32 estimates about `-6.58e-34` and
`-1.38e-32`. The graph records this as calibration only: Arb coefficient
balls and interpolation-remainder bounds remain open, and it has no edge to
`lambda_le_0_goal`.

The Arb Chebyshev interpolant-moment scout promotes that calibration into
ball arithmetic for the interpolants:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json
```

It Arb-encloses the Chebyshev interpolant coefficients and incomplete-Gamma
panel moments for degrees `16`, `20`, `24`, and `32`. All consecutive
interpolant Cauchy deltas are below the unscaled quadrature caps, but the
graph records this as interpolant arithmetic only: the true-function
interpolation remainder on each panel remains open, and it has no edge to
`lambda_le_0_goal`.

The interpolation-remainder route matrix turns that open line into explicit
Bernstein/Taylor-model proof obligations:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json
```

It records the six compact-panel Gamma masses, with `20<=y<=50` the heaviest
panel at mass upper about `0.6021615933248953`, and sizes 20
degree/rho Bernstein budgets plus 16 minimal-degree rows. At `rho=2`, degree
`128` allows only a value sup-norm budget about `1.60e-2` on the heaviest
panel, while degree `160` raises that budget above `6.85e7`. The graph records
this as route sizing only: analytic-domain, sup-norm, endpoint, and true
interpolation-remainder certificates remain open, and it has no edge to
`lambda_le_0_goal`.

The endpoint parity-repair matrix isolates the branch issue on the first
panel:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json
```

It expands the finite `n<=30` Phi truncation at `x=0` with Arb power series
through odd degree `15`. The odd coefficients are tail-sized rather than
structurally zero: the first odd coefficient has absolute upper about
`1.50e-1300`, and the largest recorded low-order odd coefficient is degree
`15` with upper about `1.54e-1255`. The graph records this as an endpoint
repair obligation only: a valid compact proof must either use exact evenness
of the infinite theta kernel plus a certified tail charge, or certify the
first panel in the `x` variable. Near-evenness is not promoted to endpoint
analyticity, and the node has no edge to `lambda_le_0_goal`.

The endpoint x-panel route matrix quantifies that second repair route:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json
```

It records the change of variables `y=T*x^2`, mapping `0<=y<=1` to
`0<=x<=0.01`, with transformed Gamma density power `x^42` for `alpha=20.5`.
The first-panel mass upper is about `1.62e-21`, only about `2.68e-21` of the
heaviest panel mass, and the matrix records 18 first-panel Bernstein budgets.
The graph records this as endpoint route sizing only: exact transformed
moments, an x-domain/sup-norm theorem, and the true x-panel interpolation
remainder remain open, and the node has no edge to `lambda_le_0_goal`.

The endpoint x-moment Taylor certificate discharges that local finite-core
obligation for the worst row:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json
```

For `T=10000`, `F_21`, it expands the finite `n<=30`
cancellation-reduced cores through degree `64`, integrates all `65` monomial
rows by exact transformed Gamma moments, and bounds the post-degree-64 tail
using a rigorous `|z|<=0.2` complex-disk Cauchy majorant. The first-panel
value and derivative balls are certified negative, and each absolute-to-cap
ratio is below `4.10e-47`. The graph records this as a genuine finite
certificate for `0<=y<=1`, but only for the finite core: the separate `n>30`
tail composition, five panels with `y>=1`, all-row coverage, and uniform
collar theorem remain open. The node has no edge to `lambda_le_0_goal`.

The compact x-moment Taylor certificate then extends the same method across
all six compact panels at once:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.json
```

For `T=10000`, `F_21`, the full range `0<=y<=200` maps into
`x<=sqrt(0.02)<0.38`. A degree-`128` Arb Taylor model is integrated by `129`
exact transformed Gamma moments, and the rigorous complex-disk remainder
radii consume less than `4.0e-10` and `1.3e-9` of the value and derivative
caps. Both compact integrals are certified negative. This removes the
six-panel interpolation-remainder obligation for the finite `n<=30` core;
the graph now records the compact source as finite-validated rather than open.
The separate `n>30` normalization/tail splice, `y>200` composition, all-row
coverage, and uniform collar theorem remain open, and the node has no edge to
`lambda_le_0_goal`.

The worst-row full-expectation certificate composes the remaining sources for
that row:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.json
```

It combines the compact degree-`128` certificate, the finite-core `y>=200`
tail, and an Arb global `n>=31` normalization correction. Explicit positive
monotonicity margins extend the omitted-`n` bounds from `0<=x<=1` to every
`x>=0`, including the `Phi(0)` denominator correction. The complete true
value and derivative expectations are negative, with first-omitted ratio
uppers about `0.9707101` and `0.9693568`; both margins below one exceed
`0.029`. No integral source remains open for `T=10000`, `F_21`, and the final
row proof does not use generalized Gauss-Laguerre quadrature. The graph still
records this as one-row finite evidence: the other grid rows, all-source
aggregation, and finite-grid-to-collar theorem remain open, and the node has
no edge to `lambda_le_0_goal`.

The all-row direct expectation certificate generalizes that construction to
the complete recorded grid:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.json
```

For each of the five recorded `T` values and `F_21` through `F_24`, it splits
at `y=min(200,(64/625)T)`, keeping the exact-moment core inside
`x<=8/25<0.38`. A degree-`384` Cauchy certificate, rowwise monotone real-tail
bounds, and the global `n>=31` normalization correction certify all `20`
value and all `20` derivative balls negative and below one first omitted
term. The worst ratios remain at `T=10000`, `F_21`, below `0.970711` and
`0.969357`. This retires every direct integration source on the recorded
grid and excludes floating quadrature as a proof input. The graph preserves
the remaining signed-stencil aggregation and real-`T` collar obligations;
the node has no edge to `lambda_le_0_goal`.

The recorded-grid stencil composition certificate then feeds those residual
balls into the degree-40 fixed-`k=22` inequalities:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.json
```

It uses the exact rational budgets `A=1/2` and `B=9/1000`, rather than the
older floating threshold decimals, and recomputes all perturbation bounds in
Arb. The 20 residual rows use less than `5.7e-4` of either budget. Positive
retained margins certify all four normalizers, the B product, companion
product, and weighted-gap derivative at each of the five recorded `T`
values. The graph records this as a fixed-`k` finite certificate only: no
interval in `T`, all-`k` theorem, or cone-entry bridge follows, and the node
has no edge to `lambda_le_0_goal`.

The finite-collar-segment certificate removes the gaps between those five T
values:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.json
```

A two-regime exact-moment/Cauchy and incomplete-Gamma majorant proves the
same rational residual budgets for every real `1156<=T<=10000`. The worst
value and derivative fractions are below `1.003e-3` and `1.031e-3` of their
budgets. Composing with the four positive Arb perturbation margins certifies
the complete fixed-`k=22` stencil system on that bounded interval. This is an
interval certificate, not sampled evidence; it still has no edge to
`lambda_le_0_goal`.

The full-kernel evenness/Cauchy lemma supplies the exact asymptotic factor
needed beyond the bounded interval:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.json
```

The Jacobi theta transformation and the covariance of
`L_x=2*x^2*d_x^2+3*x*d_x` prove the full infinite kernel exactly even. On
`|z|<pi/8`, subtracting the normalized even Taylor terms through degree `40`
therefore leaves an order-`42` zero and explicit `x^42` Cauchy factors. The
graph distinguishes this exact lemma from the rejected inference based on
finite-truncation near-evenness.

The kernel Mellin upper-wall certificate uses that exact evenness to remove
one full all-k cone wall:

```text
outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md
work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json
```

On `0<=y<=1/25`, 200 Arb subintervals and a full-kernel Cauchy tail prove
`g'(y)^2-g(y)g''(y)>0` for `g(y)=Phi(sqrt(y))`; the minimum outward lower
bound is above `0.978`. On `u>=1/5`, an analytic `n=1`-dominant estimate
retains a `D=L'-uL''` lower margin above `5.6118` after the complete
`n>=2` perturbation. Hence `g` is strictly log-concave on the half-line.
The classical Berwald-Borell Gamma-normalized Mellin theorem then proves
`x_k(lambda)<=1` for every real lambda and every `k>=1`. Together with the
separate Cauchy-Schwarz lemma, both pointwise cone walls are now exact; only
the adjacent-k wall `x_(k+1)>=x_k` remains. The certificate sharpens the
cone-entry, defect-tail, and raw-corridor targets and has no edge to
`lambda_le_0_goal`.

The log-concave Mellin monotone-wall countermodel prevents overextending that
theorem:

```text
outputs/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.md
work/rh_compute/results/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.json
```

For the log-concave density `exp(-5*y)*1_[0,1](y)`, the normalized Mellin
transform has the exact incomplete-Gamma form
`gamma_lower(p,5)/(5^p*Gamma(p))`. Arb certifies both `x_1` and `x_2` inside
the Berwald-Borell upper wall but `x_2-x_1<-0.027`. Thus generic
log-concavity does not imply `x_(k+1)>=x_k`; the remaining theorem must use
zeta-specific structure. This is an abstract countermodel gate, not a
counterexample to the zeta sequence, and it has no edge to
`lambda_le_0_goal`.

The T=1156 monotone-wall certificate supplies the corresponding actual-kernel
promotion gate:

```text
outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.json
jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate
```

Rigorous ACB integration and analytic n-series/u-tail bounds enclose
`A_119..A_122` at `lambda=-1156`. Arb then certifies
`x_121-x_120<-1.68e-8`, while both contractions remain strictly between zero
and one. Thus the fixed-`k=22`, `T>=1156` theorem cannot be promoted to an
all-`k` cone theorem at `T=1156`, and the universal all-real-lambda
monotone-contraction formulation is false. The graph redirects cone entry to
a moderate fixed negative parameter with a deep finite collar plus a
zeta-specific eventual-`k` saddle theorem. This counterexample gate has no
edge to `lambda_le_0_goal`.

The kernel summand-shift lemma starts the revised fixed-parameter tail:

```text
outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md
work/rh_compute/results/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.json
jensen_window_pf_negative_lambda_kernel_summand_shift_lemma
```

For `a_n=(log n)/2`, exact exponent bookkeeping gives
`phi_n(u)=n^(-1/2)*phi_1(u+a_n)`. After shifting the moment integral, the
relative integrand is
`rho=n^(-1/2)*(1-a_n/v)^(2k)*exp(T*(2*a_n*v-a_n^2))`; it increases in `v`
and decreases in `k`. At `T=100`, `k>=300`, and `v<=3/2`, Arb certifies the
sum over `n=2..20` below `2.122e-29` of the first-summand moment, while
`a_21>3/2` puts all `n>=21` beyond the split. This sharpens the tail target
to the dominant `n=1` adjacent wall and a shifted `v>=3/2` far-tail budget.
The later dominance theorem closes that far-tail source. The shift lemma has
no edge to `lambda_le_0_goal`.

The first-summand dominance certificate completes the all-n reduction:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json
jensen_window_pf_negative_lambda_first_summand_dominance_certificate
```

In the original `u` variable, every tail ratio satisfies an exact negative
logarithmic derivative. An adaptive split `a(k)=log(k)/8`, strict concavity
of the first-summand tilted integrand, and 15 Arb-positive endpoint and
half-line propagation gates prove

```text
M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6, k>=300,
|L_k-L_k^(1)|<=16/(k-1)^6,              k>=301.
```

Thus the full shifted far tail is no longer open. This interval theorem
sharpens the cone-entry and raw-corridor targets but has no edge to
`lambda_le_0_goal`.

The repaired lambda=-100 source also supplies a deeper finite collar:

```text
outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.json
jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate
```

Arb certifies 76 positive coefficients `A_245..A_320`, both pointwise cone
walls for `x_246..x_319`, and all adjacent walls for `k=246..318`. Nineteen
previously unused rows extend the finite handoff through `k=318`; the weakest
new log gap remains above `3.709e-6`. This finite certificate supports the
open tail target only and has no edge to `lambda_le_0_goal`.

The exact cumulant bridge sharpens the analytic tail obligation:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json
jensen_window_pf_negative_lambda_first_summand_cumulant_bridge
```

For the continuous first-summand log moment `F(t)`, it proves

```text
F'''(t)=kappa_3,t(2*log(U)),
Delta^3 F(k-1)=integral_[0,1]^3 F'''(k-1+r+s+v)drdsdv.
```

The Gamma contribution is exact, and coefficient positivity after
`k=319+n` proves that the single bound

```text
kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318,
```

would imply `L_k^(1)>=1/(4*k^2)` for all `k>=319`. The bridge itself is
exact. The later compact-plus-ray paired certificates discharge the displayed
cumulant estimate; the historical finite floating samples are not used as
the proof and the node has no proving edge to `lambda_le_0_goal`.

The leading part of that cumulant estimate is now certified:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json
jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate
```

After `x=2*log(u)`, symbolic differentiation through order seven, `40,740`
Arb compact subintervals, and exact `u>=5` ray gates prove caps `13/20`,
`1/100`, and `1/1000` through fifth order.

```text
t^2*V'''(x_t)/V''(x_t)^3<=13/20, t>=318.
```

Consequently the open cumulant floor is reduced to the signed normalized
remainder

```text
seventh-order normalized saddle remainder>=-79/1000.
```

The interval theorem sharpens the open target but does not promote its nine
floating remainder samples and has no proving edge to `lambda_le_0_goal`.

The paired compact remainder is now certified:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json
jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate
```

Exact paired standardized moments, `40,736` Arb eighth-derivative envelope
intervals, rigorous midpoint errors, and explicit two-sided tails prove

```text
seventh-order normalized saddle remainder>=-79/1000,
0.9264<=u<=5.
```

The `4,074`-block cover reaches `t=V'(x(5))>1.5241916613e10`. The same
outward-rounded blocks certify strict negative third log-cumulant throughout
the compact range; the closest upper endpoint is below `-9.13e-6`.

The remaining ray is now proved analytically:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json
jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate
```

An adaptive `sqrt(8 log q)` window, first-order Gaussian moment comparison,
and explicit tails prove `H_t>=-3/250` on `u>=5`. Together with the compact
theorem this proves the required cumulant bound for every real `t>=318`.
This node closes the adjacent-wall target but has no proving edge to
`lambda_le_0_goal`.

The resulting lambda=-100 wall closure is recorded as:

```text
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.json
jensen_window_pf_negative_lambda_first_summand_saddle_wall_target
```

Exact geometry gives a unique first-summand saddle in
`(log(k)+22/25)/8 < s_k < log(k)/4`. Nine 60-digit saddle-centered samples
through `k=20000` satisfy `L_k^(1)>1/(4*k^2)` and have increasing
`k^2*L_k^(1)`, but these are floating diagnostics. The analytic certificates
and exact cumulant bridge now prove

```text
L_k^(1)>=1/(4*k^2), k>=319,
```

which dominates the certified `16/(k-1)^6` perturbation and splices to the
finite collar. This node is `interval_validated` but has no proving edge to
`lambda_le_0_goal`; the downstream cone/all-order bridge remains open.

Full cone entry at lambda=-100 is now certified:

```text
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json
jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate
```

A precedence-merged repaired source proves `321` coefficient signs, `319`
pointwise cone rows, and `318` adjacent prefix rows. The exact all-k lower
and upper walls plus the analytic adjacent tail from `k=319` prove the full
infinite ratio cone at `lambda=-100`. This discharges cone entry; forward-flow
legitimacy is supplied by the following certificate, while the all-order
bridge remains open.

Infinite forward invariance is certified in:

```text
outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md
work/rh_compute/results/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.json
jensen_window_pf_heat_flow_infinite_cone_invariance_certificate
```

For the defects `d_k=1-x_k` and increments `h_k=x_{k+1}-x_k`, the exact heat
ODE factors as a cooperative nearest-neighbour term, a bounded zeroth-order
term, and a nonnegative source. After an exponential transform, the pointwise
walls force any negative spatial infimum on a compact heat interval to be
attained at a finite index, where the first-crossing derivative is
nonnegative. Hence the full ratio cone propagates from `lambda=-100` to every
finite `lambda>=-100`, including `lambda=0`. This closes the flow handoff but
does not prove PF-infinity, the all-shape Jensen bridge, RH, or
`lambda_le_0_goal`.

The first post-cone higher-difference diagnostic is:

```text
outputs/jensen_window_pf_defect_complete_monotonicity_scout.md
work/rh_compute/results/jensen_window_pf_defect_complete_monotonicity_scout.json
jensen_window_pf_defect_complete_monotonicity_scout
```

Arb arithmetic certifies `3284` alternating defect-difference intervals and
`3288` negative-log-contraction differences on five cached nonnegative heat
times, with both channels complete through order 8. It keeps `838`
order-9-through-12 intervals as inconclusive. The exact
Hausdorff-defect sequence `x=(1/2,1,1,...)` satisfies the full static ratio
cone and complete defect monotonicity but gives degree-3 Jensen discriminant
`-27/16`; therefore this pattern cannot be promoted directly to the all-shape
bridge.

The multiplier-specific high-precision pass resolves those finite
inconclusives for the negative-log channel:

```text
outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.json
jensen_window_pf_multiplier_complete_monotonicity_frontier_scout
```

Using the dps220 coefficient sources through `A_57` and 250-digit derived Arb
arithmetic, all 7980 intervals in the complete order-55 difference triangle
are strictly positive on the five cached nonnegative heat parameters. This
sharpens the counting-measure target but does not construct its required
unit-atomic measure or promote finite complete monotonicity to PF-infinity.

The exact Hausdorff bridge converts that finite pattern into the correct
infinite recovery problem:

```text
outputs/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.md
work/rh_compute/results/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.json
jensen_window_pf_multiplier_hausdorff_uniqueness_bridge
```

It proves that the complete integer sequence `y_k=-log(x_k)` determines a
unique finite measure `nu` on `[0,1]`. The multiplier product exists exactly
when `dnu/dr=q(-log r)*sum_j r^alpha_j` for a unit counting measure, and the
multiset is then unique. The bridge also retains the exact periodic-factor
guard `sin(2*pi*s)`, so the continuous Mellin obstruction cannot be promoted
without an additional canonical-interpolation or growth theorem.

The first conditional atom constraints are interval-certified:

```text
outputs/jensen_window_pf_multiplier_leading_atom_bound_certificate.md
work/rh_compute/results/jensen_window_pf_multiplier_leading_atom_bound_certificate.json
jensen_window_pf_multiplier_leading_atom_bound_certificate
```

If the unit counting representation exists, the order-6 atom-kernel root lies
in `(4.863538496,4.863538497)`, so its smallest atom exceeds
`4.863538496`. All other available orders through 55 are weaker. The order-3
count inequality also gives `N(11/2)<=1`. These are conditional constraints;
they do not by themselves prove that an atom occurs below `11/2` or construct
the measure. The later ratio composition uses the same cutoff to reject the
unit-atomic measure entirely.

The exact heat-flow hierarchy exposed by that guard is:

```text
outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md
work/rh_compute/results/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.json
jensen_window_pf_heat_flow_jensen_hierarchy_lemma
```

It proves `J_(d+1,n)=J_(d,n)+z*J_(d,n+1)`,
`partial_z J_(d,n)=d*J_(d-1,n+1)`, and
`partial_lambda J_(d,n)=((4*n+2)*partial_z+4*z*partial_z^2)*J_(d+1,n)/(d+1)`.
For the one-atom Hausdorff defects `d_k=(4/9)(9/25)^(k-1)`, the exact
full-cone cubic boundary is `x=5/9,y=21/25,z=589/625`, and the hierarchy gives
`(partial_lambda F)/r_0=329728/2109375>0`, so the discriminant derivative is
negative. Thus even complete defect monotonicity plus the local heat ODE does
not provide the needed coupled higher-minor invariant. No such invariant is
proved, and the node only sharpens the open `target_jensen_window_pf_bridge`.

The cubic hierarchy now has an exact reciprocal-defect cone:

```text
outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md
work/rh_compute/results/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.json
jensen_window_pf_cubic_reciprocal_defect_invariance_lemma
```

Writing `q_k=(1-x_k)^(-1/2)`, cubic hyperbolicity is exactly
`q_(k+1)-q_k<=1`. At a saturated boundary the heat field points inward iff
the next shift obeys the same inequality. This gives a conditional infinite
maximum principle and certifies 318 lambda=-100 prefix margins plus 310
nonnegative-grid margins. The next certificate closes its all-k lambda=-100
tail, and the subsequent weighted-flow certificate closes forward-uniformity.
All higher minor cones remain open.

The cubic entry theorem now holds at every shift at lambda=-100:

```text
outputs/jensen_window_pf_cubic_m100_tail_entry_certificate.md
work/rh_compute/results/jensen_window_pf_cubic_m100_tail_entry_certificate.json
jensen_window_pf_cubic_m100_tail_entry_certificate
```

Compact-plus-ray negative skewness gives a two-sided adjacent log wall for
every `k>=319`. Telescoping yields `d_m>=1/(5*m+1)`, and exact positive-
coefficient polynomial checks prove `q_(k+1)-q_k<1` for every `k>=319`, with
the increment tending to zero. Spliced to the 318 Arb prefix margins, this
proves every shifted cubic Jensen polynomial is hyperbolic at `lambda=-100`.
The next certificate closes the spatial tail uniformly on each finite forward
heat interval; no higher-degree minor cone is claimed here.

The cubic cone now propagates through lambda=0:

```text
outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md
work/rh_compute/results/jensen_window_pf_cubic_forward_uniform_tail_certificate.json
jensen_window_pf_cubic_forward_uniform_tail_certificate
```

For `g_k=q_(k+1)-q_k`, the exact increment field is increasing in `g_(k+1)`.
At a nonincreasing cubic barrier it obeys `sqrt(k)*g_k'<=7*r_k`. The entry
estimate `sqrt(k)*g_k<12` and the explicit coercive supersolution

```text
G_k^epsilon=(12+7*R_L*t)/sqrt(k)+epsilon*k*exp(5*R_L*t)
```

prove `sup_[-100,L] g_k=O_L(k^(-1/2))` for every finite `L`. This legitimizes
the infinite first-crossing principle and proves every shifted degree-3 Jensen
polynomial hyperbolic at `lambda=0`. The node only sharpens
`target_jensen_window_pf_bridge`: degree four and every all-degree conclusion
remain open.

The unchanged cubic cone cannot control the quartic frontier:

```text
outputs/jensen_window_pf_quartic_boundary_flow_obstruction.md
work/rh_compute/results/jensen_window_pf_quartic_boundary_flow_obstruction.json
jensen_window_pf_quartic_boundary_flow_obstruction
```

The normalized discriminant is `Disc(J_4)=256*x^6*y^2*Q(x,y,z)`. The exact
rational boundary polynomial

```text
(1+13*w/20)^2*(1+21*w/50)*(1+57*w/25)
```

has nondecreasing contractions satisfying every local pointwise wall and three
strict neighboring cubic inequalities. Nevertheless,

```text
Q'/r_1=-13108711376416987159336748097/
       20606742971316325673502124987495<0.
```

Thus the cubic invariant does not promote unchanged to degree four. This is an
abstract local countermodel, not a failure of the zeta trajectory; it blocks
one proof shortcut and requires a new coupled quartic condition.

That coupled quartic boundary condition is now exact locally:

```text
outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md
work/rh_compute/results/jensen_window_pf_quartic_double_root_threshold_lemma.json
jensen_window_pf_quartic_double_root_threshold_lemma
```

For `P=(1+a*w)^2*(1+b*w)*(1+c*w)`, `2*a+b+c=4`, and `p=b*c`, define

```text
U=(-a^2+2*a+p)*(3*a^2-5*a+5*p)/(6*p^2).
```

The normalized heat value at the double root has sign `u-U`, while
`P''(-1/a)=2*(a-b)*(a-c)`. Hence the exact inward condition is

```text
(3*a^2-4*a+p)*(u-U)<=0.
```

On the triple-root stratum, first-order viability requires `u=U`; there the
first variation retains `(1+a*w)^2`. This sharpens the quartic handoff but does
not construct a closed contraction-coordinate cone or prove higher-time
tangency.

The quartic threshold is exactly adjacent-degree polar contact:

```text
outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md
work/rh_compute/results/jensen_window_pf_quartic_quintic_polar_contact_lemma.json
jensen_window_pf_quartic_quintic_polar_contact_lemma
```

Normalized adjacent windows satisfy

```text
P_d=P_(d+1)-w*P_(d+1)'/(d+1).
```

For a positive-root hyperbolic quintic, a nonroot polar zero is simple and a
shared nonzero root loses exactly one unit of multiplicity. Thus a quartic
double root forces a quintic triple root. In the quartic boundary coordinates,
`P_5(-1/a)` is a nonzero prefactor times `-(u-U)`, so the heat threshold
`u=U` is precisely quintic triple contact. The identity extends to every
adjacent degree, but the resulting top-down infinite hierarchy has no
noncircular compactness or limiting closure yet.

The infinite hierarchy now has an exact cofinal-degree reduction:

```text
outputs/jensen_window_pf_cofinal_degree_polar_closure_lemma.md
work/rh_compute/results/jensen_window_pf_cofinal_degree_polar_closure_lemma.json
jensen_window_pf_cofinal_degree_polar_closure_lemma
```

Repeated polar descent proves that one hyperbolic terminal degree `D` closes
every lower degree at the same shift. Therefore, for fixed `n`, an unbounded
set of hyperbolic terminal degrees is sufficient for every finite degree.
Applying this at every shift would close the complete Jensen family.

The evidence audit keeps the distinction sharp: 1,050 Sturm rows certify
degrees 3 through 12 on shifts 0 through 20 and five nonnegative heat values.
The separate 2,875 degree-3 through degree-12 rows certify contraction
inequalities. Neither bounded dataset is a cofinal terminal sequence.

The cofinal route has an exact endpoint-equivalence guard:

```text
outputs/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.md
work/rh_compute/results/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.json
jensen_window_pf_cofinal_scaling_limit_equivalence_gate
```

For fixed shift `n`, `P_(D,n)(z/D)` converges locally uniformly to the entire
function `F_n(z)=sum_j A_(n+j)z^j/j!`. Thus any unbounded hyperbolic degree
subsequence puts `F_n` in the Laguerre-Polya class, and Jensen's theorem gives
the converse. At shift zero the cofinal terminal theorem is already endpoint-
equivalent; it is not supplied by familiar fixed-degree/large-shift asymptotics
or by the bounded degree-3 through degree-12 ladder.

The polar hierarchy now identifies how a genuine boundary avoids every fixed
degree:

```text
outputs/jensen_window_pf_polar_heat_collision_cascade_lemma.md
work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json
jensen_window_pf_polar_heat_collision_cascade_lemma
```

At a multiplicity-`m` root of `J_d`, the adjacent value `J_(d+1)(xi)`
controls the complete low heat jet. If the whole upper tower is hyperbolic,
polar multiplicity lifting creates a fixed-root cascade whose multiplicity
grows once per degree. Scaling the cascade proves that the entire coefficient
function must be `exp(-z/xi)` times a bounded-degree polynomial. Since the Xi
coefficient function is not of that form, every fixed Jensen degree is strict
at an LP endpoint; any loss of hyperbolicity must escape to `D->infinity`.

That escaping layer has an exact first asymptotic model:

```text
outputs/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.md
work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json
jensen_window_pf_scaled_double_zero_boundary_layer_lemma
```

The scaled Jensen expansion begins with `F(z)-z^2*F''(z)/(2D)` and has second
correction `(z^3*F'''/3+z^4*F''''/8)/D^2`. Near a nondegenerate double zero `rho<0`, under
`lambda=lambda_*+tau/D` and `z=rho+eta/sqrt(D)`, the local polynomial tends to
`eta^2+8*rho*tau-rho^2`. Hence the finite pair is still split by an unscaled
gap of order `D^(-3/2)` at `lambda_*`, and its collision is shifted to
`lambda_*+rho/(8D)+O(D^(-2))`. The exact `D^(-2)` term contains
`F'''(rho)/F''(rho)=3U'(rho)/U(rho)` after writing
`F=(z-rho)^2 U`, so the global root external field first appears there. Nested
finite-degree eventual thresholds exhaust the Newman boundary. This sharpens the live obligation to a degree- and
zero-height-uniform remainder estimate; neither theorem supplies that estimate
or an edge to `lambda_le_0_goal`.

The global correction now has an exact root-field interpretation:

```text
outputs/jensen_window_pf_newman_root_external_field_lemma.md
work/rh_compute/results/jensen_window_pf_newman_root_external_field_lemma.json
jensen_window_pf_newman_root_external_field_lemma
```

At a squared double root `rho=-x^2`, the regularized field is
`E_x=sum_(j!=*)m_j/(x_j^2-x^2)`, while the signed stiffness is the positive
sum `K_x=sum_y 1/(x-y)^2`. The pair equations become `q'=8-4qS` and
`q=8*(t-t_*)-16*K_x*(t-t_*)^2+...`. Finite positive-coefficient LP products
realize both signs of `E_x`, so this node sharpens the bridge to an Xi-specific
paired-cancellation and far-tail theorem; it supplies no edge to
`lambda_le_0_goal`.

That field is now compared with the classical zero-density reference:

```text
outputs/jensen_window_pf_newman_classical_field_balance_gate.md
work/rh_compute/results/jensen_window_pf_newman_classical_field_balance_gate.json
jensen_window_pf_newman_classical_field_balance_gate
```

The exact continuum field tends to `-pi/8`, matching the `-pi/4` quantile
velocity at high positive time. The published fixed-time asymptotics make all
zeros above `exp(C/t)` simple, so a hypothetical positive-boundary collision is
confined to a `t`-dependent compact region. However, exact even-lattice
deformations move every zero by less than one while forcing the collision field
to either sign of infinity. The graph therefore blocks macroscopic
Riemann-von Mangoldt control as the missing estimate and retains the
lambda-uniform weighted reciprocal-gap bound as an open handoff.

The counting discrepancy has now been integrated through the field kernel:

```text
outputs/jensen_window_pf_newman_local_odd_count_reduction_lemma.md
work/rh_compute/results/jensen_window_pf_newman_local_odd_count_reduction_lemma.json
jensen_window_pf_newman_local_odd_count_reduction_lemma
```

The published uniform `O(log X)` zero-count error controls every contribution
outside `|y-c|<log(4c)^2`. The remaining high collision field is
`-pi/8+S_H+O(1/log c)`, where `S_H` is exactly the inverse-square weighted odd
local counting discrepancy. This removes the mesoscopic and far tail from the
live root-field obligation.

An exact even backward-heat polynomial supplies the decisive scope guard: it
has collision field `-pi/8` and center drift `-pi/4` but still undergoes
positive square-root birth. Accordingly this node sharpens the bridge while
blocking promotion from field balance to no collision. It has no proving edge
to `lambda_le_0_goal`.

The natural renormalized-energy upgrade has now been audited at the boundary:

```text
outputs/jensen_window_pf_newman_boundary_energy_direction_gate.md
work/rh_compute/results/jensen_window_pf_newman_boundary_energy_direction_gate.json
jensen_window_pf_newman_boundary_energy_direction_gate
```

For any double-zero birth, the squared gap begins
`q=8*(t-t_*)-16*K*(t-t_*)^2+...`. The exact classical-field model has
`K>0`, negative quadratic gap curvature, and positive cubic gap jet, while
still producing the pair. Its Rodgers-Tao renormalized interaction is
asymptotic to `1/(8*(t-t_*))`: it decreases forward but is not locally
integrable at the birth time. Finite Xi boundary-integrated energy would
therefore be a valid collision obstruction.

The published theorem does not supply that estimate. It assumes `Lambda<0`
and integrates from `Lambda/4` to `0`, strictly after the boundary, so it is an
interior relaxation theorem in the opposite assumption direction. This node
sharpens the open bridge to an Xi-specific boundary-trace theorem and has no
proving edge to `lambda_le_0_goal`.

Positive-time asymptotics now close the possible height-escape loophole:

```text
outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md
work/rh_compute/results/jensen_window_pf_newman_positive_boundary_attainment_lemma.json
jensen_window_pf_newman_positive_boundary_attainment_lemma
```

If `Lambda>0`, nonreal zeros at times increasing to `Lambda` are trapped by
the absolute strip theorem and the uniform `exp(C/t)` high-zero theorem in the
fixed rectangle `|Re z|<exp(2C/Lambda)`, `|Im z|<=1`. Local uniform convergence
and conjugation symmetry force a finite real multiple zero of `H_Lambda`.
Thus `Lambda<=0` is exactly equivalent to simplicity of every `H_t`, `t>0`.

For a multiplicity-`m` boundary zero, the universal Hermite split gives ordered
cluster energy `m(m-1)/(8*(t-Lambda))+O((t-Lambda)^(-1/2))`. Every possible
multiplicity therefore has a nonintegrable endpoint trace. This node sharpens
the bridge to Xi positive-time simplicity or finite-truncation endpoint-energy
integrability; it does not prove either condition and has no proving edge to
`lambda_le_0_goal`.

The simplicity target now has a first-correlation/Wiener formulation:

```text
outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md
work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_correlation_target.json
jensen_window_pf_newman_strict_laguerre_correlation_target
```

For `phi_t(u)=exp(t*u^2)Phi(u)`, define
`K_(1,t)(v)=integral phi_t(s+v)phi_t(s-v)s^2 ds`. Its Fourier transform is
exactly `H_t'(xi/2)^2-H_t(xi/2)H_t''(xi/2)`. Positive-boundary attainment and
Wiener's theorem, together with the published Platt-Trudgian bound
`Lambda<=1/5`, show that `Lambda<=0` is equivalent to density in `L1(R)` of
the translates of every `K_(1,t)`, `0<t<=1/5`.

An exact triangular-Gaussian convolution is smooth, positive, even, strictly
log-concave, and positive definite, yet its nonnegative Fourier transform has
double zeros. The graph therefore blocks promotion from those generic kernel
properties and records Xi-specific translate density, strict Fourier
positivity, or a stronger total-positive factorization as the open theorem.
This target has only a `would_imply_if_proved` edge to `lambda_le_0_goal`.

A complementary Polymath-15 high-frequency route now has one proved
asymptotic node, two exact route-classification nodes, and one sharply
delimited open node:

```text
outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md
outputs/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.md
outputs/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.md
outputs/jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.md
outputs/jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.md
outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.md
jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem
jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate
jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate
jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit
jensen_window_pf_newman_polymath15_lambda01965_provenance_audit
jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target
```

For
`c_*=4911678521/1933561194=2.540223984760008...`, the exponent-pair
envelope, zeta-jet handoff, and exact Polymath-15 remainder transfer prove
first-Laguerre positivity on every fixed scaled ray `c>=c_*+epsilon` once
`L` is sufficiently large. This is an asymptotic theorem, but it supplies no
practical `L_epsilon` and has no edge to `lambda_le_0_goal`.

The wall gate then separates two analytically different burdens. On
`2<c<=c_*`, the Euler product already gives a nonzero zeta floor and the
missing input is stronger weighted exponential-sum cancellation. At `c=2`,
a strict radius-proportional cancellation gain would hand off conditionally
to published zeta 1-line bounds. For fixed `c<2`, the saddle lies on
`Re(s_*)<1`, so exponent-pair improvement alone does not supply the missing
phase-amplitude floor. This is a route boundary, not an impossibility theorem
and not a proof of the inner Wronskian target.

The Gaussian/Legendre gate then audits the natural heat-semigroup rewrite
without adding a new analytic hypothesis. The finite Gaussian-shift identity
is exact, and Legendre optimization shows that inserting the same pointwise
partial-sum exponent profile reproduces the weighted dyadic frontier exactly.
At `c_*`, equality occurs at
`q_*=4800718975/7734244776`; at `c=2`, the unchanged deficit is
`3133668399/48144906818`. Thus a same-input repackaging cannot lower `c_*`:
the outer route needs a genuinely better partial-sum profile or cancellation
beyond pointwise majorants. This node has no proving edge to the conclusion.

The ANTEDB beta-frontier audit then checks whether that contact was an artifact
of retaining the 2023 exponent-pair table. At pinned commit
`99668603896af86e6cda90ed6755cf3116aab0ac`, four Tao-Trudgian-Yang and two
Cushing post-2023 pairs improve intermediate radii, but both the exact current
rational hull and the documented one-pass direct-beta envelope retain
`alpha_*=62831/155153`, `beta_*=220633/620612`, and
`c_*=4911678521/1933561194`. Twelve finite beta-only van der Corput passes also
preserve that contact. This is a current-source finite nonpromotion result, not
a claim of infinite transform closure or that future estimates cannot improve
the frontier. It has no proving edge to the conclusion.

The Lambda 0.1965 provenance audit keeps a separate literature-status gate.
The peer-reviewed record is `0<=Lambda<=1/5`; a June 2026 Zenodo preprint
claims `Lambda<=393/2000`. Its archive integrity, exact row arithmetic, and
22 portable certificate invocations passed locally, while three FLINT/Arb
producer packages and independent theorem-to-code review remain outside this
audit. The graph therefore labels 0.1965 as a reproduced unrefereed candidate.
If accepted, its only edge sharpens the hypothetical positive-boundary range
to `(0,393/2000]`; it supplies no exclusion theorem and no proving edge.

The surviving high-frequency target is the corrected finite Riemann-Siegel
coercivity inequality on `0<=c<=c_*+o(1)`. Its edge supports the broader
strict-correlation target only: even proving it would still leave the compact
`x` and all-time uniform arguments visible rather than promoting directly to
the conclusion.

The first stronger correlation criterion has now been audited:

```text
outputs/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.md
work/rh_compute/results/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.json
jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate
```

All correlations satisfy `partial_t K_n=2v^2 K_n+2K_(n+1)`, and their
Fourier transforms are the generalized Laguerre hierarchy. At a hypothetical
double boundary root, the universal contact is
`F_1=F_1'=0`, `F_2=3F_1''>0`. Complete monotonicity in `v^2` would give a
positive Gaussian mixture and strict Fourier positivity, but the exact bound
`K_(1,t)(v)<=C exp(2tv^2+18v-pi exp(4v))` makes this impossible. The same
two-sided super-Gaussian tail rules out direct PF-infinity membership, while
smooth evenness rules out the classical decreasing-convex Polya criterion.

The graph therefore sharpens the strict-correlation target to
tail-compatible spectral-square or relative hierarchy-coercivity arguments.
This exact gate blocks three attractive generic promotions but has no edge to
`lambda_le_0_goal`.

The surviving shape information is now exact and uniform on the target window:

```text
outputs/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.md
work/rh_compute/results/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.json
jensen_window_pf_newman_positive_time_strong_logconcavity_gate
```

The published strict concavity of `log(Phi(sqrt(r)))` implies
`(log Phi)''(u)<=-kappa` globally, where an Arb origin certificate gives
`kappa=74.9076197180...`. Consequently `phi_t=exp(tu^2)Phi` is uniformly
strongly log-concave and admissible throughout `0<=t<=1/5`. The joint Hessian
bound and Prekopa marginalization propagate
`(log K_(n,t))''<=-2(kappa-2t)` to every correlation.

This strong shape theorem still cannot be promoted to transform zero-freeness:
`Fourier[K_(0,t)](xi)=2H_t(xi/2)^2`, and Hardy's theorem gives real double
zeros at `t=0`. Thus even an Xi correlation with the certified strong
log-concavity, positive definiteness, and super-Gaussian tail can have Fourier
zeros. The graph sharpens the live target to the `s^2` weighting or a relative
hierarchy-coercivity argument and gives this node no edge to
`lambda_le_0_goal`.

The generic first-weighted-correlation implication is now closed as well:

```text
outputs/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.md
work/rh_compute/results/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json
jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate
```

For
`phi_(delta,epsilon)(u)=exp(-u^2-delta*u^4-epsilon(cosh(4u)-1))*(1+u^4/10)`,
an elementary rational bound gives
`(log phi_(delta,epsilon))''<=-(2-6/sqrt(10))-12delta*u^2-16epsilon*cosh(4u)`.
At `delta=1/10`, `epsilon=1/1000`, this is a smooth positive even uniformly
strongly log-concave admissible kernel with theta-type double-exponential decay.
It also satisfies

```text
d^2/dr^2[-r-r^2/10+log(1+r^2/10)]
 =-r^2(r^2+30)/(5(r^2+10)^2)<0,
d^2/dr^2 cosh(4sqrt(r))>0               (r>0),
```

so subtracting `(cosh(4sqrt(r))-1)/1000` preserves strict root-variable
concavity, matching the stronger shape known for the Xi kernel. At
`delta=epsilon=0`, exact Gaussian differentiation supplies the cross-check

```text
L_1[F_0](3)=-743*pi*exp(-9/2)/51200<0.
```

A 256-bit Acb integral with explicit Arb tails directly certifies

```text
L_1[F_(1/10,1/1000)](21/5)
 =[-0.000213478480591774964507646766853475076598648832979809
   +/- 5.09e-55]<0.
```

The exact correlation identity yields
`Fourier[K_(1,1/10,1/1000)](42/5)<0`. Thus even root-variable concavity,
uniform strong log-concavity, theta-tail decay, and admissibility fail after
the `s^2` weighting: they do not imply positive definiteness of `K_1`, let
alone strict positivity. This countermodel blocks generic
Xi-style-shape-and-tail-plus-weighting promotion and leaves Xi arithmetic,
modular grouping, or relative hierarchy coercivity as the live input. It has
no edge to `lambda_le_0_goal`.

The stronger strict-monotonicity subroute now has a rigorous counter-gate:

```text
outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md
work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.json
jensen_window_pf_newman_strict_laguerre_monotonicity_scout
```

Writing `L_t=H_t'^2-H_tH_t''` and
`M_t=-L_t'=H_tH_t'''-H_t'H_t''`, the exact zero limit gives

```text
M_t(x)>0 for x>0
  => L_t(x)=integral_x^infinity M_t(y)dy>0.
```

Thus this strict monotonicity throughout `0<t<=1/5` would imply
`Lambda<=0`. The correlation and endpoint-subtracted theta forms are

```text
M_t(x)=4 integral_0^infinity v*K_(1,t)(v)*sin(2xv)dv,
64M_t=(A_t+1/2)A_t'''-A_t'A_t'',  A_t=8H_t-1/2.
```

For `Q_t(s)=integral_s^infinity v*K_(1,t)(v)dv`, one also has

```text
M_t(x)=8x integral_0^infinity Q_t(s)cos(2xs)ds.
```

Strong log-concavity makes `Q_t` concave then convex with exactly one
inflection, and parity blocks the elementary fixed-interior shift criterion.
Those obstructions no longer call for a harder boundary contour theorem,
because the desired sign itself is false. At
`x=1401016343/100000`, a 256-bit Arb Xi jet certifies

```text
M_0(x)/A_0(x)^2=-3.2418674697568...e-6<0,
M_0(x)/L_0(x)=-0.0001463088540165...<0.
```

Continuity in heat time preserves this sign for all sufficiently small
positive `t`, rejecting the universal condition. The exact close-pair law

```text
M[F](m)=-4a^2G(m)G'(m)+a^4(G(m)G'''(m)-G'(m)G''(m))
```

explains why all 12000 moderate rows and 20 selected rows through `x=200`
had the expected sign. The graph now gives this node a rejection edge to the
open strict-correlation target: it closes only global strict decrease, while
direct strict-Laguerre positivity, Wiener density, and corrected `C1`
double-zero transversality remain live.

The surviving arithmetic spectral-square route has also been audited:

```text
outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md
work/rh_compute/results/jensen_window_pf_newman_theta_summand_spectral_square_gate.json
jensen_window_pf_newman_theta_summand_spectral_square_gate
```

The exact theta primitive satisfies `Phi=(R''-R)/8` and `R'(0)=-1/2`. For
`C_t=Fourier_cos[exp(tu^2)R]`, it gives
`H_t=1/16+D_t[C_t]/8` with
`D_t=-4t^2 d_x^2+4tx d_x+(2t-1-x^2)`. Thus strict Laguerre positivity is one
explicit curvature inequality for `A_t=D_t[C_t]`. The bilateral
shifted-profile transform times the Dirichlet shift sum is exactly
`xi((1+iz)/2)/4`, so direct Mellin assembly reconstructs the object whose
zeros are at issue rather than a zero-free factor.

More sharply, every finite theta truncation retains
`A_N=f_(N,t)'(0)>0`. Its first Laguerre expression is therefore
`-2A_N^2/x^6+O(x^-8)<0` at sufficiently high frequency. Pairwise self/cross
forms also have explicit negative numerical witnesses. The graph closes all
finite summand-square promotions and sharpens the live route to an
noncircular curvature estimate after exact infinite modular cancellation.
This node has no edge to `lambda_le_0_goal`.

The Gasper fake-Xi benchmark has now been tested as a comparison route:

```text
outputs/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.md
work/rh_compute/results/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.json
jensen_window_pf_newman_gasper_fake_xi_remainder_gate
```

The scaled fake-Xi transform is an established real-zero model, and the exact
identity `L[H]=alpha^2 L[P]+B[alpha P,E_alpha]+L[E_alpha]` isolates its margin.
At `x=25`, however, exact scalar minimization plus Arb midpoint certificates
with explicit tails puts `inf R_alpha` strictly above one at both `t=0` and
`t=1/2`; independent 80-digit values are `2.5852...` and `2.5692...`.

The stronger direct convolution route fails exactly: `M=Phi/Psi` satisfies
`M(0)>1` but tends to one at infinity. A positive cosh transform with a finite
tail limit must be constant, so no Cardon/Polya one-factor positive multiplier
can transfer the fake-Xi zero theorem to Xi. The graph closes these two
promotions and leaves only sign-aware or multi-block Gasper structure feeding
the open strict-correlation target. This node has no edge to
`lambda_le_0_goal`.

The positive-residual two-block extension has now been closed:

```text
outputs/jensen_window_pf_newman_gasper_residual_two_block_gate.md
work/rh_compute/results/jensen_window_pf_newman_gasper_residual_two_block_gate.json
jensen_window_pf_newman_gasper_residual_two_block_gate
```

The first residual kernel `Phi-Psi_9` is exactly positive. For the full
nonnegative 9/5 family, however, tail positivity forces
`0<=beta<=pi-3/(2pi)`, while the residual Laguerre expression is one exact
convex quadratic in beta. Acb Cauchy derivative enclosures and Arb endpoint
signs prove it negative at `x=66` on `0<=beta<=23/10` and at `x=50` on the
remaining interval. Larger beta makes the residual kernel eventually
negative.

The tail-matched multiplier also has a negative quartic discriminant and
off-imaginary zeros, so the standard universal-factor proof is unavailable.
The graph therefore closes a positive 9/5 Gasper block plus independently LP
positive residual. Only signed blocks or a coupled mixed-term square remain;
this node has no edge to `lambda_le_0_goal`.

The classical three-block extension now closes the full positive 9/5/1
architecture:

```text
outputs/jensen_window_pf_newman_classical_three_block_residual_gate.md
work/rh_compute/results/jensen_window_pf_newman_classical_three_block_residual_gate.json
jensen_window_pf_newman_classical_three_block_residual_gate
```

The established Polya P2 and de Bruijn real-zero kernels both leave strictly
positive Xi-kernel residuals. The residual transforms nevertheless satisfy
`L[R_P2](86)<0` and `L[R_dB](52)<0`. More generally, exact tail and origin
constraints put every nonnegative 9/5/1 block with a globally nonnegative
residual into one rational parameter triangle. Three Acb spectral
certificates cover all `64,908` closed boxes in that triangle and give a
negative first Laguerre enclosure on every box.

Gasper's published one-shift modulus square does not sign the mixed products
created by a multi-shift linear combination. The graph therefore leaves only
signed higher blocks or a genuinely coupled square after modular
cancellation. This node has no edge to `lambda_le_0_goal`.

The full signed three-shift universal-factor extension is also closed:

```text
outputs/jensen_window_pf_newman_signed_universal_factor_residual_gate.md
work/rh_compute/results/jensen_window_pf_newman_signed_universal_factor_residual_gate.json
jensen_window_pf_newman_signed_universal_factor_residual_gate
```

The signed multiplier `cosh(9z)+beta cosh(5z)+gamma cosh(z)` has only
imaginary zeros exactly when one explicit quartic has all roots in `[0,1]`.
Endpoint conditions, the critical cubic, and the Xi origin bound put every
candidate with a globally nonnegative residual in one rational rectangle.
Two Acb residual Laguerre tests and exact critical and quartic discriminants
classify `4,094` adaptive leaves grown from `3,416` base boxes. The maximum
depth is six and no leaf is unresolved.

Thus the standard signed 9/5/1 universal-factor cone cannot be paired with an
independently LP positive residual. The graph now hands only higher shifts or
a genuinely coupled matrix square to the open strict-correlation target. This
node has no edge to `lambda_le_0_goal`.

The higher-shift option now has an exact regularization gate:

```text
outputs/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.md
work/rh_compute/results/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.json
jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate
```

Every symmetrized arithmetic theta block has an exact Bessel-K expansion with
coefficients

```text
c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3m)/m!.
```

For fixed `n` the expansion has an integrable absolute majorant and may be
transformed term by term. The coefficients are intrinsically signed, with the
first `n=1` negative coefficient at `m=7`. More decisively, the fixed-block
zero-frequency transforms equal a negative constant times `n^(-1/2)`, so
their ordinary arithmetic sum diverges while the complete Xi transform is
finite. The graph therefore closes naive termwise higher-shift summation and
hands only modularly grouped mixed-term matrices or a sign-preserving
renormalized summation theorem to the strict-correlation target. This node has
no edge to `lambda_le_0_goal`.

The divergent higher-shift transform now has an exact endpoint
renormalization:

```text
outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md
work/rh_compute/results/jensen_window_pf_newman_theta_cell_renormalization_gate.json
jensen_window_pf_newman_theta_cell_renormalization_gate
```

The continuous theta-index integral vanishes pointwise. Subtracting its cells
from the integer blocks therefore leaves `Phi` unchanged while producing
normal convergence in every polynomially weighted `L1` norm. The transformed
series is Euler's convergent sum-integral formula for zeta, every renormalized
block is positive at zero frequency, and the endpoint first Laguerre
expression becomes an absolutely convergent mixed-term matrix.

This repair is endpoint-only. Every cell block has tail

```text
r_n(u)=-(pi/2)*(3n-1)exp(-5u)+O_n(exp(-9u)),
```

so multiplication by `exp(tu^2)` is nonintegrable for every `t>0`. The graph
therefore hands a t-compatible modular grouping, not the endpoint cell blocks,
to the strict-correlation target. This node has no edge to
`lambda_le_0_goal`.

The moment integral also has an exact Laguerre scale-kernel gate:

```text
outputs/jensen_window_pf_laguerre_scale_mixture_gate.md
work/rh_compute/results/jensen_window_pf_laguerre_scale_mixture_gate.json
jensen_window_pf_laguerre_scale_mixture_gate
```

Every fixed-scale kernel is a rescaled `L_D^(-1/2)` and has simple negative
roots. Exact two-atom and exponential-density countermodels show that positive
mixing and log-concavity do not preserve this property. Half-integer Gamma
mixing is hyperbolic in every degree by an Euler-Jacobi factorization. The
graph therefore records an Xi-specific common-interlacing or total-positive
connection theorem as the live handoff, not generic mixture positivity.

The cubic multiplicity boundary extends to an exact all-degree benchmark:

```text
outputs/jensen_window_pf_rank_two_boundary_family_lemma.md
work/rh_compute/results/jensen_window_pf_rank_two_boundary_family_lemma.json
jensen_window_pf_rank_two_boundary_family_lemma
```

For `1/2<u<1`, the model coefficients `A_k=a^(k-1)(c+k*b)/u` factor every
shifted Jensen window into one simple and one repeated negative root. This is
an exact all-degree multiplier family, and finite pointwise products remain
multiplier sequences. In alpha coordinates they give
`-log x_k=sum_j -log(1-1/(k+alpha_j)^2)`. However, the equal positive mixture
of the `u=3/5` and `u=2/3` sequences has degree-3 discriminant `-937/3456`,
and a fractional pointwise power has discriminant `<-27/125`. Thus the live
canonical-product target requires a genuine counting measure; naive convex or
arbitrary positive-measure integration is blocked.

That sufficient subclass was isolated as an explicit target:

```text
outputs/jensen_window_pf_multiplier_counting_measure_target.md
work/rh_compute/results/jensen_window_pf_multiplier_counting_measure_target.json
jensen_window_pf_multiplier_counting_measure_target
```

After normalizing `B_0=B_1=1`, it asks for a convergent unit-atomic product
`B_k=product_j M_k^(alpha_j)`, equivalently
`-log x_k=sum_j -log(1-1/(k+alpha_j)^2)`, derived directly from Phi/Newman
data. The later interval certificate now rejects this sufficient subclass:

```text
outputs/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.md
work/rh_compute/results/jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.json
jensen_window_pf_multiplier_unit_atomic_obstruction_certificate
```

The order-6 magnitude bound forces every unit atom above `4.863538496`.
Exact monotonicity of the one-atom consecutive-moment ratio would then cap
`a_7/a_6`, but Arb proves the actual ratio exceeds that cap by more than
`7.68e-4`. The target node is therefore a validated rejected-route record,
not an open blocker. General multiplier sequences and the other all-order
Jensen/PF bridges remain open.

The natural continuous-index strengthening now has a rigorous obstruction:

```text
outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md
work/rh_compute/results/jensen_window_pf_mellin_multiplier_power_sum_obstruction.json
jensen_window_pf_mellin_multiplier_power_sum_obstruction
```

Eleven Arb-enclosed `Phi` log moments induce the power sums required by a
continuous elementary-multiplier product. Three shifted Hankel determinants
are strictly negative, including the shift-2 size-4 determinant near
`-2.588644974358276e-19`. This blocks promotion from the integer target to the
natural Mellin interpolation. It does not reject a product identity proved
directly at integer indices, so the counting-measure target remains open.

Finally, the full-real-T fixed-k certificate combines the bounded segment
with a full-kernel ray proof:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.json
```

For `T>=10000`, fix `x_*=1/5`, so the real tail begins at `y>=400`.
Full-kernel geometric n-tail bounds, the order-`42` compact factor, and an
upper-Gamma hazard inequality prove the scaled ray bounds decrease with T;
the worst budget fractions are below `1.1e-15` and `1.3e-16`. Together with
the bounded segment, the fixed-`k=22` normalizers, B product, companion
product, and weighted-gap derivative are certified for every real
`T>=1156`. The T=1156 counterexample proves that this fixed-k result cannot be
extended to all k at its left endpoint; the theorem itself remains valid.
This interval-certificate node has no edge to `lambda_le_0_goal`.

The first-omitted denominator certificate turns part of that target into
interval-ready arithmetic:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json
```

It Arb-certifies `r_21` negative and bounded away from zero, then proves every
recorded first-omitted value and derivative denominator positive. The weakest
denominators occur at `T=10000`, `F_21`, giving scaled absolute-radius targets
about `6.78e-28` and `1.42e-30` for a `1e-6` ratio-radius cap. The graph
records this as denominator-side support only; residual numerator enclosure,
quadrature error, rounding, and grid-to-collar coverage remain open.

The coefficient-core propagation certificate discharges the coefficient-ball
source for the recorded finite grid:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json
```

It rebuilds Arb coefficient-ratio balls `r_0` through `r_21` and propagates
their radii through exact Gamma moments in first-omitted units. The worst value
and derivative coefficient-radius ratios are about `8.61e-81` and `5.52e-83`
at `T=10000`, `F_21`, far below the proposed `1e-6` ratio-radius target. The
graph records this as coefficient-radius support only; Phi node evaluation,
Laguerre node/weight intervals, quadrature error, rounding aggregation, and
grid-to-collar coverage remain open.

The worst-row Laguerre root-bracket certificate attacks the node half of the
same intervalization source:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json
```

It Arb-certifies 320 disjoint sign-changing brackets for the roots of
`L_320^(41/2)` at `T=10000`, `F_21`, with widest bracket about `4.15e-17`.
It also records that SciPy double weights underflow to zero for 30 tail nodes
on that row. The graph records this as worst-row node support only; non-floating
Christoffel-weight intervals, Phi node evaluation, quadrature remainder,
rounding aggregation, all-row coverage, and grid-to-collar coverage remain
open.

The worst-row Christoffel-weight midpoint scout repairs that floating
underflow diagnostically:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.json
```

It evaluates the generalized Laguerre Christoffel formula at Arb midpoints of
the certified root brackets, giving 320 positive non-floating midpoint weights.
The scout repairs all 30 SciPy double underflows and matches `Gamma(43/2)` to
relative error about `1.8e-18`. The graph records this as midpoint evidence
only; direct interval denominator evaluation still contains zero on all 320
rows, which is why the follow-up certificate uses a Taylor-centered
denominator enclosure instead of promoting the midpoint table.

The worst-row Christoffel-weight interval certificate closes that one weight
subproblem:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json
```

It evaluates `L_321^(41/2)` on each certified `L_320^(41/2)` root bracket by
the exact Taylor identity
`d^m L_n^(alpha)/dx^m=(-1)^m L_(n-m)^(alpha+m)`. All 320 denominator intervals
avoid zero, all 320 Christoffel-weight intervals are positive, all 30 SciPy
double underflows are repaired as interval weights, and the independent
weight-sum interval contains `Gamma(43/2)`. The graph still records this as
worst-row weight evidence only: Phi/Phi' node evaluation, quadrature remainder,
rounding aggregation, all-row coverage, and grid-to-collar coverage remain
open.

The worst-row finite-part weighted-sum interval certificate then closes the
finite `n<=30` weighted sum for the same row:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json
```

It refines the `T=10000`, `F_21`, `N=320` node brackets to 120 bisection
steps, evaluates the finite Phi/Phi' terms on those intervals, sums with
certified Christoffel weights, and subtracts the polynomial part through exact
Gamma moments. The scaled value and derivative residual ratios are certified
below one first omitted term. The graph records this as finite-part
worst-row evidence only; n-tail composition, quadrature remainder, rounding,
all-row coverage, and grid-to-collar coverage remain open.

The worst-row finite-plus-tail budget certificate composes that finite-part
certificate with the reserved Phi-tail source cap:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json
```

It adds the full `2.0e-3` Phi-tail source budget to the certified finite
`n<=30` ratio uppers. The composed value and derivative ratio uppers are
`0.9853957992836557769419015895036210773888` and
`0.9714055674762067320093698741711260875260`, still below one first omitted
term. The graph records this as worst-row finite-plus-tail budget evidence
only; quadrature remainder, rounding aggregation, all-row coverage, and
grid-to-collar coverage remain open, and it has no edge to
`lambda_le_0_goal`.

The raw-moment obstruction matrix blocks the generic Stieltjes shortcut:

```text
outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.json
```

It gives three exact positive two-atom moment-sequence witnesses: one violates
the upper raw wall, one violates the scaled-upper corridor side while the raw
upper wall holds, and one violates the monotone-bridge lower corridor side
while the raw upper wall holds. This forces the adaptive route to use
zeta-specific all-`k` structure and has only `blocks_promotion` edges, never a
proving edge to `lambda_le_0_goal`.

The historical zeta-specific raw-corridor target names that structure:

```text
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.json
```

The target is now discharged at `lambda=-100` by

```text
outputs/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.json
```

Full ratio-cone entry supplies `B_k>=0` and `B_(k+1)<=B_k`; the scaled-
curvature theorem supplies `B_(k+1)>=((2*k+1)/(2*k+3))*B_k`; and the exact
linear-to-nonlinear calculus lemma proves the complete raw-ratio decrement
corridor for every `k>=1`. This node is `interval_validated` and is no longer
a blocker for `lambda_le_0_goal`; the all-order Jensen/PF bridge remains open.

The adaptive-defect composition carries that theorem into the old defect
coordinates:

```text
outputs/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.json
```

At lambda=-100, full cone entry gives
`0<=d_k<=2/(2*k+1)`, `d_(k+1)<=d_k`, and `0<=s_k<=1` for every `k>=1`.
The upper raw-corridor wall is exactly `s_(k+1)>=s_k`. This closes the
parameter-specific defect-tail and adaptive-defect routes while preserving
the stronger simultaneous lambda=-25,-50,-100 statements as unproved audit
history. The certificate has no edge to `lambda_le_0_goal`.

The monotone-contraction stress artifact supports that frontier scout on the
existing zeta coefficient cache:

```text
outputs/jensen_window_pf_monotone_contraction_stress.md
work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json
```

It validates 2875 finite Arb-classified rows across degrees d=3..12 and the
k<=64 cache, all satisfying adjacent log-concavity and increasing ratio
contractions. It is finite evidence only; it does not prove an analytic
monotone-contraction theorem, and it has no edge to `lambda_le_0_goal`.

The signed J-fraction theorem target then states the exact missing theorem and
acceptance gates needed to turn that finite pattern into an all-order column
recurrence theorem. It validates 7 fit/misfit rows and 0 ready-to-apply rows.
It is a hygiene node only and has no edge to `lambda_le_0_goal`.

The modified signed-model target sharpens that theorem target by classifying
dead and live model classes after the Motzkin and parity-lift obstructions:

```text
outputs/jensen_window_pf_modified_signed_model_target.md
work/rh_compute/results/jensen_window_pf_modified_signed_model_target.json
```

It rejects the raw ordinary Motzkin/J-fraction model, diagonal sign
conjugation, and global length-parity sign repairs, while leaving state-space
doubled models, modified production/Riordan/network models, oscillatory
sign-regular resolvent theorems, and Xi/Phi positive cancellation identities
live only conditionally. It validates 9 model rows and 0 ready-to-apply rows.
It is a hygiene node only and has no edge to `lambda_le_0_goal`.

The Cauchy-Binet low-degree scout sharpens the structural ansatz row:

```text
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
```

It records that the selected degree-2/3/4 hard formulas are certified by
adjacent log-concavity with nonnegative Bernstein coefficients, while `0`
kernel identities are found. It is a diagnostic node only and has no edge to
`lambda_le_0_goal`.

The log-concavity frontier scout extends that diagnostic:

```text
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
```

It records first Bernstein-certificate failures at degree 3 size 6 and degree
4 size 5, and first exact countermodel negatives at degree 3 size 8 and degree
4 size 6. It is also a diagnostic node only and has no edge to
`lambda_le_0_goal`.

The ratio-condition scout audits nearby strengthened assumptions:

```text
outputs/jensen_window_pf_ratio_condition_scout.md
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
outputs/jensen_window_pf_contraction_log_concavity_scout.md
work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json
```

It rejects adjacent log-concavity, decreasing ratio contractions,
second-order log-concavity, and selected low-degree Bernstein positivity by
exact countermodel. The contraction-log-concavity scout sharpens it by
rejecting the last standalone ratio candidate through a constructed positive
extension. Both are diagnostic only and have no edge to `lambda_le_0_goal`.

## Non-Promotion Region

```text
jensen_hankel_bridge_algebra_gate
jensen_window_pf_obligation_algebra_gate
countermodel_gates
rejected_finite_prefix_promotion
```

The sign-regular transfer gap matrix and factorial multiplier split audit sit
in the obligation-decomposition region and draw only `sharpens` edges from the
exact algebra, countermodel, and threshold gates into the theorem-machinery
audit. These nodes only sharpen the theorem targets or validate rejection of
finite prefix promotion. They do not supply a theorem proving path to the
conclusion.

## Validation Rules

The checker enforces:

```text
all ledger-backed graph nodes exist in the proof-claim ledger
ledger-backed graph roles and statuses match the ledger categories and statuses
all edge endpoints exist
all edge types are in the allowed set
no forbidden edge type appears
lambda_le_0_goal remains not_proved
would_imply_if_proved edges are sourced only from open theorem targets
finite and diagnostic nodes do not point directly to lambda_le_0_goal
the Jensen-window evidence nodes feed only into open bridge targets
the countermodel nodes block finite-prefix promotion
```

## Boundary

Passing this checker means the dependency ledger has not silently promoted a
finite signed-Hankel, finite Jensen-window, Sturm, PF, countermodel, or hygiene
artifact into the missing Newman-direction theorem. It is a proof-programme
firewall, not a proof of the firewall's target.
