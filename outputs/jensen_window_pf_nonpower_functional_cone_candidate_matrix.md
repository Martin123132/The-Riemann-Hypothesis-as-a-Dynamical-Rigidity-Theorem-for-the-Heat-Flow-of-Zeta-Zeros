# Jensen-Window PF Nonpower Functional Cone Candidate Matrix

Date: 2026-07-06

Status: nonpower functional cone-candidate matrix. This is not a proof of the
all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The low-degree nonpower-functional scout reduces the live `npt_04` route to a
positive-cone contract:

```text
K_m in C
L positive on C
L(K_m) = mu_m = [t^m]1/H(-t)
```

This matrix classifies candidate cones `C`.

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_nonpower_functional_cone_candidate_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_cone_candidate_matrix.py
```

Current result:

```text
validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows
```

## Cone Contract

A viable cone proof must:

```text
1. specify C without using the unknown target coefficients mu_m;
2. prove K_m in C for every m,d,n from Xi/Phi or zeta heat-flow data;
3. prove L is positive on C independently of the desired conclusion;
4. prove L(K_m)=mu_m=[t^m]1/H(-t) for the original reciprocal coefficients;
5. explain the low-degree forcing
   h_1^2-h_0*h_2
   h_1^3-2*h_0*h_1*h_2+h_0^2*h_3;
6. avoid endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH,
   and Lambda <= 0;
7. state all-m, all-degree, and all-shift coverage.
```

## Cone Rows

```text
npc_01_raw_g_coordinate_cone:
  rejected_by_signed_low_degree_polynomials
  Rejected because mu_2=g1^2-g2 and the composition expansion has mixed signs.

npc_02_ratio_log_concavity_cone:
  rejected_as_standalone_bridge
  Rejected as a standalone bridge by exact rational countermodels and
  constructed extensions.

npc_03_finite_zeta_prefix_cone:
  finite_evidence_only
  Finite positive recurrence prefixes are useful clues, not an all-m,d,n cone.

npc_04_tautological_pullback_cone:
  rejected_as_tautological
  A cone defined by the desired mu_m>=0 conclusion is not a proof.

npc_05_endpoint_pf_or_laguerre_polya_cone:
  circular_if_endpoint_used
  Endpoint PF, real-rootedness, or Laguerre-Polya structure is circular as an
  input to the missing bridge.

npc_06_xi_phi_kernel_test_function_cone:
  live_if_kernel_cone_constructed
  A zeta-specific Xi/Phi test-function cone remains live if it gives the exact
  readout.

npc_07_cauchy_binet_determinant_integral_cone:
  live_if_integral_identity_constructed
  A positive minor or determinant-integral cone remains live if it evaluates
  to the reciprocal coefficients exactly and covers the hard column frontier
  or all Toeplitz/Jacobi-Trudi shapes.

npc_08_order_unit_or_sos_cone:
  conditional_wrapper_only
  SOS or semidefinite cone language is only a wrapper until there is a
  noncircular zeta-derived all-m membership theorem.
```

## Evidence Anchors

```text
outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md
outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_contraction_log_concavity_scout.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
outputs/arb_jensen_window_column_recurrence_stress.md
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
outputs/jensen_window_pf_positive_readout_theorem_target.md
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
```

The Cauchy-Binet cone frontier matrix validates:

```text
validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows
```

## Kill Gates

Reject a proposed cone if it:

```text
1. uses only independent g_j nonnegativity;
2. uses ratio/log-concavity as a standalone bridge;
3. promotes finite zeta prefixes to an all-m,d,n cone theorem;
4. defines C or L by the desired mu_m>=0 conclusion;
5. assumes endpoint PF, real-rootedness, Laguerre-Polya membership, RH, or
   Lambda <= 0;
6. gives generic cone language without C, K_m, L, positivity, membership, and
   exact readout proofs.
```

## Boundary

Passing this checker means the `npt_04` cone search has explicit candidate
classes, rejection gates, and live theorem targets. It does not construct the
cone `C`, the basis `K_m`, or the functional `L`; it does not prove reciprocal
coefficient positivity; and it does not prove `Lambda <= 0`.
