# Jensen-Window PF Cauchy-Binet Cone Frontier Matrix

Date: 2026-07-06

Status: Cauchy-Binet cone frontier matrix. This is not a proof of a positive
kernel, Cauchy-Binet identity, determinant integral, Jensen-window
PF-infinity, Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The nonpower functional cone matrix leaves a live row:

```text
npc_07_cauchy_binet_determinant_integral_cone
```

The older Cauchy-Binet low-degree scout found:

```text
15 selected low-degree formulas
0 kernel identities found
```

This matrix asks what a serious Cauchy-Binet cone must do beyond those
selected formulas.

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_cone_frontier_matrix.py
```

Current result:

```text
validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows
```

## Frontier Contract

A viable Cauchy-Binet cone must:

```text
1. produce an exact Cauchy-Binet, Andreief, Gram, or determinant-integral
   identity for the actual determinant or reciprocal coefficient;
2. prove every factor, measure, kernel, or minor in the identity is
   nonnegative from Xi/Phi or zeta heat-flow data;
3. cover d3_column_recurrence_m8 and d4_column_recurrence_m6, not only the
   selected low-degree minors;
4. state whether it is column-only or all-shape;
5. avoid endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and
   Lambda <= 0;
6. prove all-degree, all-shift, and all-m coverage;
7. explain why the exact rational log-concave countermodel is outside the new
   zeta-specific hypotheses without smuggling in the target.
```

## Frontier Rows

```text
cbcf_01_selected_low_degree_bernstein_certificate:
  low_degree_certificate_only
  The 15 selected low-degree Bernstein certificates are not a kernel identity.

cbcf_02_adjacent_log_concavity_cone:
  rejected_as_standalone_bridge
  Adjacent log-concavity and nearby ratio conditions are rejected as
  standalone bridge cones.

cbcf_03_column_frontier_determinant_integral:
  live_if_frontier_identity_constructed
  A positive determinant integral for d3_column_recurrence_m8 and
  d4_column_recurrence_m6 remains live. The monotone-contraction frontier
  scout gives an exact sufficient condition for these first two hard
  polynomials, the monotone-contraction theorem target isolates the missing
  analytic inequality, and the monotone-contraction stress artifact supports
  that condition on a wider finite zeta grid. The all-order analytic zeta
  monotone-contraction theorem remains open.

cbcf_04_all_shape_cauchy_binet_kernel:
  live_if_all_shape_kernel_constructed
  A uniform Cauchy-Binet or Andreief kernel for every Toeplitz/Jacobi-Trudi
  shape remains live.

cbcf_05_gram_or_hankel_square_kernel:
  conditional_wrapper_only
  Gram language is only a wrapper unless it avoids the raw moment obstruction
  Delta_2(mu)=-g_2 and supplies the exact nonordinary kernel.

cbcf_06_finite_quadrature_or_minor_fit:
  finite_evidence_only
  Finite positive minor fits can suggest a kernel but cannot prove all
  degree, shift, and shape cases.

cbcf_07_endpoint_factorization_integral:
  circular_if_endpoint_used
  Endpoint PF, real-rootedness, or Laguerre-Polya factorization is circular as
  an input.

cbcf_08_indefinite_andreief_or_signed_kernel:
  language_only
  Signed or indefinite kernel language does not prove final nonnegativity.
```

## Evidence Anchors

```text
outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
outputs/jensen_window_pf_log_concavity_frontier_scout.md
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_contraction_log_concavity_scout.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
outputs/jensen_window_pf_monotone_contraction_frontier_scout.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_monotone_contraction_stress.md
outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/jensen_window_pf_schur_shape_contract.md
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
outputs/jensen_window_pf_positive_readout_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/arb_jensen_window_column_recurrence_stress.md
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md
outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md
```

New column-frontier sharpening:

```text
validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues
validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues
```

## Kill Gates

Reject a proposed Cauchy-Binet cone if it:

```text
1. treats selected low-degree Bernstein certificates as a kernel identity;
2. uses adjacent log-concavity or ratio conditions as standalone membership;
3. promotes finite quadrature, finite grids, or fitted minor decompositions to
   an all-order theorem;
4. assumes endpoint PF, real-rootedness, Laguerre-Polya membership, RH, or
   Lambda <= 0;
5. uses raw ordinary Gram/Hankel moment PSD language blocked by
   Delta_2(mu)=-g_2;
6. gives only signed or indefinite Andreief language without a positive final
   identity.
```

## Boundary

Passing this checker means the Cauchy-Binet cone route has explicit frontier
obligations, rejection gates, and two live theorem targets. It does not
construct a positive kernel, determinant integral, reciprocal coefficient
proof, all-shape Cauchy-Binet identity, or proof of `Lambda <= 0`.
