# Jensen-Window PF Structural Ansatz Matrix

Date: 2026-07-06

Status: structural proof-search matrix. This is not a proof of Jensen-window
PF-infinity, Jensen hyperbolicity, Laguerre-Polya membership, RH, or
`Lambda <= 0`; it converts the open `jwpf_06` bridge into candidate ansatz
families and exact low-degree tests.

## Purpose

The theorem-machinery audit says no known source theorem is currently
ready to apply to:

```text
jwpf_06_sign_regular_to_jensen_pf_conversion
```

This matrix makes the next proof-search step more concrete. A proposed bridge
must survive the exact degree-2, degree-3, degree-4, and finite countermodel
tests before it is allowed to look like a manuscript lemma.

Machine-readable matrix:

```text
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
```

Current result:

```text
validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows
```

## Hard Low-Degree Tests

The hard tests are imported from the exact algebra gate:

```text
work/rh_compute/results/jensen_window_pf_obligation_algebra.json
python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py
```

They include:

```text
hard_01_degree2_contact:
  4*(a1**2 - a0*a2)
  a1**2 - a0*a2 = -det([[a0,a1],[a1,a2]])

hard_02_degree3_2x2_off_hankel:
  -3*(a0*a2 - 3*a1**2)
  -a0*a3 + 9*a1*a2
  -3*(a1*a3 - 3*a2**2)

hard_03_degree3_3x3_cubic:
  a0**2*a3 - 18*a0*a1*a2 + 27*a1**3
  a0*a3**2 - 18*a1*a2*a3 + 27*a2**3

hard_04_degree4_banded_coefficients:
  -2*(3*a0*a2 - 8*a1**2)
  -4*(a0*a3 - 6*a1*a2)
  -a0*a4 + 16*a1*a3
  -2*(3*a0*a2*a4 - 8*a0*a3**2 - 8*a1**2*a4 + 96*a1*a2*a3 - 108*a2**3)

hard_05_countermodel_kill_gate:
  d3_first_negative_contiguous_toeplitz_minor.size=8
  d4_first_negative_contiguous_toeplitz_minor.size=6
  blocked_promotion includes selected low-order Jensen-window Toeplitz minors
```

These tests are not proof obligations for the zeta coefficients alone; they
are rejection tests for proposed general bridge mechanisms.

## Ansatz Rows

```text
ansatz_01_direct_shifted_hankel_only_transfer:
  rejected_by_countermodel
  matches only the degree-2 contact and fails the degree-3/4 off-Hankel
  Toeplitz obligations

ansatz_02_positive_cauchy_binet_kernel:
  live_structural_candidate
  would close jwpf_06 only if a uniform positive Cauchy-Binet or Andreief
  identity is proved for every d,n and every Toeplitz minor shape; the
  low-degree scout below certifies the selected hard formulas by adjacent
  log-concavity but finds no kernel identity

ansatz_03_planar_network_or_production_matrix:
  live_structural_candidate
  would close jwpf_06 only if a positive network or production matrix produces
  B^{d,n,0}_j=binom(d,j) A_{n+j}(0) and all its Toeplitz minors; the raw
  ordinary Motzkin/J-fraction production matrix is rejected as manifest
  positivity by the Motzkin obstruction scout, and a global length-parity sign
  repair is rejected by same-length mixed-sign path witnesses; an
  absolute-value sign-state cover is rejected at mu_2 by the gap 2*kappa_1

ansatz_04_binomial_multiplier_preserver:
  blocked_until_input_pf
  relevant to the binomial weights, but circular unless an input PF or
  real-rootedness theorem is proved independently

ansatz_05_positive_determinant_integral:
  live_structural_candidate
  would close jwpf_06 only if every Jensen-window Toeplitz minor gets a
  manifestly nonnegative determinant-integral formula

ansatz_06_finite_grid_interpolation:
  rejected_by_countermodel
  finite Arb/Sturm/PF rectangles and interpolation over checked minors cannot
  replace an all-degree/all-shift theorem
```

No row is `ready_to_apply`.

## Schur Shape Contract

The live positive-kernel, planar-network, production-matrix, and
determinant-integral rows now have an exact bounded Schur/Jacobi-Trudi shape
contract:

```text
outputs/jensen_window_pf_schur_shape_contract.md
outputs/jensen_window_pf_column_recurrence_contract.md
work/rh_compute/results/jensen_window_pf_schur_shape_contract.json
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
python work/rh_compute/scripts/jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
```

It validates:

```text
validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows
validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows
```

The contract records `15,709` finite-band nonzero bounded shapes on the
`N=8`, order `<=5`, degree `2..5` grid, and identifies the hard degree 3 size
8 and degree 4 size 6 frontier minors as column-shape Jacobi-Trudi
obligations with mixed-sign h-monomial expansions. The recurrence refinement
records the necessary elementary-symmetric condition `C_m = h_0^m * e_m` and
shows the two hard recurrence rows are negative at the exact rational
countermodel.

The reciprocal-positivity route matrix sharpens this column-recurrence layer:

```text
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
```

It records positive renewal/convolution, positive continued-fraction, and
production-matrix routes as live only if a new Xi/Phi-specific representation
is proved, while keeping generic ratio conditions and finite recurrence stress
out of the proof route.

The ordinary Motzkin production-matrix obstruction is:

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

These reject only the raw ordinary Motzkin/J-fraction path model, global
length-parity sign repairs, and absolute-value sign-state covers. A modified
network, genuinely modified state-space doubled model, or Xi/Phi-specific
production matrix would still be live if it supplies manifest nonnegative
weights and survives the hard frontier checks.

## Cauchy-Binet Low-Degree Scout

The live Cauchy-Binet row now has a symbolic low-degree scout:

```text
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
python work/rh_compute/scripts/jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
```

Current scout result:

```text
validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found
```

Its main finding is that the selected degree-2/3/4 formulas have nonnegative
Bernstein coefficients after adjacent-log-concavity ratio parametrization.
That is a useful low-degree certificate, but it is weaker than the desired
Cauchy-Binet theorem. The same exact countermodel remains log-concave, passes
the selected low-degree minors, and fails larger contiguous Toeplitz minors.

The larger contiguous-minor frontier is:

```text
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
```

It validates:

```text
validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues
```

This frontier shows where adjacent log-concavity stops: first Bernstein
certificate failures at degree 3 size 6 and degree 4 size 5, and first
countermodel negatives at degree 3 size 8 and degree 4 size 6.

The nearby strengthened-ratio audit is:

```text
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_contraction_log_concavity_scout.md
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py
```

It validates:

```text
validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction
validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows
```

## Best Current Direction

The viable search space is structural, not another finite promotion:

```text
positive Cauchy-Binet kernel
positive planar network or production matrix
positive determinant-integral formula
```

Each candidate must reproduce the listed degree-3 and degree-4 determinant
polynomials exactly, while also explaining all shifts, all degrees, and all
Toeplitz minor shapes.

## Integration Points

This workbench sharpens:

```text
outputs/jensen_window_pf_bridge_target.md
outputs/jensen_window_pf_bridge_obligations.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/signed_hankel_jensen_dependency_graph.md
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py
python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the ansatz workbench is internally consistent with
the exact low-degree algebra and countermodel gates. It does not prove that any
kernel, network, production matrix, determinant integral, or sign-regularity
transfer theorem exists.
