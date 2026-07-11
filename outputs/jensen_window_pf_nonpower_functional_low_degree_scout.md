# Jensen-Window PF Nonpower Functional Low-Degree Scout

Date: 2026-07-06

Status: exact low-degree nonpower-functional scout. This is not a proof of
the all-order column recurrence, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The nonordinary positive-transform ansatz matrix leaves:

```text
npt_04_nonpower_positive_functional
```

live only if someone constructs a positive functional `L` and non-power basis
`K_m` such that:

```text
L(K_m) = mu_m = [t^m]1/H(-t)
```

This scout makes the first algebraic traps explicit. For:

```text
H(t) = 1 + g_1*t + ... + g_d*t^d
E(t) = 1/H(-t)
```

the reciprocal coefficients obey:

```text
mu_m = sum_{1<=j<=m} (-1)^(j+1) g_j mu_{m-j}
```

and equivalently:

```text
mu_m = sum_{alpha composition of m}
       (-1)^(m-len(alpha)) prod_i g_{alpha_i}
```

Machine-readable scout:

```text
work/rh_compute/results/jensen_window_pf_nonpower_functional_low_degree_scout.json
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_low_degree_scout.py
```

Current result:

```text
validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows
```

## Exact Formulas

The checker recomputes:

```text
mu_0 = 1
mu_1 = g1
mu_2 = g1^2 - g2
mu_3 = g1^3 - 2*g1*g2 + g3
mu_4 = g1^4 - 3*g1^2*g2 + g2^2 + 2*g1*g3 - g4
mu_5 = g1^5 - 4*g1^3*g2 + 3*g1*g2^2 + 3*g1^2*g3 - 2*g2*g3 - 2*g1*g4 + g5
mu_6 = g1^6 - 5*g1^4*g2 + 6*g1^2*g2^2 + 4*g1^3*g3 - g2^3 - 6*g1*g2*g3 - 3*g1^2*g4 + g3^2 + 2*g2*g4 + 2*g1*g5 - g6
```

For `m>=2`, the raw composition expansion has both signs:

```text
m=2: 1 positive, 1 negative
m=3: 2 positive, 2 negative
m=4: 4 positive, 4 negative
m=5: 8 positive, 8 negative
m=6: 16 positive, 16 negative
m=7: 32 positive, 32 negative
m=8: 64 positive, 64 negative
```

## Scout Rows

```text
npf_01_exact_reciprocal_polynomials:
  exact_identity
  The first reciprocal coefficients are signed polynomials in g_j.

npf_02_signed_composition_expansion:
  exact_identity
  The composition expansion has mixed signs for every m>=2.

npf_03_positive_g_cone_obstruction:
  rejected_as_generic_positive_cone
  Positivity cannot follow from independent nonnegativity of the g_j.

npf_04_degree2_log_concavity_entry:
  first_required_inequality
  Since g_j=h_j/h_0, mu_2=(h_1^2-h_0*h_2)/h_0^2; a functional proof must
  import at least this adjacent log-concavity numerator.

npf_05_degree3_cancellation_entry:
  higher_cancellation_required
  mu_3=(h_1^3-2*h_0*h_1*h_2+h_0^2*h_3)/h_0^3 requires a genuine signed
  cancellation certificate.

npf_06_tautological_basis_trap:
  rejected_as_tautological
  A basis or functional that already contains the unknown mu_m is not a
  proof.

npf_07_positive_cone_contract:
  live_if_cone_and_functional_constructed
  A serious route must specify a cone C, prove K_m in C, prove L positive on
  C, and prove L(K_m)=mu_m exactly for all m,d,n. The cone candidate matrix
  rejects raw g-coordinate, standalone ratio/log-concavity, tautological, and
  endpoint PF/LP cones, leaving only Xi/Phi kernel and
  Cauchy-Binet/determinant-integral cone routes live.
```

## Evidence Anchors

```text
outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md
outputs/jensen_window_pf_positive_readout_theorem_target.md
outputs/jensen_window_pf_positive_spectral_moment_obstruction.md
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
outputs/jensen_window_pf_reciprocal_fraction_scout.md
outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_contraction_log_concavity_scout.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md
outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md
```

The cone candidate matrix validates:

```text
validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows
```

## Kill Gates

Reject a proposed nonpower functional if it:

```text
1. uses only independent nonnegativity of the g_j;
2. treats the signed composition expansion as termwise positive;
3. defines K_m or L using the unknown target coefficients mu_m;
4. proves transformed positivity but not L(K_m)=mu_m for the original
   reciprocal coefficients;
5. promotes finite low-degree or finite-grid positivity to an all-m,d,n
   theorem;
6. assumes endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH,
   or Lambda <= 0.
```

## Boundary

Passing this checker means the `npt_04` nonpower-functional route now has
exact low-degree algebra, signed-composition traps, and a positive-cone
acceptance contract. It does not construct the cone `C`, the basis `K_m`, or
the functional `L`, and it does not prove reciprocal coefficient positivity
or `Lambda <= 0`.
