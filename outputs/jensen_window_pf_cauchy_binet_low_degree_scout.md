# Jensen-Window PF Cauchy-Binet Low-Degree Scout

Date: 2026-07-06

Status: symbolic theorem-search scout. This is not a proof of a positive
kernel, Cauchy-Binet identity, Jensen-window PF-infinity, Jensen
hyperbolicity, Laguerre-Polya membership, RH, or `Lambda <= 0`; it tests how
much of the low-degree Jensen-window algebra is already explained by adjacent
log-concavity.

## Purpose

The live ansatz:

```text
ansatz_02_positive_cauchy_binet_kernel
```

needs a uniform identity for every degree, shift, and Toeplitz minor shape.
Before searching for a full kernel, this scout asks a smaller question:

```text
Do the degree-2/3/4 hard formulas require a kernel explanation, or are they
already certified by adjacent log-concavity of a generic positive sequence?
```

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
```

Current result:

```text
validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found
```

## Ratio Parameterization

For positive adjacent ratios, write:

```text
a0 = A
a1 = A*rho
a2 = A*rho**2*x1
a3 = A*rho**3*x1**2*x2
a4 = A*rho**4*x1**3*x2**2*x3
```

The variables `x1,x2,x3` are adjacent ratio contractions. Positive
log-concavity gives:

```text
0 <= x1,x2,x3 <= 1
```

For each hard formula, the scout extracts a positive monomial in
`A,rho,x1,x2,x3`, then certifies the remaining normalized polynomial on the
unit cube by nonnegative Bernstein coefficients.

## Finding

All selected hard formulas pass this weak certificate:

```text
15 formula rows
15 rows with nonnegative Bernstein coefficients
14 rows with strictly positive Bernstein coefficients
1 boundary row: the degree-2 Jensen threshold
cauchy_binet_identity_found=false
kernel_identity_found=false
target_closing=false
```

The degree-2 row has a zero Bernstein boundary coefficient because equality is
allowed at the log-linear boundary. The degree-3 and degree-4 selected minors
have strictly positive Bernstein coefficients after the positive monomial is
removed.

## Consequence For The Ansatz

This is useful but slightly deflationary: the current low-degree hard tests do
not force a positive Cauchy-Binet kernel. They can already be discharged by
adjacent log-concavity of a generic positive sequence.

So the live ansatz must look beyond these selected low-degree formulas. A real
kernel proof still has to explain:

```text
all degrees d
all shifts n
all Toeplitz minor shapes
the binomial weights binom(d,j)
the actual zeta heat-flow coefficient structure
```

## Countermodel Warning

The finite countermodel from the exact algebra gate also lies inside the
adjacent-log-concavity ratio box:

```text
ratio_contractions_in_unit_interval=true
adjacent_log_concavity_gaps_nonnegative=true
selected_low_degree_minors_positive=true
d3_first_negative_contiguous_toeplitz_minor_size=8
d4_first_negative_contiguous_toeplitz_minor_size=6
```

Thus adjacent log-concavity and selected low-degree Jensen-window positivity
are not enough. The countermodel passes those checks but still has negative
larger contiguous Jensen-window Toeplitz minors.

## Frontier Extension

The larger contiguous-minor frontier is recorded in:

```text
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
python work/rh_compute/scripts/jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
```

Current frontier result:

```text
validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues
```

It locates the first simple Bernstein-certificate failures at degree 3 size 6
and degree 4 size 5, and the first rational countermodel negatives at degree
3 size 8 and degree 4 size 6.

## Integration Points

This scout sharpens:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/jensen_window_pf_log_concavity_frontier_scout.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/jensen_window_pf_obligation_algebra.md
work/rh_compute/results/jensen_window_pf_obligation_algebra.json
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the low-degree Cauchy-Binet scout is reproducible
and has not overclaimed its role. It does not construct a positive kernel, a
planar network, a production matrix, a determinant integral, or any all-order
sign-regularity-to-Jensen-window PF theorem.
