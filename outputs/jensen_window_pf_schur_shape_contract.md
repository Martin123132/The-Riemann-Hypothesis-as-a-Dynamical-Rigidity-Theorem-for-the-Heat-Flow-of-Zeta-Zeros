# Jensen-Window PF Schur Shape Contract

Date: 2026-07-06

Status: exact shape-contract diagnostic. This is not a proof of Schur
positivity, Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya
membership, RH, or `Lambda <= 0`; it makes the structural target for positive
Schur, network, or determinant-integral routes explicit.

## Purpose

For a fixed Jensen-window degree `d`, define the finite-support specialization:

```text
h_j = binom(d,j) * A_{n+j}
h_j = 0 for j < 0 or j > d
```

Then each finite-band Toeplitz minor

```text
det(h_{q_j-r_i})
```

is a Jacobi-Trudi skew-Schur determinant

```text
s_{lambda/mu}(h_0,...,h_d)
```

after the same row/column reindexing used in the coefficient-PF Schur target.
This contract records which bounded shapes the Jensen-window bridge must cover
and highlights the hard frontier column shapes.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_schur_shape_contract.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
```

Current result:

```text
validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows
```

## Bounded Shape Grid

The bounded diagnostic uses `N=8`, orders `<=5`, and degrees `d=2,3,4,5`.
For each degree it scans `12,020` row/column minors.

```text
d=2:
  finite-band structural nonzero = 2907
  finite-band zero after upper support = 1581
  unique finite-band shapes = 1848

d=3:
  finite-band structural nonzero = 3960
  finite-band zero after upper support = 528
  unique finite-band shapes = 2654

d=4:
  finite-band structural nonzero = 4370
  finite-band zero after upper support = 118
  unique finite-band shapes = 2989

d=5:
  finite-band structural nonzero = 4472
  finite-band zero after upper support = 16
  unique finite-band shapes = 3078
```

Across these four rows:

```text
total tested minors = 48080
total finite-band structural nonzero = 15709
```

So a structural proof cannot merely cover contiguous windows or a handful of
selected formulas. Even a small bounded grid already contains thousands of
finite-band skew shapes.

## Hard Frontier Column Shapes

The earlier log-concavity frontier scout found exact failures at degree 3 size
8 and degree 4 size 6 for generic log-concave countermodels. In Schur shape
language these are column-shape Jacobi-Trudi determinants.

```text
d3_column_shape_m8:
  lambda = (1,1,1,1,1,1,1,1)
  mu = (0,0,0,0,0,0,0,0)
  term count = 10

d4_column_shape_m6:
  lambda = (1,1,1,1,1,1)
  mu = (0,0,0,0,0,0)
  term count = 9
```

Both hard frontier determinants have mixed-sign h-monomial expansions. This
means termwise monomial positivity cannot be the explanation; a successful
positive Schur, network, or determinant-integral route, including any
production-matrix formulation, must provide independent structural positivity
for these column shapes and then for all remaining shapes.

The column-shape recurrence refinement is:

```text
outputs/jensen_window_pf_column_recurrence_contract.md
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
```

It validates:

```text
validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows
```

This records the elementary-symmetric recurrence `C_m = h_0^m * e_m` and
shows that both hard frontier recurrence rows are negative at the exact
rational countermodel.

## Consequence

This contract gives `jwpf_06` a sharper structural target:

```text
prove s_{lambda/mu}(h_0,...,h_d) >= 0
for h_j = binom(d,j) * A_{n+j}(0)
for every d, n, and finite-band Toeplitz shape
```

or prove an equivalent positive network, production matrix, Cauchy-Binet, or
determinant-integral theorem that implies the same inequalities.

## Integration Points

This diagnostic sharpens:

```text
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
outputs/positive_schur_specialization_target.md
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
work/rh_compute/results/jensen_window_pf_schur_shape_contract.json
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the bounded Schur/Jacobi-Trudi shape contract is
reproducible and the hard frontier shapes are recorded. It does not prove
Schur positivity, a positive specialization, a planar network, a production
matrix, a determinant integral, or an all-order sign-regularity-to-Jensen-window
PF theorem.
