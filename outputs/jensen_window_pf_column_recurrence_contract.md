# Jensen-Window PF Column Recurrence Contract

Date: 2026-07-06

Status: column-shape recurrence diagnostic. This is not a proof of Schur
positivity, Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya
membership, RH, or `Lambda <= 0`; it isolates a necessary column-shape
condition for any positive Schur, network, production-matrix, or
determinant-integral route.

## Purpose

The Schur shape contract identifies the hard frontier minors as column-shape
Jacobi-Trudi determinants. Column shapes are special: after normalizing

```text
g_j = h_j / h_0
h_j = binom(d,j) * A_{n+j}
```

the column determinant is:

```text
C_m = h_0^m * e_m
```

where `e_m` is the elementary-symmetric coefficient determined by:

```text
E(t) = 1 / H(-t)
H(t) = 1 + g_1*t + ... + g_d*t^d
```

So the column-shape part of the Jensen-window PF bridge requires:

```text
C_m >= 0
for every m >= 0, every degree d, and every shift n.
```

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_column_recurrence_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
```

Current result:

```text
validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows
```

## Recurrence

For the unnormalized column determinants:

```text
C[0] = 1
C[m<0] = 0
C[m] = sum_{j=1..min(d,m)} (-1)^(j-1) h0^(j-1) h_j C[m-j]
```

The first degrees are:

```text
d=2:
  C[m] = h1*C[m-1] - h0*h2*C[m-2]

d=3:
  C[m] = h1*C[m-1] - h0*h2*C[m-2] + h0^2*h3*C[m-3]

d=4:
  C[m] = h1*C[m-1] - h0*h2*C[m-2] + h0^2*h3*C[m-3] - h0^3*h4*C[m-4]

d=5:
  C[m] = h1*C[m-1] - h0*h2*C[m-2] + h0^2*h3*C[m-3] - h0^3*h4*C[m-4] + h0^4*h5*C[m-5]
```

This is a necessary Schur-positivity condition. It is not sufficient for all
skew shapes.

## Hard Frontier Rows

The recurrence exactly reproduces the two hard column-shape frontier
determinants already recorded in the Schur shape contract and obligation
algebra.

```text
d3_column_recurrence_m8:
  source shape row = d3_column_shape_m8
  countermodel value = -435846079534239/104857600000000

d4_column_recurrence_m6:
  source shape row = d4_column_shape_m6
  countermodel value = -229760849637/28672000000
```

Both values are negative at the exact rational log-concave countermodel.
Therefore generic adjacent log-concavity, contraction conditions, or selected
low-degree positivity cannot prove the column-shape recurrence condition.

The finite zeta-window coverage map is:

```text
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py
```

It validates:

```text
validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues
```

This is finite evidence only. It does not close the all-order recurrence
theorem.

The theorem-mechanism route matrix for this recurrence target is:

```text
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
```

It classifies real-root/PF endpoint language, positive renewal or convolution
kernels, positive continued fractions, production matrices, Kaluza-type sign
theorems, generic ratio conditions, finite recurrence stress, and the missing
all-Schur lift. It records `0` ready-to-apply rows.

## Consequence

The next viable structural theorem must explain why the actual zeta
heat-flow windows make this recurrence nonnegative for all `m`, `d`, and
`n`, or must bypass this language with an equivalent positive construction.

A proof route that cannot establish the recurrence sequence

```text
C_m(h_0,...,h_d) >= 0
```

for the actual Jensen windows cannot prove Schur positivity, and therefore
cannot close the Jensen-window PF bridge through this route.

## Integration Points

This diagnostic sharpens:

```text
outputs/jensen_window_pf_schur_shape_contract.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/jensen_window_pf_log_concavity_frontier_scout.md
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json
work/rh_compute/results/jensen_window_pf_schur_shape_contract.json
work/rh_compute/results/jensen_window_pf_obligation_algebra.json
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the column-shape recurrence and the two hard
frontier recurrence rows are reproducible. It does not prove Schur positivity,
a positive specialization, a planar network, a production matrix, a determinant
integral, or an all-order sign-regularity-to-Jensen-window PF theorem.
