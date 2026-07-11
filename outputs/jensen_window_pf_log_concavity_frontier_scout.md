# Jensen-Window PF Log-Concavity Frontier Scout

Date: 2026-07-06

Status: symbolic frontier diagnostic. This is not a proof of a positive kernel,
Cauchy-Binet identity, Jensen-window PF-infinity, Jensen hyperbolicity,
Laguerre-Polya membership, RH, or `Lambda <= 0`; it locates where adjacent
log-concavity stops controlling contiguous Jensen-window Toeplitz minors.

## Purpose

The previous Cauchy-Binet low-degree scout showed that the selected degree-2,
degree-3, and degree-4 hard formulas are already certified by adjacent
log-concavity. This frontier scout asks where that weak explanation breaks.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
```

Current result:

```text
validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues
```

## Method

For the same ratio parameterization:

```text
a0 = A
a1 = A*rho
a2 = A*rho**2*x1
a3 = A*rho**3*x1**2*x2
a4 = A*rho**4*x1**3*x2**2*x3
```

the scout computes contiguous Jensen-window Toeplitz determinants:

```text
degree 3, minor sizes m = 1..8
degree 4, minor sizes m = 1..6
```

For each row it extracts a positive monomial and checks Bernstein coefficients
of the normalized ratio polynomial on the log-concavity box
`0 <= x1,x2,x3 <= 1`.

## Frontier

The exact frontier is:

```text
degree 3 size 6:
  first negative Bernstein coefficient

degree 3 size 8:
  first negative value at the rational log-concave countermodel

degree 4 size 5:
  first negative Bernstein coefficient

degree 4 size 6:
  first negative value at the rational log-concave countermodel
```

The scout records:

```text
kernel_identity_found=false
target_closing=false
```

## Interpretation

Adjacent log-concavity explains the selected low-degree minors but does not
control the larger contiguous minors. In particular, the exact rational
countermodel from:

```text
work/rh_compute/results/jensen_window_pf_obligation_algebra.json
```

is still inside the ratio box and has:

```text
d3_first_negative_contiguous_toeplitz_minor.size = 8
d4_first_negative_contiguous_toeplitz_minor.size = 6
```

So a valid Cauchy-Binet/kernel route must add structure beyond adjacent
log-concavity and beyond the selected low-degree formulas.

## Ratio-Condition Audit

The nearby strengthened-ratio candidates are audited here:

```text
outputs/jensen_window_pf_ratio_condition_scout.md
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
python work/rh_compute/scripts/jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
```

Current ratio-condition result:

```text
validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction
```

This rejects adjacent log-concavity, decreasing contraction variables,
second-order log-concavity, and selected low-degree Bernstein positivity as
standalone bridge conditions.

## Integration Points

This frontier scout sharpens:

```text
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
outputs/jensen_window_pf_obligation_algebra.md
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the log-concavity frontier is reproducible and
does not overclaim. It does not prove a positive kernel, a planar network, a
production matrix, a determinant integral, or an all-order
sign-regularity-to-Jensen-window PF theorem.
