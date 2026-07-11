# Jensen-Window PF Ratio-Condition Scout

Date: 2026-07-06

Status: ratio-condition theorem-search diagnostic. This is not a proof of a
positive kernel, Cauchy-Binet identity, Jensen-window PF-infinity, Jensen
hyperbolicity, Laguerre-Polya membership, RH, or `Lambda <= 0`; it tests
natural strengthened ratio conditions against the exact rational countermodel.

## Purpose

The log-concavity frontier scout shows that adjacent log-concavity is too
weak. This diagnostic asks which nearby strengthened ratio conditions are also
too weak by exact countermodel or constructed extension, and which are merely
tautological window conditions.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
```

Current result:

```text
validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction
```

## Countermodel Ratio Point

For the exact rational countermodel:

```text
x1 = 65/66
x2 = 44/65
x3 = 50/77
```

The following quantities are positive:

```text
1 - x1
1 - x2
1 - x3
x1 - x2
x2 - x3
(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)
```

So the countermodel satisfies adjacent log-concavity, decreasing ratio
contractions, and second-order log-concavity of the five-term prefix. It still
has negative larger contiguous minors:

```text
degree 3 size 8
degree 4 size 6
```

## Rejected Conditions

The checker records four rejected bridge candidates:

```text
rc_01_adjacent_log_concavity:
  rejected_by_countermodel

rc_02_decreasing_ratio_contractions:
  rejected_by_countermodel

rc_03_second_order_log_concavity:
  rejected_by_countermodel

rc_04_selected_low_degree_bernstein_positivity:
  rejected_by_countermodel
```

These are useful diagnostics, but they cannot be promoted into `jwpf_06`.

## Tautological Conditions

The checker also separates window-level conditions that would block the
countermodel but only by restating the target:

```text
rc_05_degree3_discriminant:
  tautological_window_condition

rc_06_degree4_discriminant:
  tautological_window_condition
```

They are allowed as finite consequences after a proof, not as bridge theorems.

## Constructed-Extension Rejection

The original countermodel does not satisfy the stronger contraction condition:

```text
rc_07_contraction_log_concavity:
  x2^2 >= x1*x3
```

At the countermodel:

```text
x2^2 - x1*x3 = -1946249/10735725
```

The follow-up contraction-log-concavity scout rejects it by a positive
constructed extension:

```text
x1 = 65/66
x2 = 44/65
x3 = 1/3
x2^2 - x1*x3 = 108703/836550
```

That extension satisfies `x2^2 >= x1*x3` but still has negative frontier
minors at degree 3 size 8 and degree 4 size 6. The row is therefore recorded
as:

```text
rc_07_contraction_log_concavity:
  rejected_by_constructed_extension
```

## Integration Points

This scout sharpens:

```text
outputs/jensen_window_pf_log_concavity_frontier_scout.md
outputs/jensen_window_pf_contraction_log_concavity_scout.md
outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md
outputs/jensen_window_pf_structural_ansatz_matrix.md
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json
work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json
work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the ratio-condition audit is reproducible and keeps
rejected and tautological candidates separated. It does not prove a
positive kernel, planar network, production matrix, determinant integral, or
all-order sign-regularity-to-Jensen-window PF theorem.
