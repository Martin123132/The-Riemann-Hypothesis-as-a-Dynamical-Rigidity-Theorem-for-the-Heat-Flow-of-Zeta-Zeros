# Jensen-Window PF Contraction-Log-Concavity Scout

Date: 2026-07-06

Status: ratio-condition rejection diagnostic. This is not a proof of a
positive kernel, Cauchy-Binet identity, Jensen-window PF-infinity, Jensen
hyperbolicity, Laguerre-Polya membership, RH, or `Lambda <= 0`; it tests the
last open ratio condition from the preceding scout.

## Purpose

The ratio-condition scout left one nearby strengthened condition unresolved:

```text
rc_07_contraction_log_concavity:
  x2^2 >= x1*x3
```

The old rational countermodel had `x2^2 - x1*x3 < 0`, so it could not reject
this condition directly. This diagnostic keeps the same bad degree-3 prefix
and chooses a smaller positive `x3`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_contraction_log_concavity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py
```

Current result:

```text
validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows
```

## Constructed Extension

Use the ratio parametrization

```text
a0 = A
a1 = A*rho
a2 = A*rho^2*x1
a3 = A*rho^3*x1^2*x2
a4 = A*rho^4*x1^3*x2^2*x3
```

with

```text
rho = 33/40
x1 = 65/66
x2 = 44/65
x3 = 1/3
```

This gives the positive rational prefix:

```text
A0 = 1
A1 = 33/40
A2 = 429/640
A3 = 4719/12800
A4 = 17303/256000
```

The following quantities are positive:

```text
1 - x1 = 1/66
1 - x2 = 21/65
1 - x3 = 2/3
x1 - x2 = 1321/4290
x2 - x3 = 67/195
(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3) = 3793/38025
x2^2 - x1*x3 = 108703/836550
```

So this extension satisfies adjacent log-concavity, decreasing ratio
contractions, second-order log-concavity, and contraction log-concavity.

## Frontier Failure

Despite satisfying `x2^2 >= x1*x3`, the constructed extension still has
negative larger contiguous minors:

```text
degree 3 size 8:
  -435846079534239/104857600000000

degree 4 size 6:
  -26359418151/4096000000
```

The degree-3 value matches the old countermodel exactly because the
constructed extension preserves `A0..A3`. The degree-4 value is also negative,
so the rejection is not confined to the prefix-only obstruction.

## Consequence

`rc_07_contraction_log_concavity` is rejected as a standalone bridge condition.
It may still appear inside a stronger theorem, but it cannot by itself close
the missing Jensen-window PF bridge.

## Integration Points

This scout sharpens:

```text
outputs/jensen_window_pf_ratio_condition_scout.md
outputs/jensen_window_pf_log_concavity_frontier_scout.md
work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json
work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json
work/rh_compute/results/proof_claim_ledger.json
python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
```

## Boundary

Passing this checker means the contraction-log-concavity candidate is
reproducibly rejected as a standalone ratio bridge. It does not prove a
positive kernel, planar network, production matrix, determinant integral, or
all-order sign-regularity-to-Jensen-window PF theorem.
