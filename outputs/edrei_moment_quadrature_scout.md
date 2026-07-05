# Edrei Moment-Recurrence Scout

Date: 2026-07-04

Status: finite Arb recurrence diagnostic. This is not a proof of PF-infinity, Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Purpose

The Edrei-log route studies:

```text
H_lambda(z) = sum d_k(lambda) z^k
log H_lambda(z) = sum ell_n(lambda) z^n
q_n = n ell_n
p_n = (-1)^(n-1) q_n
```

If:

```text
H(z) = exp(gamma z) product_j (1 + beta_j z)
gamma >= 0
beta_j >= 0
```

then for `n >= 2`, the `p_n` behave like positive power sums of the `beta_j`. The shifted sequence:

```text
a_n = p_{n+1},  n >= 0
```

should therefore behave like moments of the positive measure:

```text
x dnu(x)
```

This scout reads the Arb Edrei-log row logs and checks finite monic orthogonal-polynomial recurrence data for the shifted moment sequence using Arb interval arithmetic. It is a constructive consistency diagnostic for the Edrei zero-parameter representation target.

## Scripts

```text
work/rh_compute/scripts/edrei_moment_quadrature_scout.py
work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py
```

Default promoted input rows:

```text
work/rh_compute/results/arb_edrei_log_sign_lam0_n1_n57_dps2400_boundary_tol1e-140.jsonl
work/rh_compute/results/arb_edrei_log_sign_lam1em6_n1_n57_dps2400_boundary_tol1e-140.jsonl
work/rh_compute/results/arb_edrei_log_sign_lam1em4_n1_n57_dps2400_boundary_tol1e-140.jsonl
work/rh_compute/results/arb_edrei_log_sign_lam1em2_n1_n57_dps2400_boundary_tol1e-140.jsonl
work/rh_compute/results/arb_edrei_log_sign_lam1em1_n1_n57_dps2400_boundary_tol1e-140.jsonl
```

The scout is interval-based, but it remains a finite necessary-condition diagnostic. It is not promoted to a proof of the all-order Edrei representation.

## Positive Arb Recurrence Scout

Command:

```text
python work/rh_compute/scripts/edrei_moment_quadrature_scout.py \
  --order-min 2 \
  --order-max 12 \
  --out-jsonl work/rh_compute/results/edrei_moment_recurrence_lamgrid_order2_order12_arb.jsonl \
  --out-summary work/rh_compute/results/edrei_moment_recurrence_lamgrid_order2_order12_arb_summary.json
```

Result:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
recurrence orders = 2..12
rows = 55
positive_rows = 55
negative_rows = 0
inconclusive_rows = 0
all_positive = true
```

Interpretation:

```text
The finite shifted moment sequence is constructively consistent with a
positive Edrei zero-parameter moment interpretation through recurrence
order 12.
```

Boundary:

```text
This proves only a finite recurrence diagnostic. It does not prove the
all-n moment representation, convergence of an Edrei product, PF-infinity,
or RH.
```

## Frontier Arb Recurrence Scout

Command:

```text
python work/rh_compute/scripts/edrei_moment_quadrature_scout.py \
  --order-min 2 \
  --order-max 20 \
  --out-jsonl work/rh_compute/results/edrei_moment_recurrence_lamgrid_order2_order20_arb.jsonl \
  --out-summary work/rh_compute/results/edrei_moment_recurrence_lamgrid_order2_order20_arb_summary.json
```

Result:

```text
rows = 95
positive_rows = 55
negative_rows = 0
inconclusive_rows = 40
all_positive = false
```

The frontier scout is positive through order 12 and interval-inconclusive from order 13 through order 20 for each lambda. No negative recurrence row is certified.

Interpretation:

```text
The current persisted Edrei-log enclosures support recurrence positivity
through order 12. Beyond that, the recurrence reconstruction is too
ill-conditioned for the current row precision and becomes inconclusive.
```

This replaces the earlier midpoint-only quadrature attempt, which produced spurious negative nodes from unstable Gram-Schmidt reconstruction. The corrected Arb recurrence scout records no negative rows.

## Relation To Power-Hankel Certificates

The Arb power-Hankel manifest validates many shifted determinants of:

```text
det(p_{i+j+s})_{i,j=0}^m
```

with interval-separated positive signs. The recurrence scout asks a different finite constructive question: whether monic recurrence data extracted from the shifted moments remains interval-positive.

These are consistent:

```text
power-Hankel determinant certificates: 4,205 promoted finite positives
recurrence scout: 55 finite positive recurrence rows, then 40 inconclusive rows
```

Neither gives the all-order Edrei representation. Both are finite diagnostics that sharpen the theorem target.

## Countermodel Guard

The finite-to-infinite boundary is now protected by an exact rational gate:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

Current relevant result:

```text
recurrence order <= 12
moments 0..23 preserved
positive edited moment 24
next Hankel/moment gate breaks
```

The model starts from factorial moments `m_n = n!`, which have a genuine
positive Stieltjes moment interpretation, preserves the finite prefix needed
for recurrence order `12`, and then chooses the next even moment below the
exact Schur-complement threshold. This makes the next Hankel determinant
negative while all preserved leading Hankel determinants remain positive.

Consequence:

```text
positive recurrence data through order 12 cannot be promoted to an all-order
Edrei representation without a new infinite theorem or explicit positive
parameter construction.
```

## Manifest Guard

The scout artifacts are checked by:

```text
python work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py
```

Current manifest result:

```text
validated 1 positive Arb recurrence scout and 1 inconclusive frontier scout
```

This guard exists so future notes preserve both facts:

```text
Arb recurrence data is positive through order 12
the order 13..20 frontier is inconclusive, not negative
```

without promoting either into a theorem.

## Next Useful Upgrade

The next theorem-search upgrade would be one of:

```text
1. derive recurrence coefficients analytically from the zeta structure;
2. tighten Edrei-log enclosures enough to push recurrence positivity beyond order 12;
3. prove the all-n Edrei log-power representation directly;
4. identify a stable positive parameter reconstruction from the recurrence data.
```
