# Jensen-Window PF Negative-Lambda Adaptive Scaled-Defect Target

Date: 2026-07-06

Status: open theorem target. This is not a proof of zeta cone entry,
Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_adaptive_scaled_defect_target`.

Proof boundary: this artifact replaces the finite-rejected fixed
half-width route with an adaptive or exact-cone scaled-defect target.
It does not prove that target or the separate monotone defect bridge.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda adaptive scaled-defect target: 8 rows, 0 issues, 2 live routes, 597 exact-cone rows, 76 half-width failures
```

## Target

```text
d_k = 1 - x_k
s_k = ((2*k+1)/2) * d_k
prove 0 <= s_k <= 1 for all k >= 200
prove d_(k+1) <= d_k for all k >= 199 separately
optionally sharpen 0 <= s_k <= 1 to an explicit adaptive envelope E_lambda(k)<1
```

## k200 Frontier

```text
exact cone rows: 597 / 597
fixed half-width rows: 521 / 597
fixed half-width failure rows: 76
one-third failure rows: 418
max scaled defect: 5.376643171065356005E-1 at lambda=-25.0, k=199
min exact-cone margin: 1.346067912947116963E-2 at lambda=-100.0, k=1
next finite coefficient needed without an analytic bridge: A_201
```

First fixed half-width failures:

```text
lambda=-50.0: first failure k=191, failure rows=9, max s=5.049575716217948385E-1 at k=199
lambda=-25.0: first failure k=133, failure rows=67, max s=5.376643171065356005E-1 at k=199
```

Live routes:

```text
1. uniform saddle/Laplace control with interval-safe remainders
2. direct ratio-recurrence or comparison inequalities compatible with increasing s_k
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```

Summary:

The k200 frontier keeps the exact scaled cone 0<=s_k<=1 on 597/597 rows, while rejecting the fixed half-width buffer on 76 rows and the one-third buffer on 418 rows. The replacement target is an exact-cone or adaptive-envelope all-k theorem from k>=200, plus the separate monotone defect bridge.
