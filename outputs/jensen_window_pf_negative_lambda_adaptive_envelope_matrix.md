# Jensen-Window PF Negative-Lambda Adaptive Envelope Matrix

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
adaptive scaled-defect target, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_adaptive_envelope_matrix`.

Proof boundary: this artifact records checked monotone-envelope
patterns only. It does not prove all-`k` bounds, continuous lambda
monotonicity, or a tail theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda adaptive envelope matrix: 7 matrix rows, 0 issues, 594 k-increase rows, 398 lambda-order rows, 76 half-width failures
```

## Matrix

```text
lambdas: -25.0, -50.0, -100.0
coefficient range: A_0..A_200
checked contractions: x_1..x_199
exact cone rows: 597 / 597
k-increase rows: 594 / 594
lambda-order rows: 398 / 398
half-width failure rows: 76
one-third failure rows: 418
```

Cross-lambda order gaps:

```text
-25.0>=-50.0: 199 / 199 positive, min gap 2.185649379581140840E-2 at k=1
-50.0>=-100.0: 199 / 199 positive, min gap 1.663609726220394924E-2 at k=1
```

Extrema:

```text
max scaled defect: 5.376643171065356005E-1 at lambda=-25.0, k=199
min upper cone slack: 4.623356828934643995E-1 at lambda=-25.0, k=199
min exact-cone margin: 1.346067912947116963E-2 at lambda=-100.0, k=1
min k-increase: 4.449655594664470276E-4 at lambda=-25.0, k=198
min lambda-order gap: 1.663609726220394924E-2 for -50.0>=-100.0, k=1
```

Live monotone-envelope route:

```text
prove s_k(lambda) increases in k without crossing 1
prove s_k(lambda) decreases as |lambda| increases on the needed negative-lambda ray
prove a limiting/adaptive upper envelope below the exact cone wall
keep the one-third, half-width, and nonincrease shortcuts rejected
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md
outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```

Summary:

On the k200 negative-lambda prefix, s_k is interval-increasing in k on all 594/594 adjacent rows and interval-ordered by |lambda| on all 398/398 cross-lambda rows. The largest checked scaled defect is 5.376643171065356005E-1 at lambda=-25.0, k=199, leaving upper cone slack 4.623356828934643995E-1. This supports a monotone-envelope proof search, but supplies no all-k tail theorem.
