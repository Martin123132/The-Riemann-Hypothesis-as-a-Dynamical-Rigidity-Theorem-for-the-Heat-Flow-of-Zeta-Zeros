# Jensen-Window PF Negative-Lambda Defect-Tail Theorem Target

Date: 2026-07-06

Status: open theorem target. This is not a proof of zeta cone entry,
Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_defect_tail_theorem_target`.

Proof boundary: this artifact states the analytic all-`k` defect-tail
theorem needed after the finite negative-lambda prefix. It does not
prove that theorem, does not prove `jwpf_06`, and does not establish
`Lambda <= 0`.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_defect_tail_theorem_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_tail_theorem_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes
```

## 2026-07-10 Upper-Wall Handoff

The interval theorem certificate
`outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md` now proves
`x_k<=1`, equivalently `d_k>=0`, for every real heat parameter and every
`k>=1`. The target below is retained in its original full form for audit
history, but its nonnegative-defect clause is discharged. The live tail burden
is now `d_k<=2/(2*k+1)` from `k>=150` and
`d_(k+1)<=d_k` from `k>=149`.

## Exact Tail Statement

Let:

```text
d_k = 1 - x_k
x_k = (A_{k+1}/A_k)/(A_k/A_{k-1})
```

The finite prefix currently supplies:

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_150
finite contractions: x_1..x_149
```

The analytic tail theorem needed to complete the infinite cone entry is:

```text
0 <= d_k <= 2/(2*k+1) for all k >= 150
d_(k+1) <= d_k for all k >= 149
```

Without such a theorem, the next finite bridge needs `A_151`.

## Live Routes

Two noncircular proof routes remain live:

```text
1. uniform saddle/Laplace control of the actual zeta heat-flow coefficients
2. direct ratio-recurrence or comparison inequalities for the coefficient tail
```

The scaled-defect nonincreasing shortcut is rejected by the finite scout:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md
validated Jensen-window PF negative-lambda tail-barrier scout: 447 cone-buffer rows, 444 defect-monotone rows, 179 one-third-buffer rows, 444 scaled-defect increase rows, 1 rejected candidate, 0 issues
```

The scaled-defect frontier now separates viable finite buffers from
over-strong ones:

```text
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues
one-third first fails at lambda=-25.0, k=31; fixed half-width first fails at lambda=-25.0, k=133
```

The older defect-recurrence scout is retained as a finite diagnostic,
but its one-third buffer should no longer be promoted as an all-tail
sufficient route after the k150 frontier:

```text
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues
0 <= d_k <= 2/(3*(2*k+1))
d_(k+1) <= d_k
```

## Integration

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k150_scout.md
outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md
outputs/jensen_window_pf_negative_lambda_finite_collar_k150_contract.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
```

Summary:

The defect-tail target is now indexed precisely: finite evidence covers x_1..x_149, so an analytic proof must supply 0 <= d_k <= 2/(2*k+1) for k >= 150 and d_(k+1) <= d_k for k >= 149. Two live routes remain: uniform saddle/Laplace control and a direct ratio-recurrence comparison. The scaled-defect nonincreasing shortcut is rejected by the finite diagnostic; the one-third buffer is too strong on the k150 prefix, and the fixed half-width buffer is also finite-rejected.
