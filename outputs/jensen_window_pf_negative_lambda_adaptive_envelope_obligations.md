# Jensen-Window PF Negative-Lambda Adaptive Envelope Obligations

Date: 2026-07-06

Status: exact algebraic obligation diagnostic. This is not a proof of
the adaptive scaled-defect target, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_adaptive_envelope_obligations`.

Proof boundary: this artifact translates the adaptive route into exact
ratio inequalities and separates exact inputs from open requirements.
It does not prove the open requirements.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py
```

Current result:

```text
validated Jensen-window PF negative-lambda adaptive envelope obligations: 9 obligation rows, 0 issues, 3 exact rows, 3 open requirements, 1 rejected routes
```

## 2026-07-10 Upper-Wall Handoff

`outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md` proves
`x_k<=1`, so the obligation `0<=s_k` is no longer open. The algebraic ledger
below remains useful as a historical decomposition, while the current hard
requirement is the adjacent-`k` monotone bridge and its scaled-envelope
strengthening.

## 2026-07-11 Lambda=-100 Closure Handoff

The full cone-entry and raw-corridor theorems now discharge every obligation
listed below at `lambda=-100`: both pointwise walls, the adjacent monotone
bridge, and scaled-defect increase. The simultaneous all-three-lambda
formulation remains historical rather than proved.

## Exact Translations

```text
x_k=(A_{k-1}*A_{k+1})/A_k^2
d_k=1-x_k
s_k=((2*k+1)/2)*d_k
1-s_k=((2*k+1)*x_k-(2*k-1))/2
s_k>=0 iff x_k<=1
d_k-d_(k+1)=x_(k+1)-x_k
s_(k+1)-s_k=(2+(2*k+1)*x_k-(2*k+3)*x_(k+1))/2
```

## Obligations

```text
already exact: x_k >= (2*k-1)/(2*k+1) from the boundary-threshold lemma
still open: x_k <= 1 on the needed tail
still open: x_(k+1) >= x_k on the needed tail
still open: limiting/adaptive envelope E_lambda<1
rejected: fixed half-width and one-third buffers as global routes
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md
outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md
outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md
```

Summary:

The adaptive scaled-defect target decomposes exactly: the lower threshold s_k<=1 is already the boundary-threshold lemma, while 0<=s_k is x_k<=1, the monotone bridge is x_(k+1)>=x_k, and scaled k-monotonicity is 2+(2*k+1)*x_k-(2*k+3)*x_(k+1)>=0. The remaining proof burden is a noncircular upper-wall/monotone-envelope theorem compatible with the finite half-width and one-third rejection gates.
