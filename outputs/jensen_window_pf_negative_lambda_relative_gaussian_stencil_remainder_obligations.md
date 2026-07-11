# Jensen-Window PF Negative-Lambda Relative-Gaussian Stencil Remainder Obligations

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of a uniform Taylor-tail remainder theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations`.

Proof boundary: this artifact decomposes the Taylor-tail remainder into
three exact epsilon stencils and records finite margin budgets. It does
not prove any uniform all-`k` remainder bound.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: 9 rows, 0 issues, 4 positive baseline rows, 31 blocked baseline rows, 4 exact stencil rows, 0 ready-to-apply rows
```

## Exact Epsilon Stencils

Write:

```text
f_k = f_k^(M)+epsilon_k
```

Then the three required stencil errors are:

```text
E_B = 2*epsilon_k-epsilon_(k-1)-epsilon_(k+1)
E_U = -epsilon_(k-1)+3*epsilon_k-3*epsilon_(k+1)+epsilon_(k+2)
E_C = (2*k+1)*epsilon_(k-1)-(6*k+5)*epsilon_k+(6*k+7)*epsilon_(k+1)-(2*k+3)*epsilon_(k+2)
```

## Finite Margin Budget

```text
truncation rows: 35
positive baseline rows: 4
blocked baseline rows: 31
invalid normalizers: 9
upper-wall violations: 2
companion failures: 8
weighted-gap failures: 12
weakest positive margin: 2.332981128843164884E-8 at nlrgts_M6_T2000
weakest half-margin budget: 1.166490564421582442E-8 at nlrgts_M6_T2000
```

The half-margin budget is a finite-row sufficient target only. A proof
must replace it with uniform q/T bounds for `E_B`, `E_U`, and `E_C`.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
```

Summary:

The relative-Gaussian remainder problem can be split into three exact epsilon stencils for B_k, B_k-B_(k+1), and C_(k+1)-C_k. In the finite degree<=14 truncation matrix, only 4/35 rows have positive baseline margins for all three stencils; the weakest half-margin budget is 1.166490564421582442E-8, so a promoted proof needs sharp uniform epsilon-stencil bounds rather than finite truncation promotion.
