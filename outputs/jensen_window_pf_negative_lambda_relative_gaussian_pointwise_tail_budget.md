# Jensen-Window PF Negative-Lambda Relative-Gaussian Pointwise Tail Budget

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of a uniform Taylor-tail remainder theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget`.

Proof boundary: this artifact converts epsilon-stencil margins into
pointwise log-tail and multiplicative relative-tail budgets. It does
not prove the analytic tail estimates required by those budgets.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: 9 rows, 0 issues, 4 positive baseline rows, 4 budget rows, 0 ready-to-apply rows
```

## Pointwise Envelope

Assume `|epsilon_j|<=eta_j` for `j=k-1,k,k+1,k+2`. Then:

```text
|E_B| <= eta_(k-1)+2*eta_k+eta_(k+1)
|E_U| <= eta_(k-1)+3*eta_k+3*eta_(k+1)+eta_(k+2)
|E_C| <= (2*k+1)*eta_(k-1)+(6*k+5)*eta_k+(6*k+7)*eta_(k+1)+(2*k+3)*eta_(k+2)
```

For a uniform pointwise envelope `eta_j<=eta`, the coefficient sums are:

```text
B sum: 4
companion sum: 8
weighted-gap sum: 16*k+16
at k=22, weighted-gap sum: 368
```

## Current Finite Bottleneck

```text
positive baseline rows: 4
blocked baseline rows: 31
budget rows: 4
limiting stencil counts: {'companion': 3, 'weighted_gap': 1}
weakest half-safety eta: 1.458113205526978052E-9 at nlrgts_M6_T2000
weakest limiting stencil: companion
relative tail ratio bound: 1.458113204463930993E-9
```

The conversion from multiplicative relative tail to log-tail envelope is:

```text
F_j = F_j^(M)*(1+theta_j)
|theta_j| <= rho < 1
|epsilon_j| = |log(1+theta_j)| <= -log(1-rho)
rho <= 1-exp(-eta)
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md
```

Summary:

The exact epsilon-stencil obligations reduce to a concrete pointwise log-tail target. At k=22, the absolute coefficient sums are 4, 8, and 368; among the 4 positive finite baselines, the companion stencil is limiting in 3 rows and the weighted-gap stencil is limiting in 1 row. The weakest half-safety log-tail envelope is 1.458113205526978052E-9 at nlrgts_M6_T2000, equivalent to the multiplicative relative-tail ratio bound 1.458113204463930993E-9. This sharpens the open uniform remainder theorem but does not prove any analytic tail estimate.
