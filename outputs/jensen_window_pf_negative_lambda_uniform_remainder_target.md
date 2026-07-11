# Jensen-Window PF Negative-Lambda Uniform Remainder Target

Date: 2026-07-06

Status: open theorem target. This is not a proof of the bounded
log-curvature theorem, cone entry, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_uniform_remainder_target`.

Proof boundary: this artifact states the two-scale uniform remainder
theorem needed by the signed perturbation route. It does not prove that
theorem or close the negative-lambda tail.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_uniform_remainder_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_uniform_remainder_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_uniform_remainder_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda uniform remainder target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 open requirements, 3 leading-scale rows
```

## Scale Obstruction

The fixed-`k` signed perturbation gives a leading deficit:

```text
B_k ~ D2/T^2
D2 = [190.68658368200 +/- 5.53e-12]
```

but the target bound shrinks with `k`:

```text
B_k <= 2/(3*(2*k+1))
```

Therefore a fixed-`k` expansion cannot be the all-tail theorem by itself.

Leading-only scale diagnostics:

```text
T=25: k_max~5.925432157343440434E-1 (floor 0), tail starts at k=22
T=50: k_max~3.870172862937376174E+0 (floor 3), tail starts at k=22
T=100: k_max~1.698069145174950469E+1 (floor 16), tail starts at k=22
```

These are scale diagnostics only. They ignore all remainders and do not
contradict the certified finite prefix.

## Required Theorem

A proof must provide one of the following:

```text
1. a uniform local/mesoscopic remainder theorem, plus a global far-tail saddle theorem
2. a finite collar to K, plus an analytic tail theorem from K+1 onward
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```

Summary:

The signed Gaussian perturbation route now needs a two-scale uniform remainder theorem: the fixed-k Taylor signs give the correct local direction, but the bounded-curvature target shrinks like 1/k, so a leading fixed-k D2/T^2 estimate cannot be promoted to the all-k tail. A proof must either control the moving saddle globally in k/T or combine a finite collar with a separate far-tail theorem.
