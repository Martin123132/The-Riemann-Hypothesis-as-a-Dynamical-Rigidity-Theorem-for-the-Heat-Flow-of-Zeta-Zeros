# Jensen-Window PF Negative-Lambda Bounded Log-Curvature Target

Date: 2026-07-06

Status: open theorem target. This is not a proof of the defect-tail
theorem, cone entry, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_bounded_log_curvature_target`.

Proof boundary: this artifact states the missing quantitative
log-curvature estimate. It does not prove that estimate and does not
close the negative-lambda tail.

Retirement notice: this fixed `2/3` scaled-curvature wall is no longer a
live theorem target after the repaired k300 obstruction in
`outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md`.
The note is retained as the historical k<=22 compatibility target.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda bounded log-curvature target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes, 63 raw-threshold rows
```

## Target Statement

Let:

```text
B_k = -Delta^2 log A_k = -log(x_k)
M_k = mu_(2k)
A_k = M_k*k!/(2*k)!
```

The missing estimate is:

```text
0 <= B_k <= 2/(3*(2*k+1))
```

for the actual negative-lambda zeta heat-flow coefficient tail from the
stated starting index.

The exact raw-moment translation is:

```text
x_k = ((2*k-1)/(2*k+1)) * (M_(k+1)*M_(k-1)/M_k^2)
log(M_(k+1)*M_(k-1)/M_k^2) >= log((2*k+1)/(2*k-1)) - 2/(3*(2*k+1))
```

Ordinary moment log-convexity gives only the nonnegative left side. The
threshold on the right is strictly positive, so that standard fact is
not enough.

## Live Routes

```text
1. uniform saddle/Laplace control of the actual heat-flow moment integrals
2. zeta-specific lower bounds on raw moment-ratio growth
```

Finite diagnostics:

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_22
checked contractions: x_1..x_21
raw log-convexity rows: 63 / 63
raw curvature-threshold rows: 63 / 63
simple log-buffer rows: 63 / 63
positive threshold rows: 63 / 63
```

Extrema:

```text
min raw curvature-threshold margin: 2.034879929997686858E-3 at lambda=-25.0, k=21
min simple log-buffer margin: 2.005614607560894716E-3 at lambda=-25.0, k=21
min plain-logconvexity shortfall: 3.212417302026237154E-2 at lambda=-100.0, k=21
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```

Summary:

The bounded log-curvature target is the missing quantitative half of the negative-lambda buffered tail: prove B_k=-Delta^2 log A_k<=2/(3*(2*k+1)) for k>=22. Equivalently, for raw moments M_k=mu_(2k), prove log(M_(k+1)*M_(k-1)/M_k^2) >= log((2*k+1)/(2*k-1))-2/(3*(2*k+1)). Plain moment log-convexity only gives the nonnegative left side and is too weak.

Historical correction:

```text
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md
validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows
```

The fixed wall `C_k=(2*k+1)*B_k<=2/3` is finite-rejected on the repaired
k300 data and should not be used as a live defect-tail route.
