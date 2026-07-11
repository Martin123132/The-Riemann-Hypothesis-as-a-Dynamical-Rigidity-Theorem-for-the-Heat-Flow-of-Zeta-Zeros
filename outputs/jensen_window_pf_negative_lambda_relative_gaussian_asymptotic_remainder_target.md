# Jensen-Window PF Negative-Lambda Relative-Gaussian Asymptotic Remainder Target

Date: 2026-07-07

Status: exact theorem-search diagnostic. This is not a proof
of an analytic residual estimate, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target`.

Proof boundary: this artifact converts the degree-40 residual budget
and the formal-tail obstruction into explicit sufficient asymptotic
remainder targets. It does not prove those targets.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: 6 rows, 0 issues, 4 first-omitted rows, 4 optimized-window rows, 0 ready-to-apply rows
```

## First-Omitted-Term Target

```text
common multiplier limit: [1419.93907869608 +/- 1.36e-12]
safe common multiplier tested: 1000
safe value budget fraction: [0.67090090649120 +/- 5.01e-15]
safe derivative budget fraction: [0.704255566315065 +/- 6.49e-16]
safe multiplier closes target: True
```

Per-index first-omitted rows:

```text
F_21: value term=[4.990210169953362e-5 +/- 8.08e-21], derivative term=[9.06526068936164e-7 +/- 4.60e-22], common multiplier limit=[10275.88107470690 +/- 8.47e-12], 1000x fractions=([0.0927062514769953 +/- 8.02e-17], [0.0973152562519826 +/- 7.86e-17]), limiting=derivative
F_22: value term=[9.86436894060548e-5 +/- 5.09e-20], derivative term=[1.791970136269162e-6 +/- 3.40e-22], common multiplier limit=[5198.38689661643 +/- 2.53e-12], 1000x fractions=([0.183256543617316 +/- 4.76e-16], [0.192367367009733 +/- 1.14e-16]), limiting=derivative
F_23: value term=[0.0001907111328517060 +/- 4.72e-20], derivative term=[3.464475596787047e-6 +/- 4.81e-22], common multiplier limit=[2688.82080859471 +/- 5.95e-12], 1000x fractions=([0.354295984326812 +/- 5.90e-16], [0.371910242885484 +/- 4.18e-16]), limiting=derivative
F_24: value term=[0.000361133847314933 +/- 4.01e-19], derivative term=[6.560389959873344e-6 +/- 9.35e-22], common multiplier limit=[1419.939078696080 +/- 9.63e-13], 1000x fractions=([0.670900906491196 +/- 8.35e-16], [0.704255566315065 +/- 4.56e-16]), limiting=derivative
```

## Optimized-Window Target

This conditional route is different from summing the infinite formal
tail: it would need a theorem justifying a finite formal window and
an actual remainder bound after truncation near the least term.

```text
formal window: 21..120
least term j: 103
max value window budget fraction: [0.00093155697477913 +/- 2.65e-18]
max derivative window budget fraction: [0.00099588669254223 +/- 1.83e-18]
common least-term multiplier limit after window: [4.86406709155425e+29 +/- 6.04e+14]
```

Per-index optimized-window rows:

```text
F_21: window fractions=([0.000125542567502489 +/- 3.57e-19], [0.000134001677494933 +/- 4.24e-19]), least terms j=(103, 103), common least-term multiplier=[8.4642570199773e+31 +/- 3.37e+17]
F_22: window fractions=([0.000250228988848804 +/- 2.94e-19], [0.000267227352909834 +/- 4.71e-19]), least terms j=(103, 103), common least-term multiplier=[1.46150423923556e+31 +/- 5.65e+16]
F_23: window fractions=([0.000487828627989596 +/- 5.39e-19], [0.000521239897666488 +/- 7.31e-19]), least terms j=(103, 103), common least-term multiplier=[2.61956098454486e+30 +/- 7.04e+15]
F_24: window fractions=([0.00093155697477913 +/- 2.51e-18], [0.00099588669254223 +/- 1.41e-18]), least terms j=(103, 103), common least-term multiplier=[4.86406709155425e+29 +/- 4.63e+14]
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The degree-40 residual budget can be rephrased as a concrete asymptotic-remainder target: a common 1000x first-omitted-term theorem would still fit inside the half-safety value and derivative budgets for F_21..F_24 on 0<=u<=1/1156. A finite optimized-window route through j=120 leaves still larger least-term multiplier slack, but both routes require a real analytic remainder theorem rather than an infinite formal-tail sum.
