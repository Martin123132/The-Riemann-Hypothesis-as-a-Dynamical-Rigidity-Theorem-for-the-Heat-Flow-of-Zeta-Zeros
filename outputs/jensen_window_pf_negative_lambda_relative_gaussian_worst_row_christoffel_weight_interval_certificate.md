# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Christoffel-Weight Interval Certificate

Date: 2026-07-07

Status: worst-row Christoffel-weight interval certificate. This is not a proof
of Phi/Phi' node interval evaluation, a quadrature-remainder theorem, a
finite-grid interval certificate, a uniform collar theorem, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate`.

Proof boundary: this artifact certifies positive Christoffel-weight
intervals for the certified `L_320^(41/2)` root brackets in the single
worst row `T=10000`, `F_21`, `N=320`. It does not certify Phi/Phi'
node values, quadrature error, all recorded rows/orders, or
grid-to-collar coverage.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate: 6 rows, 0 issues, 320 interval weights, 0 Taylor denominator obstructions, 30 repaired floating underflows, 0 ready-to-apply rows
```

## Interval Weights

```text
formula: w_j=Gamma(N+alpha+1)/(Gamma(N+1)*(N+1)^2)*x_j/[L_(N+1)^(alpha)(x_j)]^2
Taylor identity: d^m/dx^m L_n^(alpha)(x)=(-1)^m L_(n-m)^(alpha+m)(x)
precision bits: 1024
quadrature order: 320
index: F_21
interval weights: 320
Taylor denominator obstructions: 0
direct recurrence denominator obstructions: 320
repaired floating underflows: 30
minimum weight interval: [2.5663925426991066e-492 +/- 6.09e-509] at root 320
maximum weight interval: [7.762348957842933e+17 +/- 46.2] at root 43
maximum relative weight width: 6.645705603047399351898944020607455890895e-16 at root 1
```

Sample intervals:

```text
root 1: weight=[1.179850762112515e-7 +/- 6.52e-23], denominator=[-4.317764415572101e+26 +/- 8.24e+10], SciPy float=1.1798507621144197e-07, underflow=False
root 2: weight=[7.358642403872687e-5 +/- 8.74e-21], denominator=[2.033338340728291e+25 +/- 1.30e+9], SciPy float=7.35864240387532e-05, underflow=False
root 43: weight=[7.762348957842933e+17 +/- 46.2], denominator=[-1074357788129241.98 +/- 8.45e-3], SciPy float=7.762348957843034e+17, underflow=False
root 100: weight=[3.1701448567708983 +/- 5.25e-17], denominator=[1.12014408851914076e+24 +/- 4.46e+6], SciPy float=3.1701448567710773, underflow=False
root 290: weight=[9.3031561782133030e-321 +/- 6.60e-338], denominator=[6.38402301363671726e+184 +/- 9.27e+166], SciPy float=9.303e-321, underflow=False
root 291: weight=[1.58785672523795015e-324 +/- 6.65e-342], denominator=[-4.9112649867334720e+186 +/- 4.06e+169], SciPy float=0.0, underflow=True
root 320: weight=[2.5663925426991066e-492 +/- 6.09e-509], denominator=[4.6448530094451616e+270 +/- 3.29e+253], SciPy float=0.0, underflow=True
```

## Mass Check

```text
weight sum interval: [1.1082798113786904e+19 +/- 2.58e+2]
mass Gamma(43/2) contained: True
relative weight-sum mass error interval: [+/- 1.26e-17]
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst-row Christoffel-weight interval certificate replaces the midpoint-only obstruction with a Taylor-centered Arb enclosure of L_321^(41/2) on each certified L_320^(41/2) root bracket. All 320 denominator intervals avoid zero, all 320 Christoffel-weight intervals are positive, and the interval weight sum contains Gamma(43/2). This certifies the worst-row weights only; Phi/Phi' node evaluation, quadrature remainder, all-row coverage, rounding aggregation, and the grid-to-collar bridge remain open.
