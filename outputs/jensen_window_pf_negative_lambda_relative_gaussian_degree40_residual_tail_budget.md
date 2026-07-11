# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-40 Residual Tail Budget

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of an infinite Taylor-tail theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget`.

Proof boundary: this artifact converts finite degree-40 collar margins
into sufficient value and derivative residual-tail targets. It does not
prove the residual estimates.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: 8 rows, 0 issues, 5 budget inequalities, 4 finite tail profile rows, 0 ready-to-apply rows
```

## Residual Target

Write:

```text
F_i(u) = P_i^(40)(u) + R_i(u),  i in {21,22,23,24}
0 <= u <= 1/1156
```

A sufficient fixed-`k=22` residual target is:

```text
|R_i(u)|  <= 5.382819486765314521E-01 * u^3
|R_i'(u)| <= 9.315354075509573936E-03 * u
```

The value-tail budget is limited by the companion product.

## Margin Data

```text
normalizer_min_lower: [0.439627934927551258358085649181 +/- 2.62e-31]
B_product_min_lower: [37.2033572903651204156102209610 +/- 2.68e-29]
companion_product_min_lower: [17.2250223756550133023672457518 +/- 2.11e-29]
weighted_gap_derivative_min_lower: [27.4244024269680055723947152492 +/- 3.52e-29]
finite_normalizer_abs_upper: [1.00000000000000000000000000000 +/- 1e-34]
finite_derivative_abs_upper: [917.618341545662737252398795656 +/- 4.88e-28]
```

Budget inequalities:

```text
normalizer_value_residual: allocated=2.198139674637756280E-01, raw threshold=3.395695876365208626E+08, half-safety=5.382819486765314521E-01, bound at half-safety=3.484466661937167012E-10, limiting=False
B_product_value_residual: allocated=1.860167864518255953E+01, raw threshold=5.375875774509763687E+03, half-safety=5.382819486765314521E-01, bound at half-safety=1.862567296783089075E-03, limiting=False
companion_product_value_residual: allocated=8.612511187827506021E+00, raw threshold=1.076563897353062904E+00, half-safety=5.382819486765314521E-01, bound at half-safety=4.306255591663002313E+00, limiting=True
weighted_gap_value_residual: allocated=6.856100606742001169E+00, raw threshold=9.044023628251490976E+00, half-safety=5.382819486765314521E-01, bound at half-safety=4.080612044257607152E-01, limiting=False
weighted_gap_derivative_residual: allocated=6.856100606742001169E+00, raw threshold=1.863070815101914787E-02, half-safety=9.315354075509573936E-03, bound at half-safety=3.428050303371000584E+00, limiting=True
```

## Finite Tail Profile

```text
F_21, degrees 42..80: value fraction=[0.000125542567495511743 +/- 1.93e-23], derivative fraction=[0.000134001677480446364 +/- 4.16e-23], largest value degree=42, largest derivative degree=42
F_22, degrees 42..80: value fraction=[0.000250228988828346650 +/- 3.53e-22], derivative fraction=[0.000267227352867346908 +/- 9.09e-23], largest value degree=42, largest derivative degree=42
F_23, degrees 42..80: value fraction=[0.000487828627931362480 +/- 3.38e-22], derivative fraction=[0.000521239897545509047 +/- 4.75e-22], largest value degree=42, largest derivative degree=42
F_24, degrees 42..80: value fraction=[0.000931556974617909435 +/- 2.12e-22], derivative fraction=[0.000995886692207194084 +/- 4.49e-22], largest value degree=42, largest derivative degree=42
```

The finite profile is a plausibility diagnostic only. A proof still
needs an analytic majorant for every term beyond degree 40, or a
stronger argument that controls the residual as a whole.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The degree-40 Arb collar margins imply a concrete sufficient fixed-k residual target on 0<=u<=1/1156: prove |R_i(u)|<=A*u^3 with half-safety A=5.382819486765314521E-01 and |R_i'(u)|<=B*u with half-safety B=9.315354075509573936E-03 for i=21..24. The finite degree-42..80 tail profile uses less than 0.1% of these budgets, but this remains a target because no analytic majorant for the infinite residual tail is proved.
