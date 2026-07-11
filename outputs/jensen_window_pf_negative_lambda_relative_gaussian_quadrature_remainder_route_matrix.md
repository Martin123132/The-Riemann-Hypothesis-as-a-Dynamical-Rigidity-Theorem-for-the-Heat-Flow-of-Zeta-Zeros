# Jensen-Window PF Negative-Lambda Relative-Gaussian Quadrature-Remainder Route Matrix

Date: 2026-07-07

Status: quadrature-remainder route matrix. This is not a proof
of a quadrature-remainder theorem, a finite-grid interval certificate,
a uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix`.

Proof boundary: this artifact records the exact classical
Gauss-Laguerre remainder prefactor and the derivative/interval
obligations needed to turn the floating quadrature ladder into a
certificate. It does not prove those obligations.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: 7 rows, 0 issues, derivative order 640, 2 derivative-sup caps, 0 ready-to-apply rows
```

## Remainder Constants

```text
T: 10000
index: F_21
quadrature order: 320
derivative order: 640
Laguerre error prefactor: 2.791244149880588131514918942024363250922057058755595558776957E-159
quadrature ratio radius cap: 0.0000010
intervalization per-source cap: 2.000000000000000000E-3
quadrature cap below per-source cap: True
value unscaled expectation error cap: 6.782032247872604818E-40
derivative unscaled expectation error cap: 1.424226772053247012E-38
value required 640th derivative sup cap: 2.4297524271256407631984588349963763207485E+119
derivative required 640th derivative sup cap: 5.1024800969638456035049425453208614326121E+120
```

## Budget Compatibility

```text
finite-plus-tail value ratio upper: 0.9853957992836557769419015895036210773888
finite-plus-tail derivative ratio upper: 0.9714055674762067320093698741711260875260
value ratio upper after quadrature cap: 0.9853967992836557769419015895036210773888
derivative ratio upper after quadrature cap: 0.9714065674762067320093698741711260875260
both ratios still below one after quadrature cap: True
```

Live routes:

```text
Cauchy derivative route: prove 640th y-derivative bounds for the cancellation-reduced cores.
Interval adaptive route: integrate the cancellation-reduced Gamma expectation with a separate far-tail bound.
Rejected route: do not promote floating N=96..320 order-spread into a remainder theorem.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The classical N=320 generalized Gauss-Laguerre remainder formula reduces the worst-row quadrature source to 640th-derivative control with normalized prefactor about 2.79e-159. A 1e-6 ratio-radius certificate would follow from value and derivative 640th-derivative supremum bounds below the recorded caps, and adding that cap to the finite-plus-tail budget still leaves both channels below one. This is a route matrix only: no derivative or interval integration bound is proved.
