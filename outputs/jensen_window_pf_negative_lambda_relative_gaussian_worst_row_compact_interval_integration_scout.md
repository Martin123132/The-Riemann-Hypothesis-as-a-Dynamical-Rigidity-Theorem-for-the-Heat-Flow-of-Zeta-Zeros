# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Compact-Interval Integration Scout

Date: 2026-07-07

Status: worst-row compact-interval integration scout. This is not a proof
of a compact interval-integration certificate, quadrature-remainder
theorem, finite-grid interval certificate, uniform collar theorem, RH,
or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout`.

Proof boundary: this artifact diagnoses only raw Arb panel hulls
on the remaining compact interval `0<=y<=200` for `T=10000`,
`F_21`. Plain interval-Riemann hulls are explicitly rejected.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout: 7 rows, 0 issues, 6 panels, plain interval Riemann rejected, 0 ready-to-apply rows
```

## Compact Scout

```text
T: 10000
index: F_21
compact interval: 0<=y<=200
far-tail split y: 200
panel count: 6
value unscaled cap: 6.782032247872604818E-40
derivative unscaled cap: 1.424226772053247012E-38
value raw interval width bound: 0.9326768489129334646963769417351161026927
derivative raw interval width bound: 0.6094311540819450756312705442967424541098
value raw width to cap ratio upper: 1.375217360851530223059919868593402789065E+39
derivative raw width to cap ratio upper: 42790317247256500123688575487970253060.03
plain interval Riemann rejected: True
target closing: False
```

Required upgrade:

```text
Use local Taylor or Chebyshev panel models for the cancellation-reduced core, with Arb coefficient balls and exact Gamma-weighted panel moments.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

A six-panel raw Arb interval scout on the remaining compact interval 0<=y<=200 is far too wide: the value and derivative width-to-cap ratios are respectively 1.375217360851530223059919868593402789065E+39 and 42790317247256500123688575487970253060.03. Plain interval-Riemann hulls are therefore rejected; the live route is a local Taylor/Chebyshev panel model with exact Gamma-weighted moments. This does not close the compact interval.
