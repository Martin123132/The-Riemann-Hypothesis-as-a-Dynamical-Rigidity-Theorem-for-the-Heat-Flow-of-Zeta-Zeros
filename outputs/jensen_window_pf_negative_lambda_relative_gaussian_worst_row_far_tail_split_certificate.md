# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Far-Tail Split Certificate

Date: 2026-07-07

Status: worst-row far-tail split certificate. This is not a proof
of a compact interval-integration certificate, finite-grid interval
certificate, uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate`.

Proof boundary: this artifact certifies only the finite `n<=30`
cancellation-reduced continuum far tail `y>=200` for `T=10000`,
`F_21`. It does not integrate the compact interval `0<=y<=200`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: 7 rows, 0 issues, split y=200, 2 tail ratios below quadrature cap, 0 ready-to-apply rows
```

## Tail Bounds

```text
T: 10000
index: F_21
split y: 200
split x: 0.1414213562373095048801688724209698078570
tail mass upper: 2.061872387924679053576321193215761896360E-59
value total tail bound: 5.874503628626422470963622648546258502383E-59
derivative total tail bound: 6.405843031697552640324141439308114278820E-59
value scaled tail bound: 5.874503628626422470963622648546258502383E-47
derivative scaled tail bound: 6.405843031697552640324141439308114278820E-51
value ratio to first omitted upper: 8.661863308699458186697067463495530095599E-26
derivative ratio to first omitted upper: 4.497769005186247039689660658162936457902E-27
quadrature ratio radius cap: 0.0000010
value tail ratio below quadrature cap: True
derivative tail ratio below quadrature cap: True
remaining compact interval: 0<=y<=200
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst-row interval-integration route now has an Arb-certified far-tail split: on y>=200, finite n<=30 Phi majorants are monotone decreasing, the polynomial part is bounded by upper incomplete-Gamma moments, and the combined value and derivative tails are far below the quadrature route matrix caps. This retires only the far-tail part; compact integration on 0<=y<=200, aggregation, all-row coverage, and grid-to-collar coverage remain open.
