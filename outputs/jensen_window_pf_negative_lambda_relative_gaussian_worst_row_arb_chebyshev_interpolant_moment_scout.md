# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Arb Chebyshev Interpolant-Moment Scout

Date: 2026-07-08

Status: worst-row Arb Chebyshev interpolant-moment scout. This is not a proof
of a compact interval-integration certificate, quadrature-remainder
theorem, finite-grid interval certificate, uniform collar theorem, RH,
or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout`.

Proof boundary: this artifact certifies Arb arithmetic for Chebyshev
interpolants and incomplete-Gamma panel moments. It does not bound the
interpolation remainder between the interpolant and the true core.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout: 7 rows, 0 issues, 4 degrees, 3 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows
```

## Arb Interpolant Ladder

```text
T: 10000
index: F_21
compact interval: 0<=y<=200
panel count: 6
degrees: [16, 20, 24, 32]
value unscaled cap: 6.782032247872604818E-40
derivative unscaled cap: 1.424226772053247012E-38
first cap-safe pair: 16->20
cap-safe pair count: 3
reference degree: 32
reference value interpolant integral ball: [-6.58338692361019331200869759429797241582343365917043174997833E-34 +/- 1.63E-94]
reference derivative interpolant integral ball: [-1.38058387415230156402674205260518855676282541898259838947598E-32 +/- 3.91E-92]
target closing: False
```

Cauchy rows:

```text
16->20: value delta/cap 0.01084142765529810925122002325367656958771; derivative delta/cap 0.01075165326468561580445350311071629542923; below caps True
20->24: value delta/cap 6.969243783765862010324414119702129783777E-10; derivative delta/cap 6.723320659192506578602117131047695988453E-10; below caps True
24->32: value delta/cap 1.136682952723511405736824962015391755683E-18; derivative delta/cap 1.294620035297571747757230219488783361619E-18; below caps True
```

Required upgrade:

```text
Complete this route by bounding the difference between each Arb Chebyshev interpolant and the true cancellation-reduced core on its panel.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The Arb Chebyshev interpolant-moment route upgrades the floating scout's arithmetic: degree 16,20,24,32 interpolant integrals are Arb-enclosed, all consecutive Cauchy deltas are below the unscaled quadrature caps, and the degree-32 interpolant balls are [-6.58338692361019331200869759429797241582343365917043174997833E-34 +/- 1.63E-94] and [-1.38058387415230156402674205260518855676282541898259838947598E-32 +/- 3.91E-92]. This still does not close the compact interval because panel interpolation remainders remain unbounded.
