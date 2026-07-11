# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Interpolation-Remainder Route Matrix

Date: 2026-07-08

Status: worst-row interpolation-remainder route matrix. This is not a proof
of an interpolation-remainder theorem, compact interval-integration certificate,
finite-grid interval certificate, uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix`.

Proof boundary: this artifact records panel Gamma masses, Bernstein
sufficient-condition budgets, endpoint route constraints, and rejected
shortcuts. It does not prove an analytic-domain bound or a true
interpolation remainder.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row interpolation-remainder route matrix: 8 rows, 0 issues, 6 panel masses, 20 Bernstein budgets, 16 minimal-degree rows, 0 ready-to-apply rows
```

## Route Constants

```text
T: 10000
index: F_21
compact interval: 0<=y<=200
panel count: 6
heaviest panel: 20<=y<=50
heaviest panel mass upper: 0.60216159332489531950482867922076390781009134801886
value unscaled cap: 6.782032247872604818E-40
derivative unscaled cap: 1.424226772053247012E-38
rho=2 degree=128 value sup budget: 1.596890001110721912195980431156E-2
rho=2 degree=160 value sup budget: 6.858590330079954287832319494469E+7
target closing: False
```

Bernstein route formula:

```text
panel_error <= gamma_panel_mass * 4*rho^(-N)/(rho-1) * M
```

Selected heaviest-panel budgets:

```text
N=32, rho=2.0: value M<=2.015558546903338616076271270437E-31; derivative M<=4.232672948497011094413989738527E-30
N=64, rho=2.0: value M<=8.656758042123121429261484948149E-22; derivative M<=1.817919188845855500425725421185E-20
N=96, rho=2.0: value M<=3.718049268030379694411485528470E-12; derivative M<=7.807903462863797359470204761067E-11
N=128, rho=2.0: value M<=1.596890001110721912195980431156E-2; derivative M<=3.353469002332516016129568533521E-1
N=160, rho=2.0: value M<=6.858590330079954287832319494469E+7; derivative M<=1.440303969316790400667270535006E+9
```

Minimal degree rows:

```text
M=1e-12, rho=2.0: N=95
M=1e-12, rho=3.0: N=59
M=1e-6, rho=2.0: N=115
M=1e-6, rho=3.0: N=72
M=1, rho=2.0: N=134
M=1, rho=3.0: N=84
M=1e6, rho=2.0: N=154
M=1e6, rho=3.0: N=97
```

Rejected shortcut:

```text
Arb Cauchy deltas do not prove the interpolation remainder.
```

Required upgrade:

```text
For every panel, certify an analytic Bernstein-ellipse domain for the transformed cancellation-reduced value and derivative cores and an Arb-safe sup norm M on that domain.
The first panel touches y=0 while the finite core is evaluated through sqrt(y/T); a proof must either work in the x-variable near zero or certify the even analytic y-core after separating the tail/parity terms.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The interpolation-remainder route matrix isolates the remaining compact-interval theorem: the heaviest Gamma panel is 20<=y<=50 with mass upper 0.60216159332489531950482867922076390781009134801886, and the Bernstein route formula panel_error <= mass * 4*rho^(-N)/(rho-1) * M gives concrete degree/rho budgets. At rho=2 and degree 128 the heaviest-panel value sup-norm budget is 1.596890001110721912195980431156E-2, while degree 160 raises it to 6.858590330079954287832319494469E+7. This is still a route matrix only: no analytic-domain, sup-norm, endpoint, or true interpolation-remainder certificate is proved.
