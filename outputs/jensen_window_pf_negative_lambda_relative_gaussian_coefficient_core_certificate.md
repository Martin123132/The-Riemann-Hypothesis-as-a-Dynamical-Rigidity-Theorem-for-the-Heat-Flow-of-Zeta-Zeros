# Jensen-Window PF Negative-Lambda Relative-Gaussian Coefficient-Core Propagation Certificate

Date: 2026-07-07

Status: coefficient-core propagation certificate. This is not a proof
of the first-omitted residual theorem, a finite-grid interval certificate,
a uniform collar theorem, scaled-curvature monotonicity, cone entry, RH,
or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate`.

Proof boundary: this artifact propagates Arb coefficient-ratio ball
uncertainty through exact Gamma moments for the recorded finite grid.
It does not enclose finite Phi node values, Laguerre nodes or weights,
quadrature error, rounding aggregation, or grid-to-collar coverage.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: 6 rows, 0 issues, 22 coefficient rows, 20 propagation rows, 2 intervalization rows, 0 ready-to-apply rows
```

## Coefficient Source

```text
polynomial M: 20
polynomial degree: 40
first omitted j: 21
coefficient rows: 22
maximum ratio-ball radius upper: 1.1634227377488650115891040471172799978436518705721939601658222013621873399754518e-79 at r_21
source: build_ratio_rows with finite n<=80 sums plus geometric Arb tail-radius bound
```

## Radius Envelopes

```text
value core radius on v<=809/1156: 2.8695439668152282534239556934615637137779356727158943664836385899231389027150510e-84
derivative core radius on v<=809/1156: 5.7270309011434727992727339456970634466222332504294291633377035961663182322989527e-83
value core radius on v<=1: 3.5690383692270561790794609418488291706489557694643835585133897337876638112409829e-81
derivative core radius on v<=1: 7.1276951231372052472421477733323668731252511254856793833907408012840068378592106e-80
```

## Scaled Grid Handoff

```text
sum_{j=0}^{20} rad(r_j)*(i+1/2)_j*u^(j-3)
sum_{j=1}^{20} j*rad(r_j)*(i+1/2)_j*u^(j-2)
maximum value coefficient ratio: 8.610446518945492184E-81 at T=10000, F_21
maximum derivative coefficient ratio: 5.523600431279558041E-83 at T=10000, F_21
ratio cap: 1.000000000000000000E-6
intervalization per-source cap: 2.000000000000000000E-3
all coefficient ratios below ratio cap: True
```

Worst-row detail:

```text
T=10000, F_21: value ratio=8.6104465189454921839562073903461162535847674039760459699213927184581294613348575e-81, derivative ratio=5.5236004312795580405642469868102621863953562447551247810676863262230377459745125e-83
T=10000, F_21: value ratio=8.6104465189454921839562073903461162535847674039760459699213927184581294613348575e-81, derivative ratio=5.5236004312795580405642469868102621863953562447551247810676863262230377459745125e-83
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The Arb coefficient-ratio balls r_0..r_21 propagate through the exact Gamma-moment scaling with negligible first-omitted impact on the recorded finite grid: the worst value coefficient-radius ratio is at T=10000, F_21 and the worst derivative coefficient-radius ratio is also at T=10000, F_21, both far below the 1e-6 ratio-radius target. This retires the coefficient-ball source for the finite-grid intervalization ledger, but not Phi node evaluation, quadrature, rounding, or grid-to-collar coverage.
