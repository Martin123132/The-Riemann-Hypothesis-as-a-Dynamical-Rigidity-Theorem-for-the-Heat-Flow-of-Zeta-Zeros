# Jensen-Window PF Negative-Lambda Relative-Gaussian Quadrature Ladder Scout

Date: 2026-07-07

Status: high-order floating quadrature ladder scout. This is not a proof
of a quadrature-remainder theorem, a finite-grid interval certificate,
a uniform residual estimate, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout`.

Proof boundary: this artifact tests the worst cancellation-reduced
finite-grid row at higher floating quadrature orders. It calibrates
a future rigorous quadrature-radius target but does not prove it.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: 5 rows, 0 issues, 7 ladder rows, 320 reference order, 0 ready-to-apply rows
```

## Ladder Summary

```text
ladder orders: [96, 128, 160, 192, 224, 256, 320]
reference order: 320
max value ratio: 0.9707100590203351233111029
max derivative ratio: 0.9693567774758094418100653
max value ratio spread: 7.767821744285352977404633e-15
max derivative ratio spread: 7.756680887393827019918227e-15
all ladder ratios below one: True
proposed quadrature ratio radius cap: 1.0e-6
intervalization per-source cap: 2.000000000000000000E-3
cap below per-source cap: True
cap keeps worst ladder below one: True
```

Worst-row ladder:

```text
T=10000, F_21
N=96: value=0.9707100590203295412978665, derivative=0.9693567774758038687433225, delta_value=-7.645824343684963955788635e-16, delta_derivative=-7.62993120142187004368876e-16
N=128: value=0.9707100590203290642671517, derivative=0.9693567774758033914379359, delta_value=-1.241613149238753592653189e-15, delta_derivative=-1.240298506779145668073936e-15
N=160: value=0.9707100590203273554893586, derivative=0.9693567774758016851291779, delta_value=-2.950390942308912476588076e-15, delta_derivative=-2.946607264801217696074909e-15
N=192: value=0.9707100590203277038344049, derivative=0.9693567774758020339809452, delta_value=-2.602045896053889653546352e-15, delta_derivative=-2.597755497438248870549006e-15
N=224: value=0.9707100590203351233111029, derivative=0.9693567774758094418100653, delta_value=4.817430801976440500816557e-15, delta_derivative=4.810073622592609323843318e-15
N=256: value=0.9707100590203305103566682, derivative=0.9693567774758048357957197, delta_value=2.044763672796231384763709e-16, delta_derivative=2.04059276980833138920819e-16
N=320: value=0.9707100590203303058803009, derivative=0.9693567774758046317364427, delta_value=0.0, delta_derivative=0.0
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst cancellation-reduced row remains below one first omitted term through the high-order floating ladder N=96..320. The largest value ratio is about 0.970710059020335 and the largest derivative ratio is about 0.969356777475809, with order-spread below 1e-14. A future rigorous quadrature radius below 1e-6 would be far inside the 2.0e-3 per-source intervalization cap, but this artifact remains floating evidence only.
