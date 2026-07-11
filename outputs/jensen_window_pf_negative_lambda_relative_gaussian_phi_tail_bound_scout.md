# Jensen-Window PF Negative-Lambda Relative-Gaussian Phi Tail Bound Scout

Date: 2026-07-07

Status: analytic padded-range tail scout. This is not a proof
of a finite-grid interval certificate, a uniform residual estimate,
scaled-curvature monotonicity, cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout`.

Proof boundary: this artifact bounds the omitted `n>30` Phi tails on
the padded range `0<=x<=1`, conditionally on later interval proofs
of the node range and `Phi(0)>=0.44`. It does not prove those
interval facts or any quadrature-remainder theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: 6 rows, 0 issues, 3 tail bounds below 1e-1000, 2 conditional requirements, 0 ready-to-apply rows
```

## Node Range

```text
max observed node x: 8.224376518739603403E-01
padded x range: 1.000000000000000000E+00
observed slack to padded range: 1.775623481260396597E-01
all observed nodes inside padded range: True
```

## Tail Bounds

```text
Phi truncation n: 30
tail start n: 31
value Phi tail bound: 1.0086831442793310997996969024203264271106670712693419867807762937805631215569392e-1300
derivative Phi-prime tail bound: 6.6507691083940314595148973974028976314436870856142236251343815451394003562089706e-1295
c0 tail bound: 1.2454210343063937866752846585279778178365094158645266365829344452372312451655629e-1304
normalized value tail bound using c0 proxy: 2.2924616915439343177265838691371055161606069801575954245017643040467343671748619e-1300
normalized derivative-core tail bound using c0 proxy: 7.5576921686295812039942015879578382175496444154707086649254335740220458593283757e-1295
denominator relative tail bound using c0 proxy: 2.8305023506963495151711014966544950405375213996921059922339419209937073753762793e-1304
per-source intervalization cap: 2.000000000000000000E-3
tail bounds below per-source cap: True
```

Conditional requirements:

```text
Replace the floating SciPy node range by interval enclosures proving x<=1 for all grid nodes.
Replace the floating c0 lower proxy by an interval-certified lower bound Phi(0)>=0.44.
```

Remaining non-tail obligations:

```text
Laguerre node/weight intervals beyond the x-range check.
Generalized Gauss-Laguerre quadrature-remainder bound.
Coefficient-ratio and first-omitted denominator interval propagation.
Rounding aggregation below the total ratio cap.
Continuum-in-T or grid-to-collar bridge.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

On the padded range 0<=x<=1, the omitted n>30 tails in Phi, Phi', and Phi(0) are far below the 2.0e-3 per-source intervalization ratio cap after normalization by the diagnostic c0 lower proxy 0.44. This narrows the nlrgit_03 tail source to two concrete certification tasks: interval-prove the grid node range x<=1 and interval-prove Phi(0)>=0.44.
