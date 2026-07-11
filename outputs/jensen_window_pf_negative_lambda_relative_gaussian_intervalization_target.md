# Jensen-Window PF Negative-Lambda Relative-Gaussian Intervalization Target

Date: 2026-07-07

Status: open numerical-certification target. This is not a proof
of a uniform residual estimate, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target`.

Proof boundary: this artifact translates the cancellation-reduced
floating grid into explicit interval-certification obligations and
ratio-error budgets. It does not prove those obligations.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: 6 rows, 0 issues, 8 obligations, 5 open requirements, 0 ready-to-apply rows
```

## Error Budget

```text
observed worst value ratio: 0.970710059020329541
observed worst derivative ratio: 0.969356777475803869
common slack to one: 2.928994097967045900E-2
proposed total ratio error cap: 1.000000000000000000E-2
proposed per-source cap for five sources: 2.000000000000000000E-3
closed if total cap met: True
```

Per-channel budget rows:

```text
value: observed=9.707100590203295410E-1 at T=10000, F_21, slack=2.928994097967045900E-2, half_slack=1.464497048983522950E-2, cap=1.000000000000000000E-2, closes=True
derivative: observed=9.693567774758038690E-1 at T=10000, F_21, slack=3.064322252419613100E-2, half_slack=1.532161126209806550E-2, cap=1.000000000000000000E-2, closes=True
```

## Certification Obligations

```text
nlrgit_01_residual_core_identity: exact_reduction / available_exact / cap=None / Use the cancellation-reduced Gamma-expectation residual cores for value and derivative channels.
nlrgit_02_laguerre_node_weight_intervals: open_requirement / not_ready_to_apply / cap=2.000000000000000000E-3 / Replace SciPy floating generalized Laguerre nodes and weights by interval enclosures for every recorded index and order.
nlrgit_03_phi_and_c0_interval_tail: open_requirement / not_ready_to_apply / cap=2.000000000000000000E-3 / Enclose Phi, Phi', Phi(0), and their n>30 tails on all node-induced v ranges used by the finite grid.
nlrgit_04_quadrature_remainder_error: open_requirement / not_ready_to_apply / cap=2.000000000000000000E-3 / Bound the generalized Gauss-Laguerre quadrature remainder for the cancellation-reduced cores, not just the spread across orders.
nlrgit_05_ratio_and_coefficient_ball_propagation: open_requirement / not_ready_to_apply / cap=2.000000000000000000E-3 / Propagate Arb coefficient-ratio balls and first-omitted-term denominators through the ratio comparisons.
nlrgit_06_rounding_and_aggregation_budget: open_requirement / not_ready_to_apply / cap=1.000000000000000000E-2 / Account for arithmetic rounding, summation, and cross-source aggregation so the total ratio error remains below the proposed cap.
nlrgit_07_finite_grid_not_uniform_collar: rejected_route / not_ready_to_apply / cap=None / An intervalized finite T grid would by itself prove the full collar residual theorem on 0<=u<=1/1156.
nlrgit_08_acceptance_gate: acceptance_gate / not_ready_to_apply / cap=1.000000000000000000E-2 / A promoted proof must state all interval sources, cap their total ratio error, and separately bridge from finite grid certification to the full collar or all-k theorem.
```

Open gap:

This target is a certification roadmap for the finite cancellation-reduced grid. It does not provide node/weight intervals, quadrature-remainder bounds, Phi n-tail bounds, or a continuum-in-T theorem.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The cancellation-reduced finite grid leaves a common ratio slack of about 2.928994097967e-2 before reaching one first omitted term. A future interval certificate with total ratio error below 1.0e-2 in both value and derivative channels would keep the finite grid below the first-omitted threshold. The remaining work is to intervalize Laguerre nodes/weights, Phi tails, coefficient propagation, quadrature error, and then separately bridge the finite grid to the full collar.
