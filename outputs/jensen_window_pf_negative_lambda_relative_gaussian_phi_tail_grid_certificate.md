# Jensen-Window PF Negative-Lambda Relative-Gaussian Phi-Tail Grid Certificate

Date: 2026-07-07

Status: finite-grid Phi-tail source certificate. This is not a proof
of a finite-grid interval certificate, quadrature-remainder theorem,
uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate`.

Proof boundary: this artifact certifies only the omitted `n>30`
Phi-tail source for the recorded finite grid by combining the padded
tail majorants with certified `x<=1` and `Phi(0)>=0.44` side
conditions. It does not certify finite `n<=30` node evaluations.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: 6 rows, 0 issues, 3 certified tail sources, 2 certified side conditions, 0 ready-to-apply rows
```

## Certified Inputs

```text
padded tail range: 0<=x<=1
node range x<=1 certified: True
certified c0 lower: 0.44
c0 lower certified by n=1 term: True
per-source intervalization cap: 2.000000000000000000E-3
finite-grid tail source certified: True
```

Tail source ratios to cap:

```text
Phi value n>30 tail: normalized/relative=2.2924616915439343177265838691371055161606069801575954245017643040467343671748619e-1300, ratio_to_cap=1.146230845771967159E-1297, certified=True
Phi prime n>30 derivative-core tail: normalized/relative=7.5576921686295812039942015879578382175496444154707086649254335740220458593283757e-1295, ratio_to_cap=3.778846084314790602E-1292, certified=True
Phi(0) n>30 denominator tail: normalized/relative=2.8305023506963495151711014966544950405375213996921059922339419209937073753762793e-1304, ratio_to_cap=1.415251175348174758E-1301, certified=True
```

Remaining obligations:

```text
individual Laguerre node and weight intervals beyond the coarse x<=1 range
finite n<=30 interval evaluation of Phi and Phi' at certified node intervals
coefficient-ratio ball propagation for the cancellation-reduced polynomial core
generalized Gauss-Laguerre quadrature-remainder or interval adaptive integration
rounding aggregation and finite-grid to full-collar coverage
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The Phi-tail source for the recorded finite grid is now certificate-composed rather than merely conditional: the padded n>30 majorants are far below the 2.0e-3 per-source cap, while the node-c0 certificate proves x<=1 and Phi(0)>=0.44 for the same grid. This certifies the omitted-n tail component only; finite n<=30 node evaluation, weights, quadrature error, coefficient propagation, rounding, and grid-to-collar coverage remain open.
