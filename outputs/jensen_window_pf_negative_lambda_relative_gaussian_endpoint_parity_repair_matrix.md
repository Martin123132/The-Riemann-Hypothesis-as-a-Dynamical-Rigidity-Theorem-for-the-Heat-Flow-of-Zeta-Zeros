# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint Parity-Repair Matrix

Date: 2026-07-08

Status: endpoint parity-repair route matrix. This is not a proof
of endpoint analyticity, an interpolation-remainder theorem, a finite-grid
interval certificate, a uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix`.

Proof boundary: this artifact audits low-order odd Taylor coefficients
of the finite Phi truncation and records admissible endpoint repair
routes. It does not prove exact evenness or an endpoint certificate.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix: 7 rows, 0 issues, 8 odd Taylor rows, order 15, 0 ready-to-apply rows
```

## Endpoint Audit

```text
T: 10000
compact interval: 0<=y<=200
finite Phi terms: 30
series order: 15
finite c0 lower: 0.44669690046712344408698466705470911322042436709482
first odd coefficient abs upper: 1.5013918057273630212637957520005850517148059822056E-1300
max low-order odd degree: 15
max low-order odd abs upper: 1.5394049948039318674157069296860766261840406690307E-1255
first-panel x max: 0.010000000000000000000000000000000000000000000000000
compact x max: 0.14142135623730950488016887242096980785696718753769
first-panel normalized odd partial sum upper: 3.4968885647053311252358566122466130159385287853700E-1285
compact normalized odd partial sum upper: 6.2387347102785735263519242935701452522435065420123E-1268
```

Odd Taylor rows:

```text
x^1: abs upper 1.5013918057273630212637957520005850517148059822056E-1300; normalized 3.3610974335333771338009276426908465616419385949398E-1300
x^3: abs upper 3.6402333177422417605157097846505247017713474376942E-1293; normalized 8.1492244829448065419853034852974827427812665567137E-1293
x^5: abs upper 2.6442866456651337708055774371741139265073336453880E-1286; normalized 5.9196440425262168998525180553132723967562665310977E-1286
x^7: abs upper 9.1346236222453862753105577054654767907510573019252E-1280; normalized 2.0449265738564684094432368531592424927767138336319E-1279
x^9: abs upper 1.8382768416486128924309833265473645449222806531567E-1273; normalized 4.1152666152961334134122383885311720286863057552644E-1273
x^11: abs upper 2.4181962937565474806170908414451027274599109597727E-1267; normalized 5.4135058721647093880132054889305968790988989433733E-1267
x^13: abs upper 2.2400621261547380224121562980549042786895815434925E-1261; normalized 5.0147250267737304574005660913282057052890360282893E-1261
x^15: abs upper 1.5394049948039318674157069296860766261840406690307E-1255; normalized 3.4461958280752183487438861717159644648360127954165E-1255
```

Accepted repairs:

```text
Use the theta functional equation to split the infinite Phi kernel into an exact even analytic core in x, hence an analytic core in y=x^2, and charge the finite truncation/tail parity defect to the existing Phi-tail source.
Alternatively handle the first panel in x with y=T*x^2 and Gamma density 2*T^(alpha+1)*x^(2*alpha+1)*exp(-T*x^2)/Gamma(alpha+1), avoiding a y-plane branch theorem.
```

Rejected shortcut:

```text
Near-evenness and low-order cancellation are not an endpoint analyticity proof.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
```

Summary:

The endpoint parity-repair matrix makes the first-panel branch obligation explicit. For the finite n<=30 Phi truncation, Arb series coefficients through odd degree 15 are tail-sized rather than exactly zero: the first odd coefficient has absolute upper 1.5013918057273630212637957520005850517148059822056E-1300, and the largest recorded low-order odd coefficient is degree 15 with upper 1.5394049948039318674157069296860766261840406690307E-1255. This supports, but does not prove, the correct repair: either use exact evenness of the infinite theta kernel plus a certified tail charge, or handle the endpoint panel in the x variable.
