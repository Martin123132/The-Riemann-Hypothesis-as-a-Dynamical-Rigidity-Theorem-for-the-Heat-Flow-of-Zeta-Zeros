# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Full-Expectation Certificate

Date: 2026-07-09

Status: worst-row full-expectation certificate. This is not a proof
of an all-row finite-grid certificate, a uniform collar theorem, RH,
or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate`.

Proof boundary: this artifact certifies the complete true relative-Gaussian
value and derivative expectations only for `T=10000`, `F_21`. It does not
certify the other finite-grid rows or bridge the finite grid to the collar.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation certificate: 8 rows, 0 issues, 3 composed sources, 2 below-one full ratios, 0 open worst-row integral sources, 0 ready-to-apply rows
```

## Global N-Tail

```text
For 0<=x<=1 use exp(a*x)<=exp(a) and exp(-pi*n^2*exp(4x))<=exp(-pi*n^2). For x>=1, every value monomial exp(a*x-pi*n^2*exp(4x)) and derivative-core monomial x*exp(a*x-pi*n^2*exp(4x)) is decreasing for n>=31 by the recorded margins.
tail start n: 31
global extension certified: True
value x>=1 monotonicity margin lower: 659333.6652360568650356637213688877114831986923290669462179724488414151
derivative x>=1 monotonicity margin lower: 659328.6652360568650356637213688877114831986923290669462179724488414151
value n-tail bound upper: 1.008683144279331099799696902420326427110667071269341986780776293780563E-1300
derivative n-tail bound upper: 6.650769108394031459514897397402897631443687085614223625134381545139400E-1295
Phi(0) n-tail bound upper: 1.245421034306393786675284658527977817836509415864526636582934445237231E-1304
```

## Normalization

```text
|Phi/Phi0-Phi30/Phi30_0| <= |Phi-Phi30|/Phi0_lower + |Phi30|*|Phi0-Phi30_0|/(Phi30_0_lower*Phi0_lower), integrated channelwise.
finite Phi_30(0) lower: 0.4466969004671234440869846670547091132204243670948249747308518817309554
full Phi(0) lower: 0.4466969004671234440869846670547091132204243670948249747308518817309554
both c0 lower bounds certified: True
value normalized n-tail correction upper: 1.106076054762397444658468363672754444199974889197454787063460429156833E-1297
derivative normalized n-tail correction upper: 7.449597682444430065414480853106939324504343307017351194459082911804034E-1295
```

## Complete Expectations

```text
finite entire value ball: [-6.58338692361019E-34 +/- 5.98E-49]
finite entire derivative ball: [-1.38058387415230E-32 +/- 3.30E-47]
full value expectation ball: [-6.58338692361019E-34 +/- 5.98E-49]
full derivative expectation ball: [-1.38058387415230E-32 +/- 3.30E-47]
full value certified negative: True
full derivative certified negative: True
value ratio / first omitted upper: 0.9707100590203297651200076500707905236731092924001689841866982789573976
derivative ratio / first omitted upper: 0.9693567774758049120828364647202533809659101015524784238953305230805692
value margin below one lower: 0.02928994097967023487999234992920947632689070759983101581330172104260244
derivative margin below one lower: 0.03064322252419508791716353527974661903408989844752157610466947691943075
both full ratios below one: True
complete worst-row expectation certified: True
quadrature needed for worst-row expectation: False
```

## Remaining Work

```text
Apply the same direct compact/global-tail certificate to the other recorded T/index rows.
Aggregate every finite-grid source and its rounding allowance.
Prove the finite-grid-to-full-collar and scaled-curvature handoff.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
```

Summary:

The worst-row full-expectation certificate composes three rigorous sources: the degree-128 finite-core compact x-moment integral, the finite-core y>=200 tail, and a new Arb global n>=31 normalization correction. For T=10000, F_21, the complete true value and derivative expectations are certified negative, with first-omitted ratio uppers below 0.971 and positive margins above 0.029. This removes generalized Gauss-Laguerre quadrature from the worst-row expectation proof. It remains a one-row certificate: the other finite-grid rows, all-source aggregation, and the finite-grid-to-collar theorem remain open.
