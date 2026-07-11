# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint X-Panel Route Matrix

Date: 2026-07-08

Status: endpoint x-panel route matrix. This is not a proof
of an x-panel certificate, an interpolation-remainder theorem, a compact interval
certificate, a finite-grid interval certificate, a uniform collar theorem, RH,
or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix`.

Proof boundary: this artifact records the endpoint x change of variables,
first-panel mass and Bernstein budgets, exact-moment obligations, and rejected
shortcuts. It does not prove the x-panel remainder.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix: 7 rows, 0 issues, x interval 0<=x<=0.01, 18 Bernstein budgets, 0 ready-to-apply rows
```

## Route Constants

```text
T: 10000
index: F_21
alpha: 20.5
y interval: 0<=y<=1
x interval: 0<=x<=0.01
x right endpoint: 0.01
x weight power: 42
first panel mass upper: 1.6155560239464103501984383230311864446832165667162E-21
heaviest panel: 20<=y<=50
heaviest panel mass upper: 0.60216159332489531950482867922076390781009134801886
first/heaviest mass ratio: 2.682927708866247294881825729652E-21
value unscaled cap: 6.782032247872604818E-40
derivative unscaled cap: 1.424226772053247012E-38
first-panel value sup budget without Bernstein: 6.996592450057070652146864508995E-20
first-panel derivative sup budget without Bernstein: 1.469284414511984837177801591129E-18
rho=2 degree=32 value sup budget: 7.512533939108907946133043813769E-11
rho=3 degree=32 value sup budget: 6.482413531562058989909765298732E-5
```

Endpoint x transform:

```text
y=T*x^2, dy=2*T*x dx, and the Gamma(alpha+1,1) density becomes 2*T^(alpha+1)*x^(2*alpha+1)*exp(-T*x^2)/Gamma(alpha+1) dx.
For monomial x^k on the first panel, the transformed Gamma moment equals T^(-k/2)*lower_gamma(alpha+1+k/2,1)/Gamma(alpha+1).
```

Selected first x-panel Bernstein budgets:

```text
N=16, rho=1.5: value M<=5.744559538809073532872278738852E-18; derivative M<=1.206357503149905442089524316534E-16
N=16, rho=2.0: value M<=1.146321707017350455647742281154E-15; derivative M<=2.407275584736435957232110126905E-14
N=16, rho=3.0: value M<=1.505901815741565772201270637718E-12; derivative M<=3.162393813057288122111162624334E-11
N=24, rho=1.5: value M<=1.472267778676809822233399250219E-16; derivative M<=3.091762335221300627167722281553E-15
N=24, rho=2.0: value M<=2.934583569964417166458220239754E-13; derivative M<=6.162625496925276050514201924877E-12
N=24, rho=3.0: value M<=9.880221813080413031412536654065E-9; derivative M<=2.074846580746886736917133797825E-7
N=32, rho=1.5: value M<=3.773261287460370798309895500268E-15; derivative M<=7.923848703666778677674775737996E-14
N=32, rho=2.0: value M<=7.512533939108907946133043813769E-11; derivative M<=1.577632127212870668931635692769E-9
N=32, rho=3.0: value M<=6.482413531562058989909765298732E-5; derivative M<=1.361306841628032388091331484753E-3
```

Rejected shortcut:

```text
Tiny first-panel Gamma mass is not a proof of the endpoint remainder.
```

Required upgrade:

```text
Certify the first-panel value and derivative cores in x, including an x-domain/sup-norm or Taylor-model remainder bound and exact transformed moments; then splice this with the y-panel route away from zero.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
```

Summary:

The endpoint x-panel route matrix quantifies the first-panel repair suggested by the endpoint parity matrix. The change y=T*x^2 maps 0<=y<=1 to 0<=x<=0.01 and gives transformed density power x^42 for alpha=20.5. The first-panel Gamma mass upper is 1.6155560239464103501984383230311864446832165667162E-21, only 2.682927708866247294881825729652E-21 of the heaviest panel mass. This makes the first panel numerically lightweight, but it remains a route matrix only: no x-domain, sup-norm, exact-moment implementation, or true x-panel interpolation-remainder certificate is proved.
