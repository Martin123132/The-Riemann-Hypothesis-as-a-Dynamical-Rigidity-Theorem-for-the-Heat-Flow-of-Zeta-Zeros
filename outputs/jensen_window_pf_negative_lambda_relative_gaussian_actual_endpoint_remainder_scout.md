# Jensen-Window PF Negative-Lambda Relative-Gaussian Actual Endpoint Remainder Scout

Date: 2026-07-07

Status: floating endpoint theorem-search diagnostic. This is not a proof
of a uniform residual estimate, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout`.

Proof boundary: this artifact samples the actual relative-Gaussian
integral at the collar endpoint `T=1156`. It is not interval-certified
and it is not a theorem on the full collar.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout: 6 rows, 0 issues, 4 endpoint rows, 5 quadrature orders, 0 ready-to-apply rows
```

## Endpoint Setup

```text
T endpoint: 1156
u endpoint: 1/1156
quadrature orders: [96, 128, 192, 256, 320]
selected quadrature order: 320
polynomial degree: 40
normalizing Phi(0) float: 4.466969004671234589E-01
```

## Endpoint Residuals

```text
max selected value residual / first omitted term: 7.928285259094620674E-01
max selected derivative residual / first omitted term: 7.859721996680071321E-01
max value residual / first omitted term over orders: 7.990755683030266177E-01
max derivative residual / first omitted term over orders: 7.863346326661018182E-01
max selected value budget fraction: 5.249275170430453248E-04
max selected derivative budget fraction: 5.452217884786860185E-04
all selected residuals below first omitted term: True
```

Per-index selected-order rows:

```text
F_21: value residual=3.936102704571453614E-05, value first omitted=4.990210169953365261E-05, value ratio=7.887649158087939316E-01, derivative residual=7.125042884581489488E-07, derivative first omitted=9.065260689361647775E-07, derivative ratio=7.859721996680071321E-01
  spreads over orders: value=5.145232293557455705E-07, derivative=3.285549610154703259E-10; budget fractions=(7.312343864120115860E-05, 7.648708601762656246E-05)
F_22: value residual=7.820753086207332672E-05, value first omitted=9.864368940605490509E-05, value ratio=7.928285259094620674E-01, derivative residual=1.401352619723184034E-06, derivative first omitted=1.791970136269163607E-06, derivative ratio=7.820178424629132818E-01
  spreads over orders: value=1.715077431185818568E-07, derivative=2.628439688123762608E-10; budget fractions=(1.452910153393800773E-04, 1.504347133092229032E-04)
F_23: value residual=1.498120136140812519E-04, value first omitted=1.907111328517061011E-04, value ratio=7.855441440357477934E-01, derivative residual=2.695464900170918554E-06, derivative first omitted=3.464475596787048638E-06, derivative ratio=7.780296973864356813E-01
  spreads over orders: value=4.287693577964546421E-07, derivative=1.971329766092821956E-10; budget fractions=(2.783151357433080947E-04, 2.893572137271088881E-04)
F_24: value residual=2.825590067878636091E-04, value first omitted=3.611338473149329790E-04, value ratio=7.824218330370267260E-01, derivative residual=5.078934009361546487E-06, derivative first omitted=6.560389959873349845E-06, derivative ratio=7.741817240174541093E-01
  spreads over orders: value=1.715077431185818568E-07, derivative=1.971329766092821956E-10; budget fractions=(5.249275170430453248E-04, 5.452217884786860185E-04)
```

Scope warning:

```text
The endpoint quadrature is stable over N=96..320 at T=1156, but larger T values suffer double-precision cancellation after subtracting the degree-40 polynomial. This scout is therefore kept endpoint-local and is not promoted to a uniform theorem.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

At the collar endpoint T=1156, direct generalized Gauss-Laguerre quadrature of the actual relative-Gaussian multiplier gives value and derivative residuals for F_21..F_24 below one first omitted formal term and below 0.1% of the degree-40 half-safety budgets. This is useful evidence that the asymptotic-remainder target is numerically plausible, but it remains a floating endpoint scout rather than a uniform analytic or interval-certified remainder theorem.
