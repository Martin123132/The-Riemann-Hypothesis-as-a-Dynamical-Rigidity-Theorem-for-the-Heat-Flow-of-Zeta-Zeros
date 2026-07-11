# Jensen-Window PF Negative-Lambda Relative-Gaussian Cancellation-Reduced Remainder Grid Scout

Date: 2026-07-07

Status: floating cancellation-reduced theorem-search diagnostic. This is not a proof
of a uniform residual estimate, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout`.

Proof boundary: the residual is subtracted inside the Gamma expectation
and evaluated with `mpmath`, but the Laguerre nodes and weights are
floating SciPy values. This is not interval-certified.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: 6 rows, 0 issues, 20 grid rows, 5 T values, 0 ready-to-apply rows
```

## Cancellation-Reduced Setup

```text
T grid: [1156, 1500, 2000, 5000, 10000]
indices: [21, 22, 23, 24]
quadrature orders: [64, 96, 128, 192]
selected quadrature order: 192
polynomial degree: 40
mpmath dps: 80
Phi n terms: 30
normalizing Phi(0): 0.446696900467123444086984667055
```

The cancellation-reduced cores are:

```text
Phi(sqrt(v))/Phi(0)-sum_{j=0}^{20} r_j*v^j
sqrt(v)*Phi'(sqrt(v))/(2*Phi(0))-sum_{j=1}^{20} j*r_j*v^j
|E[value_core(uY)]|/u^3
|E[derivative_core(uY)]|/u^2
```

## Grid Result

```text
grid rows: 20
max value ratio / first omitted: 0.970710059020329541 at {'T': 10000, 'index': 21}
max derivative ratio / first omitted: 0.969356777475803869 at {'T': 10000, 'index': 21}
max value ratio spread over orders: 9.63500159925082018e-15
max derivative ratio spread over orders: 9.62202755450426673e-15
max selected value budget fraction: 0.000524800202052571167
max selected derivative budget fraction: 0.000545203459703802311
all grid ratios below one: True
```

Selected-order rows:

```text
T=1156, F_21: value ratio=0.793632465394960028, derivative ratio=0.785867221301462357, value budget fraction=0.0000735746909172129669, derivative budget fraction=0.0000764768700209853109
T=1156, F_22: value ratio=0.789795521254646766, derivative ratio=0.781924180209349361, value budget fraction=0.000144735197389563294, derivative budget fraction=0.000150416695748116516
T=1156, F_23: value ratio=0.785995563508640381, derivative ratio=0.778020532715779366, value budget fraction=0.000278475071849800744, derivative budget fraction=0.000289353805292219004
T=1156, F_24: value ratio=0.782232065831107238, derivative ratio=0.774155698273732754, value budget fraction=0.000524800202052571167, derivative budget fraction=0.000545203459703802311
T=1500, F_21: value ratio=0.832919449999716672, derivative ratio=0.826314924294672874, value budget fraction=7.10085059796274152e-7, derivative budget fraction=5.69890546524686566e-7
T=1500, F_22: value ratio=0.8296560372984043, derivative ratio=0.82294890969633497, value budget fraction=1.39815692953385305e-6, derivative budget fraction=1.12193888513778676e-6
T=1500, F_23: value ratio=0.826418113628439029, derivative ratio=0.819610190634395384, value budget fraction=2.69255391384505958e-6, derivative budget fraction=2.16028183913626596e-6
T=1500, F_24: value ratio=0.823205385923381321, derivative ratio=0.816298441562352474, value budget fraction=5.0788446898323392e-6, derivative budget fraction=4.07421723079900983e-6
T=2000, F_21: value ratio=0.869131337699485177, derivative ratio=0.86372870600152942, value budget fraction=4.17729861341652421e-9, derivative budget fraction=2.51876247320243147e-9
T=2000, F_22: value ratio=0.86646180203696939, derivative ratio=0.86096589018337621, value budget fraction=8.23208799927314396e-9, derivative budget fraction=4.9630228712944418e-9
T=2000, F_23: value ratio=0.86380860321109603, derivative ratio=0.858220656184892816, value budget fraction=1.58666355642320752e-8, derivative budget fraction=9.5645828343107917e-9
T=2000, F_24: value ratio=0.861171594502097685, derivative ratio=0.855492839977983097, value budget fraction=2.99536096989327354e-8, derivative budget fraction=1.80540897602240188e-8
T=5000, F_21: value ratio=0.943112359075079744, derivative ratio=0.940559976278100585, value budget fraction=3.11496658753248021e-16, derivative budget fraction=7.53938966690198406e-17
T=5000, F_22: value ratio=0.941851181693046795, derivative ratio=0.939245687824777413, value budget fraction=6.149257982854366e-16, derivative budget fraction=1.48826194099532929e-16
T=5000, F_23: value ratio=0.940593357066984882, derivative ratio=0.937935043935528828, value budget fraction=1.18726884768808249e-15, derivative budget fraction=2.87329136321343015e-16
T=5000, F_24: value ratio=0.939338872218428673, derivative ratio=0.936628029915317562, value budget fraction=2.24523399410952134e-15, derivative budget fraction=5.43333151068675714e-16
T=10000, F_21: value ratio=0.970710059020327704, derivative ratio=0.969356777475802034, value budget fraction=1.22303691212322666e-21, derivative budget fraction=1.48205195740428974e-22
T=10000, F_22: value ratio=0.970041378727742466, derivative ratio=0.968658149750003292, value budget fraction=2.41596570356468307e-21, derivative budget fraction=2.92752616324411605e-22
T=10000, F_23: value ratio=0.969373613014350293, derivative ratio=0.967960519859675693, value budget fraction=4.6676516538975674e-21, derivative budget fraction=5.65580765382364008e-22
T=10000, F_24: value ratio=0.968706760064949248, derivative ratio=0.967263885740108982, value budget fraction=8.832664258757127e-21, derivative budget fraction=1.07022257812921509e-21
```

Scope warning:

```text
The value and derivative remainders are subtracted inside the Gamma expectation. This removes the far-T double-precision cancellation seen when F_i and P_i^(40) are evaluated separately, but the quadrature is still floating-node diagnostic evidence, not an interval enclosure.
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

Subtracting the degree-40 polynomial inside the Gamma expectation removes the far-T catastrophic cancellation from separate floating subtraction. On the grid T=1156,1500,2000,5000,10000 and F_21..F_24, all sampled value and derivative residuals are below one first omitted formal term; the worst sampled value ratio is about 0.9707100590 at T=10000, F_21, and the worst sampled derivative ratio is about 0.9693567775 at the same row. This is strong diagnostic support for the first-omitted remainder target, but it is still finite floating evidence rather than an interval-certified or analytic uniform collar theorem.
