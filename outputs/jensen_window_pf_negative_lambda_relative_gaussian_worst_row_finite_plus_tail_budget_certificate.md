# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Finite-Plus-Tail Budget Certificate

Date: 2026-07-07

Status: worst-row finite-plus-tail budget certificate. This is not a proof
of a quadrature-remainder theorem, a finite-grid interval certificate,
a uniform collar theorem, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate`.

Proof boundary: this artifact composes only the worst-row `T=10000`,
`F_21`, `N=320` finite `n<=30` weighted-sum certificate with the
reserved full `n>30` Phi-tail per-source budget. It does not prove
quadrature remainder, cover all rows, aggregate rounding, or bridge
the grid to a full collar.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate: 6 rows, 0 issues, 2 composed ratios, 3 tail sources, 0 ready-to-apply rows
```

## Budget Composition

```text
per-source intervalization cap: 2.000000000000000000E-3
tail budget ratio reserved: 0.002000000000000000000
tail source certified: True
tail actual max ratio to cap: 3.778846084314790602E-1292
value finite ratio upper: 0.9833957992836557769419015895036210773888
derivative finite ratio upper: 0.9694055674762067320093698741711260875260
value composed ratio upper using full tail cap: 0.9853957992836557769419015895036210773888
derivative composed ratio upper using full tail cap: 0.9714055674762067320093698741711260875260
value remaining margin after full tail cap: 0.0146042007163442230580984104963789226112
derivative remaining margin after full tail cap: 0.0285944325237932679906301258288739124740
both composed ratios below one: True
```

Remaining sources:

```text
generalized Gauss-Laguerre quadrature-remainder or interval adaptive integration beyond the N=320 sum
rounding and cross-source aggregation under the total intervalization cap
all recorded rows and quadrature orders, not just T=10000, F_21, N=320
finite-grid to full-collar coverage
final acceptance-gate wording separating finite-grid evidence from theorem claims
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst-row finite-plus-tail budget certificate composes the certified finite n<=30 weighted-sum ratio uppers with the full 2.0e-3 Phi-tail source cap. The resulting value and derivative ratio uppers remain strictly below one first omitted term. This retires only the worst-row finite-plus-tail numerator-source budget; quadrature remainder, rounding aggregation, all-row coverage, and finite-grid to collar coverage remain open.
