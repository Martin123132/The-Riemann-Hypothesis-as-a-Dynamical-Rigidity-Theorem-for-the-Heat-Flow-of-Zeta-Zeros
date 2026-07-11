# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Laguerre Root-Bracket Certificate

Date: 2026-07-07

Status: worst-row Laguerre root-bracket certificate. This is not a proof
of a Christoffel-weight interval certificate, a quadrature-remainder
theorem, a finite-grid interval certificate, a uniform collar theorem,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate`.

Proof boundary: this artifact certifies individual root brackets for
`L_320^(41/2)` at the worst recorded row `T=10000`, `F_21`. It does
not certify weights, Phi/Phi' node values, quadrature error, rounding
aggregation, all recorded rows/orders, or grid-to-collar coverage.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: 6 rows, 0 issues, 320 root brackets, 30 zero floating weights, 2 intervalization rows, 0 ready-to-apply rows
```

## Root Brackets

```text
T: 10000
index: F_21
quadrature order: 320
alpha: 41/2
root brackets: 320
bisection steps: 60
widest bracket: 4.145787935554796563E-17 at root 320
narrowest bracket: 1.710100983788340277E-19 at root 2
degree-count argument: The degree-320 polynomial has 320 disjoint sign-changing brackets, so each bracket contains exactly one real root and there are no additional roots.
```

Sample brackets:

```text
root 1: [0.5093662041285411437875587192911190637101270795028540305793285369873046875, 0.509366204128541144314004159229826440569155465709627605974674224853515625], signs=(positive,negative), width=0.0000000000000000005264454399387073768590283862067735753953456878662109375
root 2: [0.7045343332865549646222140603339066655674827188704512082040309906005859375, 0.70453433328655496479322415871274069322505173573745196335949003696441650390625], signs=(negative,positive), width=0.00000000000000000017101009837883402765756901686700075515545904636383056640625
root 319: [1252.2151918412503443865977181163801935925850017383709200657904148101806640625, 1252.21519184125034441014756583938154731183711732001029304228723049163818359375], signs=(positive,negative), width=0.00000000000000002354984772300135371925211558163937297649681568145751953125
root 320: [1282.1894468699346277823258131839975586252489136995791341178119182586669921875, 1282.189446869934627823783692539545524258226549818573403172194957733154296875], signs=(negative,positive), width=0.0000000000000000414578793555479656329776361189942690543830394744873046875
```

## Weight Boundary

```text
floating weight count: 320
zero floating weights: 30
first zero floating weight index: 291
last positive floating weight index: 290
last positive floating weight: 9.303e-321
first zero floating weight node: 886.9174722103623
```

This underflow is a diagnostic boundary, not a theorem. The next
certificate needs non-floating Christoffel-weight intervals on these
root brackets before the quadrature row can become intervalized.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst relative-Gaussian ladder row now has certified Arb brackets for all 320 nodes of L_320^(41/2). The brackets are sign-changing, ordered, and disjoint, with the widest width 4.145787935554796563E-17 at root 320. The same row exposes the remaining weight problem: SciPy double weights underflow to zero for 30 tail nodes, so the next proof step must certify Christoffel weights non-float.
