# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Christoffel-Weight Midpoint Scout

Date: 2026-07-07

Status: worst-row Christoffel-weight midpoint scout. This is not a proof
of Christoffel-weight intervals, a quadrature-remainder theorem, a
finite-grid interval certificate, a uniform collar theorem, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout`.

Proof boundary: this artifact evaluates the Christoffel formula at
Arb midpoints of the certified `L_320^(41/2)` root brackets. It repairs
floating tail-weight underflow diagnostically, but it does not certify
weights on the whole root brackets.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout: 6 rows, 0 issues, 320 midpoint weights, 30 repaired floating underflows, 320 direct interval obstructions, 0 ready-to-apply rows
```

## Midpoint Weights

```text
formula: w_j=Gamma(N+alpha+1)/(Gamma(N+1)*(N+1)^2)*x_j/[L_(N+1)^(alpha)(x_j)]^2
quadrature order: 320
index: F_21
midpoint weights: 320
repaired floating underflows: 30
minimum midpoint weight: 2.566392542699106565813978687580860385146e-492 at root 320
maximum midpoint weight: 776234895784293262.1491715698567429748912 at root 43
relative weight-sum mass error: 1.795129455585788358215558840375043474388e-18
```

Sample weights:

```text
root 1: midpoint weight=[1.179850762112514740051681486057186567279e-7 +/- 5.07e-48], SciPy float=1.1798507621144197e-07, underflow=False
root 2: midpoint weight=[7.358642403872686701124085852697694628209e-5 +/- 3.40e-46], SciPy float=7.35864240387532e-05, underflow=False
root 100: midpoint weight=[3.170144856770898335923451568515402245142 +/- 1.89e-40], SciPy float=3.1701448567710773, underflow=False
root 290: midpoint weight=[9.303156178213303039839888749661506853929e-321 +/- 5.93e-362], SciPy float=9.303e-321, underflow=False
root 291: midpoint weight=[1.587856725237950152170645753206264752246e-324 +/- 3.82e-364], SciPy float=0.0, underflow=True
root 300: midpoint weight=[1.899032545773128525896029534306442000241e-361 +/- 3.34e-401], SciPy float=0.0, underflow=True
root 320: midpoint weight=[2.566392542699106565813978687580860385146e-492 +/- 2.79e-532], SciPy float=0.0, underflow=True
```

## Interval Boundary

```text
direct interval denominator contains zero rows: 320
current root brackets are too wide for naive interval recurrence
next step: certified weight intervals via a sharper recurrence, Taylor/Lipschitz enclosure, or refined roots
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The worst-row Christoffel-weight midpoint scout repairs the N=320 double-weight underflow without using floating arithmetic for the weights: all 320 Arb midpoint weights are positive, including 30 weights that SciPy reports as zero. The smallest midpoint weight is at root 320, and the midpoint-weight sum matches Gamma(43/2) to relative error about 1.8e-18. This is still not a weight interval certificate because direct interval evaluation over the current root brackets contains zero in the Christoffel denominator on all 320 rows.
