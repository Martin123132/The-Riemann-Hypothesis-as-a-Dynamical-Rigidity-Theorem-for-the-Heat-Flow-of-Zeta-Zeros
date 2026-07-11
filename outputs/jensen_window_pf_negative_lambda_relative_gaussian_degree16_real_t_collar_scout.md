# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Real-T Collar Scout

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof
of an infinite Taylor-tail theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout`.

Proof boundary: this artifact certifies a real-`T` collar only for the
rationalized degree-16 finite surrogate at fixed `k=22`. It does not
prove the corresponding analytic zeta-tail theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified surrogate stencil rows, 0 ready-to-apply rows
```

## Real-T Surrogate Collar

```text
real T interval: [1156, infinity)
u interval: (0, 1/1156]
positive normalizer rows: 4
certified product stencil rows: 2
certified derivative stencil rows: 1
certified surrogate stencil rows: 3
root-count failures: 0
real-T surrogate collar certified: True
```

Endpoint log-stencils at `T=1156`:

```text
B: 1.259416498980593755E-4
companion: 8.851471227548859301E-10
weighted gap: 2.518416978813492713E-4
```

Normalizer root counts:

```text
F_21: degree=8, roots=0, endpoint sign=positive, certified=True
F_22: degree=8, roots=0, endpoint sign=positive, certified=True
F_23: degree=8, roots=0, endpoint sign=positive, certified=True
F_24: degree=8, roots=0, endpoint sign=positive, certified=True
```

Stencil root-count certificates:

```text
B: degree=16, zero order=2, stripped degree=14, roots=0, certified=True
companion: degree=32, zero order=3, stripped degree=29, roots=0, certified=True
weighted_gap derivative: degree=30, zero order=1, stripped degree=29, roots=0, certified=True
```

Interpretation:

The integer scan threshold `T=1156` is no longer just a sampled
integer threshold for the degree-16 surrogate: within the rationalized
finite model, all three structured stencil signs persist on the full
real half-line `T>=1156`. The live proof route is now to replace this
surrogate statement with interval coefficient control and signed
infinite-tail stencil bounds.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

For the rationalized degree-16 finite surrogate at k=22, the integer collar threshold T=1156 extends to the whole real half-line T>=1156: the four normalizers are positive, the B and companion product stencils have no stripped roots on 0<u<=1/1156, and the weighted-gap derivative numerator is positive there. This is a finite-surrogate real-T collar, not an infinite Taylor-tail theorem.
