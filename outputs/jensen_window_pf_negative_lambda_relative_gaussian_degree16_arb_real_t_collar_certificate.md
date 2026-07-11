# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Arb Real-T Collar Certificate

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof
of an infinite Taylor-tail theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate`.

Proof boundary: this artifact certifies a real-`T` collar only for the
degree-16 finite surrogate at fixed `k=22`, using Arb coefficient-ratio
balls. It does not bound the infinite residual Taylor tail.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified stencil rows, 0 ready-to-apply rows
```

## Arb Finite-Degree Collar

```text
real T interval: [1156, infinity)
u interval: [0, 1/1156]
ratio ball rows: 9
positive normalizer rows: 4
certified stencil rows: 3
failed Bernstein rows: 0
Bernstein subdivisions: 1
finite-degree Arb collar certified: True
```

Key ratio balls:

```text
c12/c0: [-3304359.71281471153631522353332 +/- 4.76e-24] (radius [1.674061151e-65 +/- 1.96e-75], sign negative)
c14/c0: [30122498.1078104844171418342933 +/- 4.05e-23] (radius [9.034890952e-64 +/- 4.91e-74], sign positive)
c16/c0: [-237753170.949415650656506482065 +/- 4.12e-22] (radius [4.114298378e-62 +/- 4.67e-74], sign negative)
```

Bernstein certificates:

```text
F_21: degree=8, positive Bernstein coefficients=9/9, certified=True
F_22: degree=8, positive Bernstein coefficients=9/9, certified=True
F_23: degree=8, positive Bernstein coefficients=9/9, certified=True
F_24: degree=8, positive Bernstein coefficients=9/9, certified=True
B_product: degree=16, zero order=2, Bernstein degree=14, positive coefficients=15/15, certified=True
companion_product: degree=32, zero order=3, Bernstein degree=29, positive coefficients=30/30, certified=True
weighted_gap_derivative: degree=31, zero order=1, Bernstein degree=30, positive coefficients=31/31, certified=True
```

Interpretation:

The rationalized real-`T` collar has been upgraded to an Arb
coefficient-ball finite-degree collar. The remaining proof gap is now
sharper: prove signed infinite-tail stencil bounds beyond degree 16 on
the same real collar, and then remove the fixed-`k` limitation.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The degree-16 real-T collar now survives Arb coefficient-ratio balls for the finite surrogate: Bernstein coefficients certify F_21..F_24, the stripped B product, the stripped companion product, and the stripped weighted-gap derivative numerator on 0<=u<=1/1156. This upgrades the midpoint surrogate to an interval finite-degree surrogate, while leaving the infinite residual Taylor-tail theorem open.
