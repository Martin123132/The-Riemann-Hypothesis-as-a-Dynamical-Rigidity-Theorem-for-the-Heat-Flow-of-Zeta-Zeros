# Jensen-Window PF Suzuki Jordan-Totient Sign Scout

Date: 2026-07-23

Status: finite double-precision reconnaissance report. It is not a proof
of eventual sign, the all-time determinant gate, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_suzuki_jordan_totient_sign_scout.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_jordan_totient_sign_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py
```

## Target

Suzuki's scalar criterion uses

```text
c_omega(n)=n^omega*product_(p|n)(1-p^(-2*omega))>0,
h_omega^<1>(x)=x^(-1)*sum_(n<=x)c_omega(n)*g_omega^<1>(n/x).
```

Eventual one-sign behavior of `h_omega^<1>` implies innerness of
`Theta_omega`. On a sequence `omega->0`, that is RH-strength.

## Finite Grid

The deterministic grid has 600 linear points on `[1.001,20]` and
600 geometric points on `[20.01,5000]`. Reported values are
`sqrt(x)*h_omega^<1>(x)`.

| omega | samples | negative | minimum | x at minimum | maximum | x at maximum |
|---:|---:|---:|---:|---:|---:|---:|
| 1/2 | 1200 | 0 | 0.089316135 | 1.001000 | 1.125984150 | 5.378065 |
| 1/4 | 1200 | 0 | 0.310359981 | 1.001000 | 1.095393242 | 5.251194 |
| 1/8 | 1200 | 0 | 0.563176163 | 1.001000 | 1.076343111 | 113.184729 |
| 1/16 | 1200 | 0 | 0.752655758 | 1.001000 | 1.051977141 | 113.184729 |
| 1/32 | 1200 | 0 | 0.868221382 | 1.001000 | 1.030326261 | 113.184729 |

All 6,000 sampled values are positive. The smallest sampled value
is separated from zero by more than `0.08`; every sampled arithmetic
coefficient is positive.

Anchor values:

| omega | x=10 | x=100 | x=1000 | x=5000 |
|---:|---:|---:|---:|---:|
| 1/2 | 0.930274556 | 0.934014845 | 1.004921234 | 0.993340880 |
| 1/4 | 0.957454905 | 0.917619524 | 1.001660070 | 0.991804699 |
| 1/8 | 0.976644377 | 0.934266973 | 0.997260955 | 0.994465957 |
| 1/16 | 0.987779045 | 0.958533571 | 0.996673727 | 0.997556372 |
| 1/32 | 0.993750915 | 0.976727740 | 0.997641016 | 0.999122386 |

## Signed-Weight Guard

The positivity is not termwise. At `omega=1/2`, the explicit
primitive weight has opposite signs at two probes:

```text
g_(1/2)^<1>(exp(-4))=-39.803349588<0
g_(1/2)^<1>(1/2)=1.174060048>0
```

Thus the observed positive summatory values arise after a
genuinely signed arithmetic convolution; positivity of
`c_omega(n)` alone is not a proof.

## Interpretation

The grid makes the scalar target empirically credible and supplies
well-separated regression anchors. It gives no control beyond
`x=5000`, no uniformity as `omega->0`, and no rigorous rounding
certificate. Finite positivity cannot be promoted to eventual sign.

## Source

- Masatoshi Suzuki, `A canonical system of differential equations arising from the Riemann zeta-function`: https://arxiv.org/abs/1204.1827

## Proof Boundary

This finite grid is numerical reconnaissance only. It does not prove positivity or eventual one-sign behavior on an unbounded interval, the L2 residual, meromorphic innerness, any all-time determinant condition, RH, or Lambda<=0.
