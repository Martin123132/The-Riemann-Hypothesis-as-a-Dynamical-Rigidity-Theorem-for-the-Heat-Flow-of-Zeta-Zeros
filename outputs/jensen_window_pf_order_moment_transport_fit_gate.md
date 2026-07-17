# Jensen-Window PF Order-Moment Transport Fit Gate

Date: 2026-07-12

Status: primary-source theorem-fit rejection gate. This is not a proof
of order-four entry, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_order_moment_transport_fit_gate`.

```text
work/rh_compute/results/jensen_window_pf_order_moment_transport_fit_gate.json
python work/rh_compute/scripts/jensen_window_pf_order_moment_transport_fit_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_order_moment_transport_fit_gate.py
```

## Primary Theorem

[Order-Moment Transport and Hankel Determinants in Special-Function Inequalities](https://arxiv.org/abs/2606.31647v2) records a positive order-moment
transport criterion: a Gamma-normalized average of a completely monotone
function becomes a positive Mellin moment family and hence a totally
nonnegative Hankel kernel. The paper explicitly treats reciprocal
sign-regularity as a separate problem.

## Exact Fit Test

For the first Newman summand at `lambda=-100`, put `y=100*u^2`. Then

```text
phi_1(u)=pi*exp(5*u)*(2*q-3)*exp(-q), q=pi*exp(4*u)
A_t^(1)=sqrt(pi)/(10*400^t*Gamma(t+1/2))*integral_0^infinity y^(t-1/2)*exp(-y)*g(y)dy
g(y)=phi_1(sqrt(y)/10)
```

This is formally the right Gamma-average shape. However, the positive
transport theorem requires `g` to be completely monotone. At the origin

```text
d_u log(phi_1)(0)=5+8*pi/(2*pi-3)-4*pi,
  = [0.0886165024991730927364637276247675413757025436842714931604458 +/- 2.11E-63],
  > 0.0886165024991730927364637276247675413757025436842714931604458 > 0.
```

Thus `phi_1` increases immediately to the right of zero, and because
`y -> sqrt(y)/10` is increasing, `g'(y)>0` for sufficiently small
positive `y`. The transformed kernel is not completely monotone.

## Consequence

The criterion would in any event produce ordinary positive Hankel moment
orientation. The Xi programme needs the alternating reciprocal-Gamma
sign orientation. Therefore this theorem cannot be promoted into the
missing bridge.

The useful handoff is precise: seek a signed reciprocal-Gamma transport
theorem, or continue the direct stable-gap curvature certificate. This
gate rejects one theorem application, not total-positivity methods as a
whole.

```text
outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/jensen_window_pf_theorem_machinery_fit_matrix.md
```
