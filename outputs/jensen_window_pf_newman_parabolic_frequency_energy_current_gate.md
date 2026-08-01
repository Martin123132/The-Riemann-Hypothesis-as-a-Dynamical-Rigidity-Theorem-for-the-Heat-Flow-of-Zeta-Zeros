# Parabolic-Frequency Energy/Current Gate

Date: 2026-07-28

Status: exact multiplier identity with positive bulk and an
explicit pointwise nonpromotion guard. This is not a proof of
Q209, the descendant theorem, Lambda<=0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_parabolic_frequency_energy_current_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_parabolic_frequency_energy_current_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_parabolic_frequency_energy_current_gate.py
```

Current result:

```text
validated Newman parabolic-frequency energy/current gate: 18 rows, 8 exact identities, 3 coefficient bounds, 2 pointwise bridge rows, 3 nonpromotion guards, 1 open Xi bulk-margin target, 0 pointwise contact exclusions
```

## Exact Local Balance

H is real and sufficiently smooth, H_t=-H_xx, and s(t,x)>0 is C1 on a fixed spatial interval.

Set

```text
E_s=H^2+s^2 H_x^2
J_s=H H_x+s^2 H_x H_xx
```

Direct differentiation gives

```text
partial_t E_s=-2H H_xx+2s s_t H_x^2-2s^2 H_x H_xxx
partial_x J_s=H_x^2+H H_xx+2s s_x H_x H_xx+s^2 H_xx^2+s^2 H_x H_xxx
```

and exact square completion gives

```text
partial_t E_s+2 partial_x J_s=2(s H_xx+s_x H_x)^2+2(1+s s_t-s_x^2)H_x^2
```

No term involving `s_t` or `s_x` has been discarded.

## Parabolic-Frequency Coefficient

```text
L=log(x/(4*pi)), q=2tL^2, s_pf=sqrt(2t)/sqrt(1+q)
partial_t_log_s_pf=1/[2t(1+q)]
partial_x_log_s_pf=-2tL/[x(1+q)]
s_pf_times_s_pf_t=1/(1+q)^2
s_pf_x_squared=4q t^2/[x^2(1+q)^3]
A_pf=1+[1-r(t,x)]/(1+q)^2, r(t,x)=4q t^2/[x^2(1+q)]
```

On every ray-aligned outer collar, `x>=38` and
`0<t<=1/4`. Therefore

```text
0<=r<=(4t^2/x^2)<=1/5776
A_pf>=1+5775/[5776(1+q)^2]>1
```

The requested two-chart audit is

```text
q<=1 implies A_pf>=28879/23104>1
q>=1 implies A_pf>1, but the certified surplus above 1 tends to 0 as q tends to infinity
```

Thus the bulk is genuinely coercive in both charts. The
frequency chart has no claimed uniform surplus above one as
`q` tends to infinity.

Pi provenance:

pi is the ordinary circle constant already present in the established Xi coordinate L=log(x/(4*pi)); it is not introduced by a polygon, curvature image, or chosen circle.

## Integrated Identities

For a fixed finite interval `[a,b]`,

```text
d/dt integral_a^b E_pf dx+2[J_pf(t,b)-J_pf(t,a)]=2 integral_a^b B_pf dx
```

and on a time slab,

```text
integral_a^b(E_pf(t_+,x)-E_pf(t_-,x))dx+2 integral_t_-^t_+[J_pf(t,b)-J_pf(t,a)]dt=2 integral_t_-^t_+ integral_a^b B_pf dxdt
```

If J_pf(t,b)=J_pf(t,a), then the spatially integrated energy is nondecreasing in Newman time.

This is integrated derivative control, not a pointwise lower
bound for `E_pf`.

## Exact Pointwise Bridge

The bulk has the additional exact interpretation

```text
V_pf=(H,s_pf H_x), E_pf=|V_pf|^2
partial_x V_pf=(H_x,s_pf H_xx+s_pf,x H_x)
B_pf=|partial_x V_pf|^2+(A_pf-1)H_x^2>=|partial_x V_pf|^2
```

Consequently, a contact at any unknown `x_0` forces

```text
If V_pf(t,x_0)=0 for some x_0 in [a,b], then integral_a^b B_pf dx>=[sqrt(E_pf(t,a))+sqrt(E_pf(t,b))]^2/(b-a)
```

At a contact x_0, Cauchy-Schwarz on [a,x_0] and [x_0,b] gives integral |V_x|^2>=E(a)/(x_0-a)+E(b)/(b-x_0)>=[sqrt(E(a))+sqrt(E(b))]^2/(b-a), with endpoint cases read by continuity.

Combining this floor with the fixed-slice balance gives the
strict conditional criterion

```text
If both endpoint jets are certified and d/dt integral_a^b E_pf dx+2[J_pf(t,b)-J_pf(t,a)]<2[sqrt(E_pf(t,a))+sqrt(E_pf(t,b))]^2/(b-a), then V_pf(t,x) has no zero on [a,b].
```

This is a real pointwise bridge, but its Xi-specific strict
upper bound and endpoint margins have not been proved.

## Contact Cancellation

At `H=H_x=0`,

```text
E_pf=0, partial_t E_pf=0, J_pf=0, partial_x J_pf=s_pf^2 H_xx^2
2 partial_x J_pf=2s_pf^2 H_xx^2; the positive bulk is exactly canceled by the flux divergence
```

There is no contradiction: the positive curvature term is
carried entirely by the local flux divergence.

The exact shifted quadratic model makes this obstruction
concrete:

```text
H(t,x)=(x-x_0)^2-2(t-t_0)
H_t=-H_xx
H(t_0,x_0)=H_x(t_0,x_0)=0 and H_xx=2
E=E_t=J=0, J_x=4s_0^2, 2J_x=8s_0^2=positive bulk
```

A genuine backward-heat double contact obeys the coercive balance, so that balance alone cannot exclude a first-jet collision.

## Nonpromotion Decision

E_pf(t,x)>0 is exactly the desired pointwise noncontact statement, so it cannot be assumed.

A positive or monotone spatial integral of E_pf does not rule out E_pf=0 at an isolated point.

The flux contains H_xx and the local balance is not a closed scalar parabolic inequality for E_pf.

No kappa is defined by dividing through E_pf or ||(H,s_pf H_x)||.

The identity is still useful:

The exact contact-observability floor converts the multiplier into a pointwise criterion once Xi-specific endpoint margins and a strict upper bound for the integrated bulk or equivalent energy-current numerator are supplied.

Route decision:

Retain the balance as an exact multiplier component, reject energy positivity alone as the descendant theorem, and return the pointwise step to the Xi contact-normal arithmetic hierarchy unless an independent trace/observability theorem is proved.

Open handoff:

On C_j^out=[t_(j+1),t_j]x[38,R_j], derive source-level endpoint lower margins and prove the strict bulk-budget inequality D_j(t)<2[sqrt(E_pf(t,38))+sqrt(E_pf(t,R_j))]^2/(R_j-38), where D_j=d/dt integral_38^R_j E_pf dx+2[J_pf(t,R_j)-J_pf(t,38)]. Do not define D_j by a relative quotient or assume the desired contact exclusion.

## Boundary

This artifact proves the local and integrated multiplier identities, exact parabolic-frequency derivatives, and positive outer-domain bulk coefficient. It also proves that bulk positivity or integrated energy control alone does not exclude a first-jet contact. It does not prove Xi endpoint lower margins or the strict bulk-budget inequality, certify Q209, prove the ray-aligned descendant theorem, prove Lambda<=0, prove PF-infinity, prove RH, or establish a Clay-prize conclusion.
