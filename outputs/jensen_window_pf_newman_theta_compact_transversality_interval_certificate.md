# Newman Theta Compact-Transversality Interval Certificate

Date: 2026-07-24

Status: rigorous compact contact exclusion. This is not a proof
of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_compact_transversality_interval_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_compact_transversality_interval_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_compact_transversality_interval_certificate.py
```

Current result:

```text
validated Newman theta compact-transversality interval certificate: 8 rows, 0 issues, 4 exact identities/inequalities, 1900 certified boxes, 10 adaptive subdivisions, 0 unresolved boxes, 1 compact no-contact theorem, 1 origin composition, 1 open high-frequency handoff
```

## Certified Theorem

For every

```text
0<=t<=1/5 and |x|<=38
```

the exact Newman heat-flow transform satisfies

```text
(H_t(x),H_t'(x))!=(0,0).
```

For `|x|<=1/4` this follows from the independent rational
moment margin `248/371925`. Evenness reduces the remaining
region to the certified rectangle `1/4<=x<=38`.

## Interval Method

H_(1,t)(x)=integral_0^infinity exp(tu^2)*phi_1(u)*cos(xu)du

J_1=16x^4H_(1,t), J_1'=64x^3H_(1,t)+16x^4H_(1,t)'

The retained integral is evaluated by Arb on `0<=u<=2`.
For every transform moment used in the Taylor model,

```text
t*u^2+9u-pi*exp(4u)<-2979-u
integral_2^infinity u^m*exp(tu^2)*phi_1(u)du<20*9!*exp(-2979)<10^-800
```

The omitted tail is therefore included with radius `1e-800`.
A second-order time and third-order frequency Taylor model
uses exact mixed-jet integrals and positive moment remainders.

## Partition

```text
precision bits=160
initial boxes=1890
evaluated boxes=1910
certified leaf boxes=1900
adaptive subdivisions=10
maximum depth=1
unresolved boxes=0
value/derivative branches=1115/785
minimum certified ratio lower=[1.7521852560914668686294040482508965391487117311 +/- 2.31e-47]
minimum box=t:[9/50,1/5], x:[757/20,38], branch=derivative
```

Every leaf proves `|J_1|>B_0` or `|J_1'|>B_1`; the exact
componentwise tail theorem then excludes full contact.

## Proof Boundary

No value beyond `|x|=38` is certified here. The next theorem
must use the existing corrected Riemann-Siegel partition:
the dominant and oscillatory-zeta regions are already closed,
while the scaled critical phase-transversality layer remains
open. The compact theorem does not imply `Lambda<=0` or RH.
