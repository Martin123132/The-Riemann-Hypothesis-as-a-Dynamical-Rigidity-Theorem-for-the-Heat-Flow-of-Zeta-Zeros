# Jensen-Window PF Compound Order-Four First-Summand Curvature Bridge

Date: 2026-07-12

Status: exact reduction with a proved gap floor, a compact curvature
certificate, a closed full-kernel perturbation budget, and one open analytic
ray. This is
not a proof of order-four entry, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_first_summand_curvature_bridge`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_first_summand_curvature_bridge.py
```

## Stable Coordinate

For the continuous first-summand family put

```text
g(t)=d(t)^2-x(t)^2*d(t-1)*d(t+1), x(t)=exp(-B(t)), d(t)=1-x(t)
ell(t)=log(d(t)),
J(t)=2*B(t)-ell(t-1)+2*ell(t)-ell(t+1).
```

Then the cancellation in the order-three gap factors exactly:

```text
g(t)=d(t)^2*(1-exp(-J(t))).
```

With `phi(z)=1/(exp(z)-1)` and
`chi(z)=exp(z)/(exp(z)-1)^2`, differentiation gives

```text
ell''=phi(B)*B''-chi(B)*(B')^2
K(t)=(log g)''=2*ell''+phi(J)*J''-chi(J)*(J')^2
J^(r)=2*B^(r)-ell^(r)(t-1)+2*ell^(r)(t)-ell^(r)(t+1), r=0,1,2
```

This replaces a difference of two terms of size `t^-2` by positive
coordinates of sizes `d~t^-1` and `J~t^-1`.

## Proved Gap Floor

The existing first-summand cumulant wall, negative `H'''`, and continuous
scaled-curvature growth imply the following. Two extra interval-Simpson
blocks cover `0.925<=u<=0.9264`, with minimum shifted-buffer margin
`5.59711570013466157765892076575866435221556140951846016072971E-4`. Combined with the
existing compact and ray theorem, they give `B_1'<0` and increasing
`(2t+1)B_1(t)` early enough for the complete real tail.

```text
1/(2*m+1)<=B_1(m)<=3/(2*m-1), integer m>=319
1/(2*t+3)<=B_1(t)<=3/(2*t-3), real t>=319
J_1(t)>=2/(2*t+3)-2/(2*t-1)+(t-3)/(6*t^2)>=1/(7*t), t>=319
```

After `t=319+n`, the final numerator is

```text
4*n**3 + 3412*n**2 + 955637*n + 87486770.
```

Every coefficient is positive for `n>=0`. Independently, the repaired Arb
coefficients and certified moment perturbation give the cross-check

```text
J_1(319)-1/(7*319) > 2.69590680165062467909584378012170371110591865709957217006698E-3.
```

Thus `J_1(t)>=1/(7*t)` for every real `t>=319`. This is already a
theorem and is not part of the remaining curvature ceiling.

## Localized Curvature Form

For `r=0,1,2`, the tent representation and Taylor's theorem give

```text
j_0=2*B-ell''; j_1=2*B'-ell'''; j_2=2*B''-ell''''
E_r=(1/12)*sup_|s|<=1 |ell^(r+4)(t+s)|; |J^(r)-j_r|<=E_r, r=0,1,2
```

Because both `phi` and `chi` decrease on the positive half-line,

```text
U=2*ell''+phi(j_0-E_0)*max(j_2+E_2,0)-chi(j_0+E_0)*max(abs(j_1)-E_1,0)^2
j_0>E_0 and U<=7/(2*t^2) imply J>0 and K<=7/(2*t^2)
```

This is the intervalization coordinate: all central terms share the same
real parameter, while shifted dependence is confined to three explicit
derivative envelopes. It avoids independent enclosure of nearly equal
adjacent moments on the large-parameter tail.

## Compact Curvature Theorem

A deterministic cache of outward-rounded Arb quadrature tiles now proves

```text
K_1(t)<=7/(2*t^2), 319<=t<=V'(2).
```

The cache contains `107452` derivative tiles and
assembles `1073` positive central blocks.
Its weakest margin is `1.14426650120038260530285332880335545546623246507456181669743E-10`, and the
largest certified `t^2 U(t)` upper bound is `3.33723355756790554314260494045514876254074160915177414272085E+0`.
This is a real-parameter interval theorem, not finite sampling.

## Remaining Analytic Ray

For `k=n+3`, the centered second difference has the exact tent form

```text
P_n=Delta^2 log(g)(k)=integral_[-1,1](1-|s|)*K(k+s)ds, k=n+3
```

so the complete first-summand theorem would follow after proving only

```text
K_1(t)<=7/(2*t^2) on the mode ray u>=2.
```

The checked Gaussian-cumulant target expands the standardized tilted
partition exactly through `epsilon^6`, records the alternating factorial
cumulant signature through order eight, and gives explicit candidate
corridors that clear seven conditional `t+-2` ray collars. A 1.8-million-block
finite certificate and coefficient-positive asymptotic theorem now prove the
formal corridors for every `u>=2`. Only the exact-minus-formal central
remainder and two-tail estimates remain open.

It would give

```text
K(t)<=7/(2*t^2) => P_n^(1)<=7/(2*(k^2-1))<=18/(5*k^2), k>=320
```

using `-log(1-z)<=z/(1-z)` and the exact margin

```text
18/(5*k^2)-7/(2*(k^2-1))=(k^2-36)/(10*k^2*(k^2-1)).
```

## Full-Kernel Budget

Write `M_j=M_j^(1)*(1+delta_j)` with `0<=delta_j<=2/j^6`.
The exact Lipschitz transfer through `B`, `ell`, and `J` gives

```text
a_j=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)); |B_j-B_j^(1)|<=a_j
a_j<=8/(j-1)^6; |ell_j-ell_j^(1)|<=32*j/(j-1)^6; |J_j-J_j^(1)|<=176*j/(j-2)^6<=1/(14*j)
|log(g_j)-log(g_j^(1))|<=2528*j^2/(j-2)^6, j>=319
|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2), k>=320
```

After `k=320+n`, the last comparison is certified by the
coefficient-positive polynomial

```text
n**6 + 1902*n**5 + 1482055*n**4 + 604691300*n**3 + 135889991935*n**2 + 15877422036942*n + 748002501678169.
```

Consequently the sole continuous target would imply

```text
P_n<=18/(5*k^2)+2/(5*k^2)=4/k^2, k=n+3>=320.
```

## High-Precision Scout

These saddle-centered rows use finite-upper mpmath quadrature. They are
not interval enclosures and are not promoted to the continuous theorem.

| t | saddle | t J(t) | t^2 (log g)'' | margin below 7/2 |
|---:|---:|---:|---:|---:|
| `319` | `0.926843126422534373746647515760644309` | `1.00285142495154167220802555038171167` | `1.79490588460442263240447097960303823` | `1.70509411539557736759552902039696177` |
| `400` | `0.97888094701689181843228254728997258` | `1.05780877651205911009808887022116525` | `1.94559368840472305388749756244874879` | `1.55440631159527694611250243755125121` |
| `700` | `1.10694529940786730130387507834158606` | `1.17469456249194627832611471213263743` | `2.25147275293958900562968385660131873` | `1.24852724706041099437031614339868127` |
| `1000` | `1.18793754114521619465795554468871379` | `1.23618419484572746294853754785757159` | `2.39813564902961774575711336224201266` | `1.10186435097038225424288663775798734` |
| `2000` | `1.34400782600033672514953074675098446` | `1.33216398202528785440896450054418576` | `2.59807241547915554978310950735453838` | `0.901927584520844450216890492645461619` |
| `5000` | `1.5485738239457640734903380663464682` | `1.42442320146611216441780421169711298` | `2.74905394348829534408302000954709661` | `0.750946056511704655916979990452903385` |
| `20000` | `1.85736880805336508625944338043977827` | `1.51923688461353714866705160022343311` | `2.85739663453953075903928820213676513` | `0.642603365460469240960711797863234874` |
| `100000` | `2.21836163013573509263544082755573846` | `1.5932949879963637618540811102735149` | `2.91156137414444764804708033458057033` | `0.588438625855552351952919665419429673` |

The sampled scaled curvature rises from about `1.795` toward the
formal limit `3`, leaving more than `0.58` at the largest sample.
The compact range is now rigorously interval-certified. The formal cumulant
model and its first omitted coefficient layer are closed globally. The next
proof task is the cancellation-preserving exact-density central and two-tail
theorem at the explicit finite and asymptotic remainder budgets.

```text
outputs/jensen_window_pf_compound_order4_condensation_gate.md
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md
outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/signed_hankel_jensen_dependency_graph.md
```
