# Newman Positive-Boundary Diagonal Exhaustion Gate

Date: 2026-07-24

Status: exact independent compact-exhaustion theorem and route guard.
The growing Xi shell family remains open; this is not a proof of
`Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.py
```

## Diagonal Exhaustion

Positive-boundary attainment gives

```text
If Lambda>0, then 0<Lambda<=1/5 and H_Lambda has a finite real multiple zero c
```

The fixed finite witness changes the quantifier geometry:

```text
For any sequences 0<delta_j<=1/5 with delta_j->0 and R_j>0 with R_j->infinity, Lambda<=0 iff (H_t,H_t')!=(0,0) on [delta_j,1/5]x[-R_j,R_j] for every j
If Lambda>0 with boundary collision c, then every sufficiently large j satisfies delta_j<=Lambda and R_j>=|c|, so the same fixed point (Lambda,c) lies in every later rectangle
```

A particularly simple choice, composed with the certified
`|x|<=38` core, is

```text
With delta_j=1/(5j) and R_j=38+j, Lambda<=0 iff (H_t,H_t')!=(0,0) on [1/(5j),1/5]x[-(38+j),38+j] for every j>=1
Under Lambda>0, J=ceil(max(1,1/(5*Lambda),|c|-38)) places (Lambda,c) in every linear-exhaustion rectangle j>=J
Because (H_t,H_t')!=(0,0) is already certified for 0<=t<=1/5 and |x|<=38, the linear exhaustion only leaves 1/(5j)<=t<=1/5 and 38<|x|<=38+j
```

| j | delta_j | radius | remaining shell |
|---:|---:|---:|:---|
| 1 | 1/5 | 39 | 38<|x|<=39 |
| 2 | 1/10 | 40 | 38<|x|<=40 |
| 5 | 1/25 | 43 | 38<|x|<=43 |
| 10 | 1/50 | 48 | 38<|x|<=48 |
| 100 | 1/500 | 138 | 38<|x|<=138 |

For each fixed j, no common zero on its compact rectangle is equivalent to min[H_t(x)^2+H_t'(x)^2]>0 there

## Quantifier Correction

```text
No coupling between delta_j and R_j is logically required: R_j may diverge arbitrarily slowly, for example R_j=38+log(1+j), independently of delta_j->0
The radius 4*pi*exp(25/delta) remains a valid stronger per-delta localization that contains every possible contact above that time floor, but it is not required by the cofinal positive-boundary contradiction
```

The earlier `4*pi*exp(25/delta)` radius remains useful when one
wants a single fixed-delta rectangle containing every possible
contact. It is not the minimal cofinal contradiction. One fixed
finite boundary contact is eventually captured by any spatial
radius tending to infinity.

## Arbitrary-Height Field Guard

Fix `c>0` and a target field `b<0`. Put

```text
a_c^2=c^2*(3-b*c)/(1-b*c)
P_c(z)=(z^2-c^2)^2*(z^2-a_c^2)
B_c=1/c+2*c/(c^2-a_c^2)=b at the positive double zero c
K_c=1/(2*c^2)+1/(c-a_c)^2+1/(c+a_c)^2=b^2-3*b/c+5/(2*c^2)
```

The shifted backward-heat flow satisfies

```text
F_(lambda+tau)=exp(-tau*d_z^2)P_c solves partial_tau F=-partial_z^2 F
z_+/-=c+/-sqrt(2*tau)+2*b*tau+O(tau^(3/2))
For tau>=0, exp(-tau*D^2) preserves real-rootedness because it is the coefficientwise limit of [(1-sqrt(tau/n)D)(1+sqrt(tau/n)D)]^n, and each first-order factor preserves real-rootedness by derivative interlacing
For tau<0 sufficiently close to zero the two roots near c are nonreal, while every root is real for tau>=0; after the time shift, this polynomial flow has Newman-style boundary lambda
```

At the classical field value:

```text
For b=-pi/8, a_c^2=c^2*(24+pi*c)/(8+pi*c), B_c=-pi/8, and K_c=pi^2/64+3*pi/(8*c)+5/(2*c^2)->pi^2/64
```

| c | a_c | a_c-c | B_c | K_c |
|---:|---:|---:|---:|---:|
| 1 | 1.56078839045 | 0.560788390447 | -0.392699081699 | 3.83230981386 |
| 10 | 11.8571804569 | 1.85718045689 | -0.392699081699 | 0.297022293277 |
| 100 | 102.453153983 | 2.45315398291 | -0.392699081699 | 0.166243541218 |
| 1000 | 1002.53679334 | 2.53679334432 | -0.392699081699 | 0.155393166012 |

Thus local field, center drift, positive stiffness, evenness,
and backward-heat evolution do not bind collision height to
boundary time. Any such bound must use global Xi structure.

## Live Handoff

```text
Use the independent diagonal exhaustion as the minimal compact endgame: prove delta-dependent C1 separation on its growing finite shells, or derive an Xi-specific nodal-flux/resolvent continuation across them. Do not impose an exponential radius/time-floor coupling, and do not infer any radius bound from local field, drift, or stiffness data alone
```
