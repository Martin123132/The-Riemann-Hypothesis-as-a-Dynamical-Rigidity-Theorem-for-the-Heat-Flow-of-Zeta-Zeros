# Jensen-Window PF Newman Polymath-15 q=1 Finite-Height Real-Edge Remainder Gate

Date: 2026-07-31

Status: effective finite-height sign for the retained first-order q=1 real-edge model. This is not a proof of an Xi-level edge theorem, Lambda <= 0, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate.py
```

## Effective Domain

On q=1 and L>=50, t<=1/5000 and the saddle identity gives

```text
h=1/a<exp(-25)<1/72000000000.
```

All estimates below are uniform for -1<=p<=1. The proof enlarges the physical curve to the rectangular parameter box, so no monotonicity in p or L is assumed.

## Stable Terminal Logarithm

Write Z=z_N/S_a=-Q exp(W), with Q=exp[-pi i(p^2/2-p+3/8)]. The continuous branch W->0 is evaluated only after the h=0 and y=0 quotients have been removed. The exact stable identities are saved in the JSON artifact. Their analytic majorants give

```text
|W|<6h,                 |W_x|<2h^2,
|Z+Q|<12h,              |Z_x-i h theta Q/2|<10h^2.
```

The constant t*pi^2/64 is not bounded in isolation. It cancels algebraically against the -pi^2/16 part of Re(w^2). This is the essential finite-height cancellation.

The endpoint frame has the parallel cancellation

```text
|mu|<h^2,               |mu_x|<h^4.
```

Here the -pi/8 in K_rate cancels the +pi/8 from i(chi)/2 before absolute values are taken.

## Endpoint Derivatives

At p=+-1/2 the common numerator/denominator zero of C_0 is factored on complex disks of radius 1/4. The saved disk bound is

```text
[1.06656932105132389235282563401131576761036690288863598515 +/- 6.23e-57]
<3.
```

Cauchy on radius 1/8 and a 4096-box Arb cover of the complementary real intervals prove

```text
|F^(k)| < 3 k! 8^k, k=0,...,5,
(|F|,...,|F^(5)|) < (3,24,384,9216,294912,11796480).
```

Consequently, for C_1=F'''/(12pi^2),

```text
(|C_1|,|C_1'|,|C_1''|)<(86,2731,109227),
|H-F|<86h,
|H_x+hF'/(4pi)|<232h^2,
|H_xx-h^2F''/(16pi^2)|<800h^3.
```

## Four-Jet Budget

With A, B, M as in the cofinal gate, the complete normalized edge satisfies

```text
c=A+h e_0,                       |e_0|<100,
d_edge=hB+h^2 e_1,               |e_1|<250,
c_x=hB+h^2 e_2,                  |e_2|<250,
d_(edge,x)=h^2M+h^3 e_3,         |e_3|<900.
```

No quotient by A, c, or H is introduced. Direct expansion gives

```text
h^-2(c*d_x-d*c_x)-(A*M-B^2)=h[e_0M+A e_3-B(e_1+e_2)]+h^2[e_0e_3-e_1e_2]
```

Using |A|<4, |B|<3, and |M|<3, the linear constant is 5400 and the quadratic constant is 152500. Hence

```text
|a^2 J_edge/S_a^2-K_edge(p)|<5500h<1/10000000
```

## Finite-Height Sign

The preceding cofinal gate proved K_edge<-3/8000 on the whole cell. Combining the two strict inequalities yields

```text
a^2 J_edge/S_a^2<-3749/10000000<0
```

for q=1, L>=50, and -1<=p<=1 in the retained first-order model, including both removable C_0 points and every zero real projection.

## Boundary and Handoff

This closes the finite-height remainder requested by Formal Core 11.148 inside the first-order endpoint/Dirichlet model. It does not bound the omitted higher-order Xi approximation error, does not extend the sign to q>1, and does not splice the signed edge through adjacent cutoffs. Those are the next three gates before the edge block can be inserted into the cumulative Abel/contact scalar. No winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
