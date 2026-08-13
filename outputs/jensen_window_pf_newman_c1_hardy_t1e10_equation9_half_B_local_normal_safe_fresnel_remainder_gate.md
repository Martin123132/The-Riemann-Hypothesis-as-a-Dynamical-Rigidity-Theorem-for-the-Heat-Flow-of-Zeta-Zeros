# Normal-safe Fresnel remainder for local B modes

Date: 2026-08-13

Status: saved-height interval certificate for two specified arcs; not a
proof of the complete local B current

The exact local-domain gate locates the balance points between the outer
Gaussian and normal Fresnel derivatives.  Choose the exact rational split
points

```text
x_621=24238523659/100000000000000 inside ['0.00024238523658', '0.00024238523660'],
x_622=48583287419/200000000000000 inside ['0.00024291643709', '0.00024291643710']. (NF1)
```

On `0<x<=x_621`, `q_B(621,x)` is negative and increasing, while on
`x_622<=x<=1/2`, `q_B(622,x)` is positive and increasing.  The certified
split values are

```text
q_621=[-18.202517652441687900608545544598399434940741531089673408624756274480881877836554322927284666333097463684370558 +/- 1.25e-106],
q_622=[14.529724498597999445466432218645430050656146139278559294075788340148207219989786471784717994115228039171758335 +/- 1.25e-106].                       (NF2)
```

Thus the third-order Fresnel expansion is uniformly valid on both arcs.
Its completed-current remainder from Section 11.364 is

```text
|R_ext,m(x)|<=15*m*x^2/[4*pi^4*|Bx/2-m|^7].          (NF3)
```

Restoring the exact Kummer weight and the real half-domain projection gives

```text
E_m^N=2*(pi/(32t))^(1/4)
 integral_arc [x(1-x)]^(-1/4)|R_ext,m(x)|dx.          (NF4)
```

This is the normalization used by equation (9), not the raw completed
current.  Put `d=|Bx/2-m|` and
`h(x)=x^(7/4)(1-x)^(-1/4)`.  Since

```text
h'(x)/h(x)=(7-6x)/[4x(1-x)]>0 on 0<x<=1/2,           (NF5)
```

each dyadic `d` panel is bounded at the correct monotone endpoint, while
`integral d^(-7)dd` is evaluated exactly.  Arb outward rounding gives

```text
E_621^N <= [4.0054166841610652718047561300773683707259239441567150642283569872789896222883816218382739353452218055482083445e-11 +/- 3.22e-120],
E_622^N <= [1.5474209308082053751687683331843082485197181315184665661496285054284174667275841323992414476090060749317139280e-10 +/- 1.96e-119],

E_621^N+E_622^N
 <= [1.9479625992243119023492439461920450855923105259341380725724642041563164289564222945830688411435282554865347624e-10 +/- 1.92e-119] < 1.97e-10. (NF6)
```

This uses only 12 and 23
analytic panels respectively; there is no oscillatory numerical quadrature.
The combined bound is below `2.30e-5` of the reference physical target
`8.6e-6`.

The result closes only the positive local B-tail remainders on the two
normal-safe arcs.  On the complementary outer-safe arcs, the exact 621/622
currents must remain combined with the pole-subtracted cotangent background;
its raw absolute value is known to be too large to triangle away.  Negative,
zero, outer-positive, and A-fold sectors are not reassigned here.

Pi provenance: every `pi` in (NF3)--(NF6) comes from the equation-(9)
Fresnel phase or the paper normalization `(pi/(32t))^(1/4`; Arb evaluates
that same constant directly.  No geometric or fitted value is inserted.

Proof boundary: the two stated saved-height positive-B normal-safe remainder
arcs only.  No complementary outer-safe local-current estimate, complete B
face estimate, A-fold splice, complete paired residual, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
