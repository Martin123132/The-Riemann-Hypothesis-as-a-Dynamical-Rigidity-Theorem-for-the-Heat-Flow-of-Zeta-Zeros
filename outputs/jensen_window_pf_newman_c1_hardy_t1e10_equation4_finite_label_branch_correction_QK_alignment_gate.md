# Finite-label branch correction and exact Q_K alignment

Date: 2026-08-27

Status: exact corrected finite-roster identity certified; not a proof of a
branch-sum or complete `Q_K` enclosure

The finite geometric contour prefix from the preceding gate cannot be
identified directly with the equation-(4) Kummer roster.  The obstruction is
an explicit branch term, not an unresolved sign convention.

For `s=1/2+it`, odd `alpha`, and

```text
r=exp(-i*pi/4),
z=1/2+q*r,
v=exp(i*pi/4)z=q+(1+i)/(2sqrt(2)),
c_alpha=pi*alpha*(1+i)/sqrt(2),
d=exp(-pi*t/4-i*pi/8),                              (BC1)
```

the exact label in the finite csc-prefix is

```text
T_alpha=-integral_C exp(-i*pi*z^2+i*pi*alpha*z)z^(-s)dz
       =-d integral_(R+i0) exp(-pi*v^2+c_alpha*v)v^(-s)dv. (BC2)
```

Put

```text
I_alpha(gamma)=integral_0^infinity x^(-s)exp(-pi*x^2+gamma*x)dx,
a=1/4-it/2,  b=3/4-it/2,  xi=i*pi*alpha^2/4,
A_alpha=Gamma(a)pi^(-a) 1F1(a;1/2;xi)/2,
B_alpha=c_alpha Gamma(b)pi^(-b) 1F1(b;3/2;xi)/2.    (BC3)
```

Then `I_alpha(c_alpha)=A_alpha+B_alpha` and
`I_alpha(-c_alpha)=A_alpha-B_alpha`.  The upper boundary value of `v^(-s)`
on the negative real axis gives `exp(-i*pi*s)=-i exp(pi*t)`, hence

```text
T_alpha=-d[(A_alpha+B_alpha)-i exp(pi*t)(A_alpha-B_alpha)]. (BC4)
```

After simplifying the odd-label phase and the gamma reflection in Appendix
A, the finite A33 source label is

```text
S_alpha=i exp(pi*t)d(A_alpha-B_alpha).               (BC5)
```

Therefore the exact relation is

```text
T_alpha=S_alpha-d I_alpha(c_alpha),
S_alpha-T_alpha=d I_alpha(c_alpha).                  (BC6)
```

The missing term in (BC6) is precisely what is lost if the integration
constant obtained on one side of the `Log(v^2)` cut is carried across the
whole horizontal endpoint without a jump.  This diagnosis does not require
assigning an unstated branch convention to the source: (BC2)--(BC6) compare
the explicit exact RSI label and explicit A33 Kummer label themselves.

Let `Hardy_t[X]=2 Re[exp(i*theta(t))X]`.  Odd-label Kummer reflection gives
both `exp(-i*pi/8)Phi1` and `exp(-i*pi/8)Phi2` real, and exact phase algebra
then yields

```text
Hardy_t[T_alpha]
 =-exp(-pi*t/4)|Gamma(1/4+it/2)|pi^(-1/4)
   Re[exp(-i*pi/8)Phi1],                              (BC7)

Hardy_t[S_alpha]
 =2pi^(5/4)alpha exp(-3pi*t/4)
   Re[exp(-i*pi/8)Phi2]
  /[(1+exp(-2pi*t))|Gamma(1/4+it/2)|]
 =H(t)K_t(alpha).                                    (BC8)
```

Thus the exact finite RSI prefix projects through `Phi1`, while the physical
Kummer label projects through the thermal `Phi2` term.  At the independent
rigorous witness `t=5, alpha=3`,

```text
Hardy_t[T_alpha]=[0.52744285057894241284578776810412557038591177996147660613597443224 +/- 7.61e-67],
Hardy_t[S_alpha]=[-0.52744224874646436532866272823041177777625373079527672835527672240 +/- 2.69e-66],
Hardy_t[S_alpha-T_alpha]=[-1.0548850993254067781744504963345373481621655107567533344912511546 +/- 4.20e-65]. (BC9)
```

The last ball excludes zero, so the uncorrected identification is false.

For the actual roster, define

```text
B_W=sum_(alpha=A,A+2,...,B) d I_alpha(c_alpha),
P_W=sum T_alpha,  S_W=sum S_alpha.                   (BC10)
```

There are `2481423` labels and, exactly,

```text
S_W=P_W+B_W,
H(t)Q_K=Hardy_t[S_W]
       =Hardy_t[-sum_(n=79789)^(2561211) n^(-s)
                 +C_79788-C_2561211+B_W]. (BC11)
```

Equation (BC11) is the corrected finite-contour route to the physical roster.
It uses no infinite A21 interchange, but the new finite branch sum `B_W` must
be retained and enclosed.  The next stage is to combine the released
Dirichlet block with the exact classical upper component before norms, then
derive a common-contour or summation representation for `B_W+C_n--C_m`.

Pi provenance: every `pi` comes from the original sine/Gaussian RSI kernel,
the rotation by `pi/4`, the Riemann-Siegel phase, or standard gamma/Kummer
normalization.  No fitted or unexplained circle constant is inserted.

Proof boundary: exact finite-label reduction, explicit branch correction,
corrected `H(t)K_t(alpha)` projection, and the augmented actual-window
identity only.  No enclosure of `B_W`, `Q_K`, `D_K`, `Delta_KU`, or `J_Z`, no
released-pole cancellation theorem, non-A bound, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
