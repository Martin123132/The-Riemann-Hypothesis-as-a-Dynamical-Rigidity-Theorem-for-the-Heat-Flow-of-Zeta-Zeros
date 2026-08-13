# Universal logistic Morse phase and characteristic fold reduction

Date: 2026-08-10

Status: exact phase and fold-obstruction reduction validated; not a proof of the fold-uniform remainder

For a positive Poisson mode define

```text
w=(1-x)/x,  r=t/(2*pi*m^2),  v=w/r,
x=1/(1+r*v).
```

The reduced Kummer phase becomes exactly

```text
psi_m(x)-psi_m(x_m)=-(t/2)(v-1-log v).                  (LM1)
```

It is independent of `m`.  The signed global Morse coordinate

```text
s=sgn(v-1)sqrt(2(v-1-log v))                            (LM2)
```

maps `(0,infinity)` monotonically onto the real line and gives

```text
psi_m(x)-psi_m(x_m)=-t*s^2/4,
dv/ds=v*s/(v-1),   (dv/ds)|_(s=0)=1.                    (LM3)
```

The removable local chart is

```text
v(s)=1+s+s^2/3+s^3/36-s^4/270+s^5/4320+O(s^6).         (LM4)
```

The alpha current can be kept grouped without separating its boundary and
Fresnel series.  With `q=sqrt(x/2)(alpha-2m/x)`, it is exactly

```text
J_m=(-1)^m e^(-i*pi*m^2/x)
 int_(q_A)^(q_B)[m*sqrt(2)x^(-3/2)+q*x^(-1)]
 e^(i*pi*q^2/2)dq.                                     (LM5)
```

For a full interior saddle, the complete-line Fresnel integral and the
universal Gaussian in `s` give the raw leading main

```text
2^(5/4)(t/pi)^(1/4)/sqrt(m).                            (LM6)
```

The paper normalization `(pi/(32t))^(1/4)` cancels the prefactor in (LM6)
exactly, leaving `1/sqrt(m)`.  Integer parity also cancels
`exp(-i*pi*m^2)`, leaving the classical phase carrier

```text
exp(i{t/2[log(t/(2*pi))-1]-t*log(m)}).                 (LM7)
```

After the paper factor `exp(-i*pi/8)`, this is the leading Riemann--Siegel
phase `theta_0(t)-t*log(m)`.  This identifies the correct classical carrier,
but not yet its real/conjugate assembly or remainder.

The endpoint transition has a second scale.  Let `kappa_m(C)` be the change
of the inner endpoint coordinate across one outer Gaussian width:

```text
kappa_m(C)=(2/sqrt(t)) partial_s q_C|_(s=0),
Delta_m(C)=1-pi*kappa_m(C)^2/2.                         (LM8)
```

When the joint saddle lies on `alpha=C`,

```text
Delta_m(C)=(2*pi*m^2-t)/(2*pi*m^2+t).                   (LM9)
```

At the exact coalescence
`m=sqrt(t/(2*pi))`, `C=sqrt(8t/pi)`, one has

```text
kappa=-sqrt(2/pi),   Delta=0.                           (LM10)
```

Thus the full two-dimensional Hessian and reduced x saddle remain
nondegenerate, but the boundary tangent is characteristic.  A plain
one-dimensional Gaussian remainder that freezes the Fresnel endpoint is not
uniform through this fold.

At `t=10^10`, 130-digit reconstruction gives

```text
39853..39894: Delta_m(A)<0,
39895..39936: Delta_m(A)>0,

Delta_39894(A)=[-6.2665739085848621125804364673098811428099067711555368331915967423319074248146462706834440160035120757084743292e-6 +/- 2.62e-110],
1/|Delta_39894(A)|=[159576.83011287155125040879588326167809116390545664550595069633726549018036314721481056089437426917728684875842 +/- 6.66e-100].
```

The sign change is exactly the certified `42+42` split.  Any proposed bound
that divides by `Delta` loses more than `159000` already at the central
height and cannot be the desired uniform theorem.

The next route is a fold-uniform estimate of (LM5) after (LM1)--(LM4), or an
equivalent two-variable endpoint normal form, retaining both members of the
grouped current.  Away from the characteristic chart, ordinary stationary
phase and integration by parts remain available.

Proof boundary: exact universal Morse algebra, leading classical carrier,
and a certified characteristic-fold obstruction at `t=10^10`.  No
fold-uniform remainder, complete `T_upper` identification, source-aligned
height-uniform error, `Lambda<=0`, RH, or prize-level conclusion is proved.
