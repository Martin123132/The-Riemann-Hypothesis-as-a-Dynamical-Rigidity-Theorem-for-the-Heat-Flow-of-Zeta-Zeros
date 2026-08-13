# Finite Airy-Fresnel detuning ODE for the lower fold

Date: 2026-08-11

Status: exact finite-core detuning ODE and representative rigorous residuals
validated; not a proof of the fold-uniform source bound

Define the finite lower-fold transform

```text
G_Y(d)=integral_0^Y Ai(-lambda-y)
       exp(i*y^2/(4beta)-i*d*y)dy.                     (DO1)
```

The Airy equation and two integrations by parts give the exact inhomogeneous
ODE

```text
G_Y''/(4beta^2)+i(1+d/beta)G_Y'
 +(lambda-d^2+i/(2beta))G_Y=R_0(d)+R_Y(d),             (DO2)

R_0(d)=-Ai'(-lambda)+i*d*Ai(-lambda),                 (DO3)

R_Y(d)=exp(i[Y^2/(4beta)-dY])
 {Ai'(-lambda-Y)+i[Y/(2beta)-d]Ai(-lambda-Y)}.       (DO4)
```

Both endpoint currents remain explicit in (DO3)--(DO4); no endpoint or
Fresnel component is frozen or discarded.  For `Y=64` and all 84 detunings,

```text
1/(4beta^2)=[5.386082772039459846225327316998964406087808234379546019543736335661519e-8 +/- 3.63e-77]<5.4e-8, (DO5)

[0.9989660164058730268146412076928379403046804990694147652857241331770869 +/- 3.16e-71]
 <1+d/beta<
[1.001046516728601239526999504941188266479505191850955964832024665208645 +/- 3.69e-70],                 (DO6)

lambda-d^2>=[0.02648392089227129526190669315193748699964603836822134305616786739158831 +/- 1.84e-64]>0.026.       (DO7)
```

Thus the new equation is a nondegenerate first-order detuning transport with
a small but retained second-derivative term; it does not inherit the
`1/Delta` singularity of the frozen one-dimensional saddle expansion.

The identities

```text
G_Y'=-i integral_0^Y y f_d(y)dy,
G_Y''=-integral_0^Y y^2 f_d(y)dy                       (DO8)
```

were evaluated by rigorous complex ball quadrature at modes `39853`,
`39894`, and `39936`.  At all three representatives, the real and imaginary
balls of the independently assembled left side minus right side contain
zero.

The exact grouped compact-core value is `2pi*sum_m G_64(d_m)`, so (DO2)
provides a detuning-space route to all 84 modes while preserving their finite
sum.  The next obligation is an interval ODE/energy estimate or a sourced
variation-of-constants formula that bounds this complete discrete sum and
matches its outer-Morse limit.

Pi provenance: the chirp in (DO1) comes from the Kummer quadratic and the
detuning Fourier factor comes from finite Poisson summation.  No fitted
geometric normalization is introduced.

Proof boundary: exact canonical finite-core ODE, coefficient margins, and
three rigorous saved-height residual checks only.  No all-detuning ODE
enclosure, outer-chart join, complete `T_upper`, height-uniform source
theorem, `Lambda<=0`, RH, or prize-level conclusion is proved.
