# Exact signed-pair half-domain logistic-Morse transform

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of a quantitative amplitude enclosure

For `m>0`, introduce two genuinely different Fresnel coordinates

```text
q_D^(sigma)(m,x)=sqrt(x/2)(D-2 sigma m/x),
sigma in {+1,-1}, D in {A,B}.                       (PM1)
```

Completing the square for the `+m` and `-m` Fourier characters gives

```text
P_(sigma m)=[E_B^(sigma)-E_A^(sigma)]/(i*pi)
 +sigma m sqrt(2/x)[F(q_B^(sigma))-F(q_A^(sigma))],

J_(sigma m)=(-1)^m exp(-i*pi*m^2/x)P_(sigma m)/x.     (PM2)
```

Although the Fresnel coordinates in (PM1) are different, the parity and
outer phase in (PM2) are identical.  Hence the symmetric pair is exactly

```text
J_m+J_-m=(-1)^m exp(-i*pi*m^2/x)(P_m+P_-m)/x.         (PM3)
```

Use the global logistic coordinate

```text
r=t/(2*pi*m^2), v=(1-x)/(r*x),
s=sgn(v-1)sqrt(2[v-1-log v]).                          (PM4)
```

The phase remains globally exact,

```text
psi_m(x)-psi_m(x_m)=-t*s^2/4.                         (PM5)
```

At the half-boundary `x=1/2`,

```text
v_H=2*pi*m^2/t,
s_H=sgn(v_H-1)sqrt(2[v_H-1-log v_H]).                 (PM6)
```

Since `x:0->1/2` maps monotonically to `s:+infinity->s_H`, the pair in the
normal form of Section 11.354 has the exact one-phase representation

```text
K_m+K_-m=(-1)^m exp(i psi_m(x_m))
 integral_(s_H)^infinity exp(-i*t*s^2/4)A_m^pair(s)ds,

A_m^pair=r*x*(dv/ds)[x(1-x)]^(-1/4)(P_m+P_-m).        (PM7)
```

This is the phase-adapted cancellation object that was missing from the
tangent B calculation.  It keeps the negative completion in the same
Gaussian integral as its positive partner; the `1/m` endpoint current has
already cancelled before an absolute value is taken.

At `t=10^10`, `sqrt(t/(2*pi))` lies strictly between `39894` and `39895`.
The boundary coordinates are

```text
s_H(39894)=-0.00001143224830869833396525762486762794157524439977986905477397573085287735,
s_H(39895)=0.00003870020326729481772904816690773252795751152173839883173741857094783794.                       (PM8)
```

Thus modes `1..39894` integrate across the Gaussian saddle, while every
outer pair `m>=39895` starts on its nonstationary side.  The negative partner
does not create another saddle; it changes only the exact pair amplitude.

Equation (PM7) also explains why midpoint roster parity could not cancel the
Fourier pairs coefficientwise.  The cancellation is distributed through
the signed Fresnel currents along the complete `s` contour.

Pi provenance: every `pi` comes from completing the equation-(9)
Kummer/Fourier phase, the logistic saddle scale, or the exact reflection
phase.  No fitted constant is used.

Proof boundary: exact signed-current algebra, half-domain orientation, and
one-phase pair transform only.  No bound for `A_m^pair`, nonlinear B
crossing, outer-pair sum, or A-fold splice is proved.  No theorem establishing
`T_upper`, `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion follows here.
