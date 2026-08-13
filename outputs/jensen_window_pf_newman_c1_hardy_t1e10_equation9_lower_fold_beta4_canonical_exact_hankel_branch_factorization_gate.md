# Beta^-4 canonical fold: exact Hankel branch factorization

Date: 2026-08-13
Status: exact branch identity; not a proof of the ordinary-carrier splice

For `X>0`, put

```text
xi=2X^(3/2)/3,       A(X)=Ai(-X),       A_X=dA/dX.
```

The exact Bessel connection formulae give

```text
A(X)=sqrt(X)/(2sqrt(3))
 {e^(i*pi/6) H^(1)_(1/3)(xi)
  +e^(-i*pi/6)H^(2)_(1/3)(xi)},                       (HB1)

A_X(X)=-X/(2sqrt(3))
 {e^(-i*pi/6)H^(1)_(2/3)(xi)
  +e^(i*pi/6) H^(2)_(2/3)(xi)}.                       (HB2)
```

The `pi/6` phases are forced by the Bessel reflection identity at orders
`1/3` and `2/3`; they are not fitted phases.  If `U,V` are the real rational
weights from the beta-minus-four Airy reduction, define

```text
C_1(X)=1/(2sqrt(3)){
 U sqrt(X)e^(i*pi/6)H^(1)_(1/3)(xi)
 -V X e^(-i*pi/6)H^(1)_(2/3)(xi)},

C_2(X)=1/(2sqrt(3)){
 U sqrt(X)e^(-i*pi/6)H^(2)_(1/3)(xi)
 -V X e^(i*pi/6)H^(2)_(2/3)(xi)}.                    (HB3)
```

Then, exactly,

```text
U A+V A_X=C_1+C_2,             C_2=conj(C_1)          (HB4)
```

for real `X,U,V`.  The two Hankel functions are independent because

```text
W(H^(1)_nu,H^(2)_nu)=-4i/(pi*xi).                     (HB5)
```

Their standard leading forms identify the phase contract, without yet
asserting a quantitative remainder:

```text
C_1 ~[U X^(-1/4)+i V X^(1/4)]/(2sqrt(pi))
      *exp(i[xi-pi/4]),
C_2 ~[U X^(-1/4)-i V X^(1/4)]/(2sqrt(pi))
      *exp(-i[xi-pi/4]).                              (HB6)
```

After restoring the normal factor `exp(i[y^2/(4beta)-d*y])`, the two phase
derivatives are therefore

```text
y/(2beta)-d+sqrt(lambda+y),
y/(2beta)-d-sqrt(lambda+y).                           (HB7)
```

These are exactly the two ordinary-saddle equations already visible in the
cubic fold geometry.  Unlike a large-`X` WKB replacement, (HB1)--(HB4) remain
exact for every `X>0`.  The raw Hankel functions diverge at `xi=0`, but the
`sqrt(X)` and `X` prefactors in (HB3) give finite one-sided branch limits.  At
`X=0` the written products need either the regular Airy basis or separately
certified limiting values.  This distinction matters for a stable fold collar.

The remaining obligation is to choose and certify that collar, match the
appropriate exact Hankel branch to the endpoint-retaining logistic/Gamma
carrier on each ordinary corridor, and bound the branch remainder and
overlaps.  This gate proves no such quantitative splice, no `Q_K-T` bound,
no complete `T_upper`, no height-uniform theorem, no `Lambda<=0`, no
PF-infinity, no RH, and no prize-level conclusion.
