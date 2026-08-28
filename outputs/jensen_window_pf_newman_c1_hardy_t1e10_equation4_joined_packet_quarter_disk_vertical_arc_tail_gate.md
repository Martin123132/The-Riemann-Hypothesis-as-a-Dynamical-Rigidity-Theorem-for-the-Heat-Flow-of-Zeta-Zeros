# Quarter-disk vertical/arc/tail representation for the joined packet

Date: 2026-08-27

Status: exact contour identity and rigorous production upper-arc enclosure certified; lower-cell and complete `D_K` assembly open.

Let

```text
F_s(u)=u^(-s)exp(i*u^2/(4*pi))Delta_q(u),
M_s=integral_0^infinity F_s(u)du,
E_s=exp(pi*t/2+i*pi/4)(2*pi)^(s-1),
R=2*pi*Y.
```

Close the fourth-quadrant quarter disk with the principal branch.  The small
circle at zero vanishes because `Re(1-s)=1/2`.  If `A_Y^cw` is the clockwise
arc from `R` to `-iR` and `T_Y` is the positive-real tail from `R`, contour
orientation gives

```text
M_s=integral_[0 to -iR]F_s(u)du-A_Y^cw+T_Y.         (QD1)
```

On `u=-2*pi*i*y`, with `0<y<Y`, the principal logarithm gives exactly

```text
E_s u^(-s)(-2*pi*i*dy)=y^(-s)dy,
exp(i*u^2/(4*pi))=exp(-i*pi*y^2),
Delta_q(-2*pi*i*y)=D_W(y).                          (QD2)
```

Thus, writing `C_Y=-E_s A_Y^cw`,

```text
S_W=V_Y+C_Y+E_s T_Y,
V_Y=integral_0^Y y^(-s)exp(-i*pi*y^2)D_W(y)dy.      (QD3)
```

Combining (QD3) with the already-certified two-boundary join at
`L=621.5`, `U=39936.5` cancels the whole interior current before any norm:

```text
U_unowned=integral_0^L g_W(y)dy+C_U+E_s T_U.        (QD4)
```

For the production upper boundary put `u=-2*pi*i*U*exp(i*delta)`.  Since `U`
and `q_-` are half integers and the label count `N=2481423` is odd, the
endpoint carrier separates exactly:

```text
C_U=exp(i*[pi-t*log(U)-pi*U^2]) C_U^rot,

C_U^rot=sqrt(U) integral_0^(pi/2)
  exp(X(delta)) (1+w(delta)^N)/(1+w(delta)) ddelta. (QD5)
```

The stable coordinates used by Arb are

```text
X(delta)=(t+2*pi*U^2-q_-R+i/2)delta
 -i*pi*U^2[expm1(2i*delta)-2i*delta]
 +i*q_-R[expm1(i*delta)-i*delta],

w(delta)=exp(i*R*expm1(i*delta)).                   (QD6)
```

They expose the small real endpoint slope
`206.8358840647606...` instead of subtracting three unresolved
`10^10`-scale interval terms.

Arb proves that the real action has its unique maximum at

```text
delta_* in [0.000143589738918530181449813713082292598821747619539069167538519417662969409176984 +/- 6.24e-82],
A(delta_*) in [0.01979967368031121930651718904437006522489619007936360954454363153257538 +/- 9.18e-72].
```

The complete signed arc, including the tiny initial geometric-numerator
layer, the finite Taylor-coordinate error, and the omitted compact tail,
satisfies

```text
C_U^rot in
  [0.0413826555194463161241427 +/- 2.97e-26]
 +i*[0.0401581142560263590705370 +/- 5.01e-26],

|C_U|=|C_U^rot| in [0.05766453258668822459759736 +/- 8.77e-27] < 0.057665.    (QD7)
```

The separate positive-real tail has

```text
log10 |E_s T_U| <= [-1873216239.548398360977478922000587809792808100805199583484181462058805664365019375755 +/- 2.15e-76],              (QD8)
```

so it is negligible on every displayed corridor scale.  An altered
`t=3.25`, five-label contour replay checks (QD1)--(QD3) directly.  The
independent checker changes the altered height, roster, endpoint, precision,
Taylor degree, initial split, and Arb paneling, and requires rigorous overlap
with (QD7).

The value `0.057665` is usefully comparable with the width of the stored
`D_K` corridor, but `C_U` is not `D_K`.  The lower finite cell, ordinary
Gamma-subtracted packet, grouped natural A lift, and all remaining equation-
(4) ownership terms must still be assembled with their signs intact.

Pi provenance: every `pi` in (QD1)--(QD8) comes from the inherited
Riemann-Siegel Gaussian, the fixed Mellin scaling `R=2*pi*Y`, the principal
quarter turn, or the odd Fourier roster.  No fitted geometric constant is
inserted.

## Proof boundary

Exact finite quarter-disk identity, exact vertical-current normalization,
cancellation-preserving unowned-cell reduction, and rigorous production
upper-arc/tail enclosure only.  No lower-cell enclosure, complete joined
packet, `J_Z`, `D_K`, non-A, all-height, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion.
