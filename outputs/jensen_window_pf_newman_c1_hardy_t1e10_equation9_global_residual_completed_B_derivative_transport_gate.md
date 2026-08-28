# Exact derivative transport for the crossing-completed B current

Date: 2026-08-13

Status: exact derivative reduction; not a proof artifact

No quantitative derivative norm or grouped-remainder bound is asserted.

After odd-endpoint parity, factor the common B phase.  For `s=+1,-1` put

```text
a=m sqrt(2/x),                 b=B sqrt(x/2),
a_s=s a,                       q_s=b-s a,
rho=1/(i*pi),

Y_+^tau=e^(-i*pi*q_+^2/2)[P_B,++(1-tau_m)P_bulk,+],
Y_-    =e^(-i*pi*q_-^2/2)P_B,-.                       (DT1)
```

The added non-target bulk term in `Y_+^tau` is precisely the crossing
completion from Section 11.417.  Let

```text
K_s=1/(2x)+i*pi*q_s q_s'.                             (DT2)
```

Using only `F'(q)=exp(i*pi*q^2/2)`, the phase-stripped endpoint tail obeys
the closed transport equation

```text
Y_s'=a_s q_s'+(rho-Y_s)K_s,                           (DT3)

q_s' =(b+s a)/(2x),
q_s''=(-b-3s a)/(4x^2),
K_s' =-1/(2x^2)+i*pi[(q_s')^2+q_s q_s''].             (DT4)
```

The bulk completion is not an extra forcing term.  Its phase-stripped value

```text
V_+=(1-tau_m)a(1+i)e^(-i*pi*q_+^2/2)
```

satisfies exactly `V_+'=-K_+V_+`; therefore the completed `Y_+^tau` obeys
the same inhomogeneous equation (DT3).  Differentiating once more gives

```text
Y_s''=a_s(q_s''-q_s'/(2x))-Y_s'K_s+(rho-Y_s)K_s'.     (DT5)
```

The complete phase-stripped mode pair is

```text
C_m=(Y_+^tau+Y_-)/x,
C_m'=(Y_+'+Y_-')/x-(Y_++Y_-)/x^2,
C_m''=(Y_+''+Y_-'')/x-2(Y_+'+Y_-')/x^2
       +2(Y_++Y_-)/x^3.                               (DT6)
```

Finally, for `w(x)=[x(1-x)]^(-1/4)`,

```text
w'/w =(2x-1)/[4x(1-x)],
w''/w=(12x^2-12x+5)/[16x^2(1-x)^2].                  (DT7)
```

Thus, at every common finite cutoff and regulator,

```text
A_(M,epsilon)=w sum_(m=1)^M w_(m,epsilon)C_m
```

has explicit first and second derivatives obtained from (DT3)--(DT7).
The normalized derivatives are `A_xi=A_x/sqrt(H_B)` and
`A_xixi=A_xx/H_B`.  No numerical differencing of Fresnel tails and no
distributional crossing term is required.

The remaining theorem is quantitative: obtain sum-first bounds for these
finite derivative expressions that are uniform in `M` and `epsilon`, include
the two B-window edge terms, and combine the tangential estimate with the
common-regulator joined remainder below `1.4058e-4`.

Pi provenance: `pi` in (DT1)--(DT5) is inherited from the exact equation-(9)
Fresnel phase.  Equations (DT6)--(DT7) are calculus identities and introduce
no fitted geometric constant.

Proof boundary: exact phase-stripped first/second derivative transport and
finite-sum differentiation only.  This gate does not prove uniform derivative
bounds, a completed exterior-current estimate, a joined-remainder estimate,
complete `Q_K-T` or `T_upper`, or a height-uniform theorem.  It makes no claim
of `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
