# Hardy block-20 cubic vertical-contour admissibility gate

Date: 2026-08-07

Status: exact cubic far-field obstruction to the stated vertical nonsaddle contours; repair target, not a proof of RH

## Exact far-field calculation

Write

```text
F(z)=Phi1*z+Phi2*z^2+Phi3*z^3,   xi=F'(N).
```

For equation (82), the (67b) phase is `g_n(z)=n*z-F(z)`.  On the
paper's upper horizontal/vertical geometry,

```text
Im g_n(x+iR)=R*(n-F'(x))+Phi3*R^3.                    (1)
```

The bound asserted below equation (82) drops the cubic term.  If `Phi3<0`,
then at `x=N`, `n=ceil(xi)`, (1) becomes negative once
`R^2>(ceil(xi)-xi)/|Phi3|`; consequently
`|exp(2*pi*i*g_n)|=exp(-2*pi*Im g_n)` grows super-exponentially and the
vertical integrand in equation (83) does not even tend to zero.

For (67c), the phase in equation (90) is `h_n(z)=n*z+F(z)`, and

```text
Im h_n(x+iR)=R*(n+F'(x))-Phi3*R^3.                    (2)
```

If `Phi3>0`, then at `x=0`, `n=1`, (2) becomes negative once
`R^2>(1+Phi1)/Phi3`.  The stated vertical integral again cannot converge as
an ordinary improper integral.

## Exact roster

All 374 recursive block-20 cubics have nonzero `Phi3`:

```text
Phi3 < 0: 187 calls; stated (67b) vertical contour fails
Phi3 > 0: 187 calls; stated (67c) vertical contour fails
Phi3 = 0: 0 calls
minimum crossover radius  >= 3.04017519589532305257245793795369011624623853806167E+2
maximum crossover radius  <= 1.14527362847685623456996243710184836037245870626458E+4
```

At the exact witness `R=2*R_cross`, both cases satisfy
`Im phase=-3*linear*R<0`.  The two independent 90/150-digit evaluations
overlap on every call, and the maximum identity gap is
`1.62993238224986736114733121483479073702073509112752E-145`.

This does not show that the local W2--W4 formulas are unusable.  It shows
that their printed infinite vertical-contour derivation is incomplete for
every nonzero cubic on this roster.  In particular, the missing operation is
a contour deformation with a controlled connector and tail, not another
free fitted endpoint term.

## Repair target

The ray `z=z0+r*exp(i*pi/6)` is a valid cubic-only decay candidate for both
failing signs because `sin(3*pi/6)=1`.  Requiring the retained quadratic
endpoint model to decay as well sharpens the contract: use `5*pi/6` for the
negative-`Phi3` `(67b)` branch and `pi/6` for the positive-`Phi3` `(67c)`
branch, while the opposite-sign families retain `pi/2`.  The companion
sector-contour gate proves the resulting connector decay and summability;
the full-cubic versus quadratic endpoint remainder remains the next bound.

## Pi provenance and boundary

Pi appears only in the Fourier modulus identity
`|exp(2*pi*i*w)|=exp(-2*pi*Im(w))`; its positivity scales the growth or
decay and does not determine the cubic sign.  This gate identifies a missing
contour justification.  It does not yet supply the repaired contour, bound
the exact (67b)--(67c) aggregate, prove a height-uniform recurrence, control
the outer Hardy representation, or prove `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
