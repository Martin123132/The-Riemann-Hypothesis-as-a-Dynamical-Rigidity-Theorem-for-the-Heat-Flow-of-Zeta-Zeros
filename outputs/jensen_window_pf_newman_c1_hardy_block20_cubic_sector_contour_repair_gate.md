# Hardy block-20 cubic sector-contour repair gate

Date: 2026-08-07

Status: exact sign-aware contour lemma on the saved cubic roster; not a proof of the full recurrence or RH

## Sector choice

Let `F(z)=Phi1*z+Phi2*z^2+Phi3*z^3`, with `F''>0` on `[0,N]`.
For `(67b)` put `g_n(z)=n*z-F(z)` and for `(67c)` put
`h_n(z)=n*z+F(z)`.  The sign-aware rays are

```text
                         Phi3 < 0       Phi3 > 0
(67b), g_n                 5*pi/6          pi/2
(67c), h_n                   pi/2          pi/6
```

For the two slanted cases, exact cubic Taylor expansion at any
`x in [0,N]` gives

```text
Im[g_n(x+r*exp(5*pi*i/6))-g_n(x)]
 = (n-F'(x))*r/2 + sqrt(3)*F''(x)*r^2/4 + |Phi3|*r^3,  Phi3<0,

Im[h_n(x+r*exp(pi*i/6))-h_n(x)]
 = (n+F'(x))*r/2 + sqrt(3)*F''(x)*r^2/4 + Phi3*r^3,    Phi3>0.
```

On the retained vertical rays the corresponding identities are
`(n-F'(x))*r+Phi3*r^3` for positive `Phi3` in `(67b)` and
`(n+F'(x))*r+|Phi3|*r^3` for negative `Phi3` in `(67c)`.  Every displayed
coefficient is positive.  Thus the translated far connector is bounded by
a quadratic polynomial in `r` times `exp(-2*pi*|Phi3|*r^3)` and tends to
zero.

## Exact roster checks

All 374 recursive calls satisfy the hypotheses:

```text
strict F''>0 calls                         374
negative / positive Phi3 calls            187 / 187
67b slanted / vertical contracts          187 / 187
67c slanted / vertical contracts          187 / 187
minimum F'' on [0,N]                      >= 1.18510704665172836425957259968383850220841321910285E-2
minimum 67b linear decay coefficient      >= 7.31783587828757338558964516060058600943803949358837E-4
minimum 67c linear decay coefficient      >= 2.54841047223187318260667510347044972172920215070425E-1
minimum cubic decay coefficient           >= 1.09815274373050800131128053091714181465651610988342E-8
minimum slanted quadratic coefficient     >= 5.13166404302173312625030750109516841569450549694352E-3
```

The 90- and 150-digit `sqrt(3)` interval evaluations overlap on every row.
Since the integrands are entire, the finite parallelogram deformation crosses
no singularity and creates no residue.  A saddle is not a singularity; Stokes
switching concerns a later asymptotic decomposition, not this exact contour
identity.

## Quadratic rotation and summation

The endpoint quadratic model is compatible with the paper's vertical
error-function contour.  For `(67b)`, rotate from `5*pi/6` to `pi/2`:
`sin(theta)>=1/2` and `-sin(2*theta)>=0`.  For `(67c)`, rotate from `pi/6`
to `pi/2`: `sin(theta)>=1/2` and `sin(2*theta)>=0`.  Hence the large
quadratic-model arc decays and the ray rotation is exact.

The outer `1/n` sums may also be interchanged with the repaired endpoint
rays on this finite roster.  After taking absolute values, their common
kernel is bounded near `r=0` by a constant times
`1+|log(r)|`, because `sum(exp(-pi*n*r)/n)=-log(1-exp(-pi*r))`; at infinity
the certified cubic factor dominates every polynomial amplitude.  The
logarithmic singularity is integrable.

## Pi provenance and remaining boundary

Here `pi` comes from the original Fourier exponential
`exp(2*pi*i*phase)`.  The angles are selected by the exact trigonometric
values at `pi/6`, `pi/2`, and `5*pi/6`; no circle-derived numerical fit is
inserted.

This repairs the existence of the exact nonsaddle contours and preserves the
paper's quadratic endpoint-ray value.  It does not bound the difference
between the full cubic endpoint integrals and their quadratic models, justify
the other amplitude truncations in W2--W4, establish a height-uniform
recurrence, control the outer Hardy remainder, or prove `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
