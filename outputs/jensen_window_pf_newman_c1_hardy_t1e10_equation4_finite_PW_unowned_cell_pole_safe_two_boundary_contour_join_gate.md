# Pole-safe two-boundary contour join for the unowned cells

Date: 2026-08-27

Status: exact paired-contour join certified; quantitative Mordell compression open

Put `s=1/2+it`, use the principal logarithm, and define

```text
D_W(z)=sum_(j=0)^(M-1) exp[i*pi*(A+2j)z],
g_W(z)=z^(-s)exp(-i*pi*z^2)D_W(z),
A=159577, M=2481423.                                (UJ1)
```

The finite RSI prefix and the complete branch correction are

```text
P_W=-integral_C g_W(z)dz,
B_W=integral_0^infinity[Fresnel] g_W(y)dy,
S_W=P_W+B_W.                                        (UJ2)
```

For `L=621.5` and `U=39936.5`, the lower and upper unowned cells therefore
join before any norm as

```text
U_unowned
 =P_W+integral_0^L g_W+integral_U^infinity[Fresnel]g_W
 =S_W-integral_L^U g_W.                             (UJ3)
```

Deforming `C` to the upper boundary of the real axis is legal for the finite
source.  The indentation at zero vanishes because `Re(s)=1/2`.  Hence

```text
U_unowned
 =-integral_(-infinity)^0[upper] g_W(z)dz
  -integral_L^U g_W(y)dy.                           (UJ4)
```

Writing `z=-y+i0` in the first term and rotating
`y=exp(-i*pi/4)x` through the decaying quadrant gives exactly

```text
-integral_(-infinity)^0[upper]g_W(z)dz
 =K_0 integral_0^infinity x^(-s)exp(-pi*x^2)G_W(x)dx
 =S_W,

K_0=exp(3*pi*t/4+3*pi*i/8),
G_W(x)=sum_(j=0)^(M-1)exp[-pi(1+i)(A+2j)x/sqrt(2)]. (UJ5)
```

Thus (UJ3)--(UJ5) are the same identity in cell, real-boundary, and decaying
half-line normalizations.

The finite source also has the exact csc form

```text
D_W(z)=(i/2)[exp(2*pi*i*n_-*z)-exp(2*pi*i*n_+*z)]csc(pi*z),
n_-=79788, n_+=2561211.                             (UJ6)
```

If

```text
G_N(z)=(1/(2i))z^(-s)exp(-i*pi*z^2+2*pi*i*N*z)csc(pi*z),
```

then `g_W=G_(n_+)-G_(n_-)`.  This difference must remain inside every real-
boundary integral.  At each nonzero integer `k`, both `G_N` terms have the
same residue `k^(-s)/(2*pi*i)`, so the poles cancel.  At zero each separate
term behaves as `z^(-s-1)/(2*pi*i)`, whereas the paired difference behaves
as `M z^(-s)` and is locally integrable.  Consequently separate objects
formed directly from the unregularized `G_N` are not legitimate ordinary
real-tail integrals.

There is, however, a canonical common subtraction.  Put

```text
R_N(z)=G_N(z)-G_0(z)
      =z^(-s)exp(-i*pi*z^2)
       sum_(r=0)^(N-1)exp[i*pi*(2r+1)z].             (UJ7)
```

Every `R_N` is a finite-source packet: its integer poles are removable and it
behaves as `N z^(-s)` at zero.  Therefore the individual Fresnel objects

```text
mathcal J_N(L,U)
 =integral_C R_N(z)dz
  -integral_0^L R_N(y)dy
  -integral_U^infinity[Fresnel]R_N(y)dy
```

are legitimate, and the common subtraction cancels in their difference:

```text
U_unowned=mathcal J_(n_-)(L,U)-mathcal J_(n_+)(L,U). (UJ8)
```

If `S_N` is the finite physical Kummer prefix and `B_N[L,U]` is the
corresponding finite-source interior Fresnel integral, then the same
definitions give the endpoint-complete identity

```text
mathcal J_N(L,U)=B_N[L,U]-S_N.                       (UJ8a)
```

This is the pole-safe finite/infinite contour normalization required for a
Mordell transform.  It is exact, but it is not yet a compressed or bounded
evaluation of either `mathcal J_N`.

Because `A=1 mod 4`, `M` is odd, and both production boundaries are half
integers,

```text
D_W(k+1/2)=i(-1)^k,

D_W^(r)(k+1/2)
 =i(-1)^k(i*pi)^r 2^(r-1)
   [E_r(A/2)+E_r(A/2+M)].                           (UJ9)
```

Here `E_r` is the Euler polynomial; the plus sign uses odd `M`.  This gives
all source endpoint jets in constant time for the next endpoint-complete
transformation.  If `f(z)=z^(-s)exp(-i*pi*z^2)`, then

```text
g_W^(r)(h)=sum_(q=0)^r binom(r,q)f^(r-q)(h)D_W^(q)(h),
                                      h in {L,U},  (UJ10)

(log f)' =-s/z-2i*pi*z,
(log f)''= s/z^2-2i*pi,
(log f)^(q)=(-1)^q s(q-1)!/z^q,       q>=3.
```

Together with (UJ9), this is a constant-work recurrence for every full
integrand endpoint jet.  The altered `5`-label replay at
`t=3.25` compares the hypergeometric source, direct RSI contour,
rotated negative boundary, branch correction, and finite interior integral.
Its largest join discrepancy is
`7.1748930760982397e-58`.

There is also an exact cancellation guard for the final assembly.  With
`a=39852.5`, write

```text
O_join= integral_L^a g_W
        -C_G sum_(m=622)^39852 m^(-s),

T_join= integral_a^U g_W
        -C_G sum_(m=39853)^39936 m^(-s)-mathcal A_A^nat.
```

Then (UJ3) gives, before projection,

```text
U_unowned+O_join+T_join
 =S_W-C_G sum_(m=622)^39936 m^(-s)-mathcal A_A^nat. (UJ11)
```

The two positive Fresnel subintervals cancel the `-integral_L^U g_W` in
`U_unowned` exactly.  Therefore the three packets must be assembled before a
final norm; (UJ11) is the cancellation-preserving complete target, not three
unrelated absolute-value estimates.

Pi provenance: every `pi` in (UJ1)--(UJ11) descends from the original RSI
Gaussian/sine kernel, the odd Fourier roster, the fixed quarter-turn, or the
Riemann-Siegel branch phase.  No fitted circle constant is inserted.

Proof boundary: (UJ3)--(UJ11) certify the exact pole-safe contour join and
remove the separate infinite upper-tail obligation.  They do not compress or
bound the two-boundary packet, bound the ordinary Gamma-subtracted packet,
enclose `J_Z` or `D_K`, prove a non-A bound or all-height theorem, or prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
