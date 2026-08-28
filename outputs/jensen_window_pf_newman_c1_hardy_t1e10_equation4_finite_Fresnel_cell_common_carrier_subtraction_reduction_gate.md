# Finite Fresnel cells and common-carrier subtraction

Date: 2026-08-27

Status: exact cell-amplitude reduction certified; joined amplitude bound open

From Section 11.486, write

```text
D_W(y)=sum_(j=0)^(M-1)exp[i*pi*(A+2j)y],
B_W=integral_0^infinity[Fresnel]
    y^(-s)exp(-i*pi*y^2)D_W(y)dy.                    (CC1)
```

For every integer `m` and `-1/2<=u<=1/2`, oddness of `A` gives

```text
D_W(m+u)=(-1)^m D_W(u),
exp[-i*pi*(m+u)^2]D_W(m+u)
 =exp(-i*pi*u^2-2*pi*i*m*u)D_W(u).                  (CC2)
```

Define the central half-cell and positive mode amplitudes

```text
B_0=integral_0^(1/2)y^(-s)exp(-i*pi*y^2)D_W(y)dy,

beta_m(t)=integral_(-1/2)^(1/2)
 (1+u/m)^(-s)exp(-i*pi*u^2-2*pi*i*m*u)D_W(u)du,
                                                               (CC3)
B_m=m^(-s)beta_m(t),       m>=1.                     (CC4)
```

Contiguous partition of the Fresnel/Abel limit, without exchanging an
infinite series with an integral, proves

```text
B_W=B_0+sum_(m=1)^infinity[Fresnel]m^(-s)beta_m(t).  (CC5)
```

The altered finite-roster quadrature witnesses have maximum covariance and
cell-factorization discrepancies `1.146333401319095635377025e-84`
and `2.010764683385948796148028e-87`.

The exact global-Morse Gamma carrier has

```text
G_m=2 Re[B(t)exp(i*theta_0)m^(-s)],
B(t)=I_bulk(t)/I_G(t).                               (CC6)
```

Put

```text
C_G(t)=H(t)exp[-i*theta(t)]B(t)exp(i*theta_0).       (CC7)
```

Then, exactly and mode by mode,

```text
Hardy_t[C_G(t)m^(-s)]=H(t)G_m.                       (CC8)
```

Thus the ordinary owned cells are governed by the explicit coefficient
defect `beta_m-C_G`; this is subtraction before norms, not comparison of two
separately accumulated large values.

For `39853<=m<=39936`, let `a_m` denote the already-defined real physical
projection of `A_m+A_-m` and use the canonical Hardy lift

```text
mathcal A_m=(H(t)/2)exp[-i*theta(t)]a_m,
Hardy_t[mathcal A_m]=H(t)a_m.                        (CC9)
```

Because `A_transition=A_endpoint+G_extra`, the exact joined target is

```text
H(t)J_Z=Hardy_t{
 P_W+B_0
 +sum_(m=1)^621 m^(-s)beta_m
 +sum_(m=622)^39852 m^(-s)[beta_m-C_G]
 +sum_(m=39853)^39936
       [m^(-s)(beta_m-C_G)-mathcal A_m]
 +sum_(m=39937)^infinity[Fresnel]m^(-s)beta_m}.     (CC10)
```

Here

```text
P_W=-sum_(n=79789)^2561211 n^(-s)
    +C_79788-C_2561211.                              (CC11)
```

Equation (CC10) is the requested amplitude-level common lattice.  It does not
say that matching saddles makes `beta_m-C_G` or the A-corrected coefficient
small.  The next obligation is twofold: derive the natural integral lift of
`mathcal A_m` on the same cells, then join `P_W` to the lower and upper
unowned-cell packages through an exact Mordell/contour transformation before
attempting an interval norm.

The mode ownership is now algebraically disjoint:

```text
1..621:         unowned lower cells,
622..39852:     39231 ordinary Gamma-subtracted cells,
39853..39936:   84 Gamma-plus-A-subtracted cells (42+42),
39937..infinity:unowned upper cells.                 (CC12)
```

Pi provenance: every `pi` in (CC1)--(CC12) descends from the finite RSI
Gaussian, odd Fourier roster, Riemann-Siegel phase, or the exact global-Morse
Gamma normalization.  No fitted constant is introduced.

Proof boundary: exact finite-source chirp covariance, contiguous Fresnel-cell
partition, Gamma common-coefficient lift, A canonical Hardy lift, and joined
`H J_Z` reduction only.  No bound for any coefficient defect, natural
cell-integral A lift, Mordell join of `P_W` with unowned cells, actual-height
`J_Z` or `D_K` enclosure, non-A, all-height, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
