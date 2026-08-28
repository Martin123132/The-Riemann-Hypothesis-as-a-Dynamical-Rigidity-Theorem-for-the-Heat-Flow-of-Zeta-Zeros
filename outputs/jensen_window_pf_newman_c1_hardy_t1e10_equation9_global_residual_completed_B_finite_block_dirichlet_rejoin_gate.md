# Dirichlet rejoin of the finite completed-B block

Date: 2026-08-13

Status: exact common-regulator compression and revised quantitative target;
not a proof of a bound for the compressed residual

For `M>=B=5122421`, split the smooth crossing-completed current at the analytic
tail threshold:

```text
Ctr_(M,epsilon)=Ctr_<B,epsilon+Ctr_>=B,M,epsilon,
O_(M,epsilon)=eta_delta E_W Ctr_>=B,M,epsilon.        (DR1)
```

Do **not** norm the remaining `5122420` positive modes or construct their
`1,280,347` crossing charts.  Rejoin them to the common remainder:

```text
J_(M,epsilon)
 =eta_delta E_W Ctr_<B,epsilon+R_join,comp(M,epsilon).
```

The completed allocation then gives the finite identity

```text
G_(M,epsilon)=chi_W Btr_(M,epsilon)
               +O_(M,epsilon)+J_(M,epsilon),

J_(M,epsilon)=G_(M,epsilon)-chi_W Btr_(M,epsilon)
               -O_(M,epsilon).                       (DR2)
```

Thus the finite crossing roster is represented by the already-established
global kernel, not by millions of local Fresnel charts.  With
`w_m=exp(-pi*epsilon*m^2)`, put

```text
D_(M,epsilon)(u)=sum_(m=-M)^M w_m exp(-2pi*i*m*u),
G_(T,epsilon)(u)=sum_(m=622)^39894 w_m exp(-2pi*i*m*u),
C_(M,epsilon)=D_(M,epsilon)-G_(T,epsilon).            (DR3)
```

The compressed global integrand is exactly

```text
G_(M,epsilon)(x)
 =H_x+integral_0^2481422 f_x(u)C_(M,epsilon)(u)du
  +sum_(m=622)^39894 w_m(P_A,m+P_B,m).                (DR4)
```

At `epsilon=0`, (DR3) becomes the two removable Dirichlet quotients

```text
D_M=sin((2M+1)pi*u)/sin(pi*u),
G_T=exp(-40516pi*i*u)sin(39273pi*u)/sin(pi*u).        (DR5)
```

The analytic outer gate supplies its own dominated Abel limit and
`|E_outer|<8e-10`.  Consequently

```text
R_end=E_Btr,win+E_outer+R_Dir,                        (DR6)

R_Dir < 1.405792e-4
        ==> R_end < 8.6e-6,

R_Dir < 1.319792e-4
        ==> R_end < 0.                                (DR7)
```

This is the key compression: the next quantitative object is one
Dirichlet/Abel kernel plus the fixed `39273`-mode target-tail block,
with the certified B window and analytic outer block removed exactly.  The
remaining object still contains the endpoint half-current, A-fold sector,
negative modes, and nonlocal steps in their valid common grouping.

No estimate for `R_Dir` is proved here.  In particular this does not prove
the working target, a complete `T_upper`, a height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.

Pi provenance: `pi` in (DR3)--(DR5) is forced by the Gaussian Abel weight
and integer Fourier character already present in equation (9); no fitted or
geometric value is inserted.
