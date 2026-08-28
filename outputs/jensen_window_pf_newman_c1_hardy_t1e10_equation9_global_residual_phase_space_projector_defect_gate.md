# Endpoint-completed phase-space projector defect

Date: 2026-08-13

Status: exact signed-projector reduction and transition audit; not a bound for
the compressed endpoint defect

Put `U=[0,L]`, `L=2481422`, `T={622,...,39894}`, and

```text
K_x(u,m)=f_x(u)exp(-2*pi*i*m*u),
I_m(x)=integral_U K_x(u,m)du,
P_m(x)=Abel-integral_R K_x(u,m)du.                   (PD1)
```

For `M>=39894` and the common weight `w_m=exp(-pi*epsilon*m^2)`, the
Gamma-normalized finite residual integrand is exactly

```text
Delta_(M,epsilon)
 =H_x+sum_(m=-M)^M w_m I_m
     -sum_(m=622)^39894 w_m P_m.                    (PD2)
```

Let `F_M={-M,...,M}`.  Splitting the first sum by the target projector gives

```text
Delta_(M,epsilon)
 =H_x+sum_(m in F_M minus T) w_m I_m
     -sum_(m in T) w_m(P_m-I_m),                    (PD3)

F_M minus T=[-M,621] union [39895,M].                (PD4)
```

Thus the residual is one endpoint-completed, oriented phase-space rectangle:
the finite source interval over frequencies outside `T`, minus the two
full-line spatial tails over frequencies inside `T`, with the Poisson
half-current retained.  The full-line integrals and the limit in (PD1)--(PD3)
use the already-certified common Abel prescription.  No divergent piece is
split independently.

After the inherited physical projection and common limit,

```text
lim Delta_(M,epsilon)=R_KGamma=Q_K-G
 =E_Btr,win+E_outer+R_Dir.                           (PD5)
```

Therefore (PD3) is an exact representation of the whole Gamma-normalized
residual, not of `R_Dir` alone.  The latter is obtained only after subtracting
the certified B window and outer block.

The saved-height stationary audit identifies three occupancy interfaces:

```text
B face:                 621|622,
A face:                 39852|39853,
half-boundary corner:   39894|39895. (PD6)
```

By contrast, `39694|39695` is only the certified ordinary/fold proof-ownership
boundary.  Both the A-face gap and the half-boundary gap retain the same sign
across it, so it is not an admissible place to split norms.  At the upper
corner, `4*39894=A-1` and `4*39895=A+3`; this is the lattice bracket of the
same A null face identified by the exact bi-Morse chart.

The old height-uniform `C^2` interchange majorant remains rigorous but is
quantitatively useless here:

```text
C_AB/(2*pi^2) = [25519454103834020117297802.049102364899541820144413028107520521758596699227478761 +/- 4.10e-55],
[C_AB/(2*pi^2)]/39894 = [639681508593623605486.98556296942810697202136021489517490150202432939036515462880 +/- 1.57e-60],
[C_AB/(2*pi^2)]/5122421 = [4981912674462723801.3622468846473893691170288706088445497784195712528703180544437 +/- 2.79e-62],
M needed by that majorant for the working target
  > [181530796190574657530674282496.00000000000000000000000000000000000000000000000000 +/- 3.51e+13]. (PD7)
```

This is not evidence that the actual projector defect is large.  It proves
that a raw derivative norm destroys too much cancellation.  The quantitative
theorem must estimate (PD3) as one signed object, use the exact hyperbolic
triangle coordinates away from (PD6), and install compatible B-face, A-face,
and half-corner charts before taking absolute values.

Pi provenance: `pi` in (PD1)--(PD7) comes from the equation-(9) Fourier
character, Gaussian Abel weight, Kummer phase, and the exact stationary
equations.  No fitted or geometric constant is introduced.

Proof boundary: exact finite/common-Abel projector algebra, mode counts,
saved-height transition classification, and failure of one inherited raw
majorant only.  No phase-adapted signed estimate, bound for `R_Dir`, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
