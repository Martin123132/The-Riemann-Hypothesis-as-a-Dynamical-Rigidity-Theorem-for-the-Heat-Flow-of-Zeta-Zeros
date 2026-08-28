# Finite-regulator equivalence for `R_after_A`

Date: 2026-08-23

Status: exact common-kernel/mode-mask equivalence and physical limit order
certified; no quantitative bound for `R_after_A`

Let `T={622,...,39894}`, `w_m=exp(-pi*epsilon*m^2)`, and take one
finite symmetric cutoff

```text
M>=B=5122421,       epsilon>0.                              (RA1)
```

The stronger cutoff in (RA1) is deliberate.  The A deletion itself only
needs `M>=39936`, but the finite analytic B exterior
`O_(M,epsilon)=eta_delta E_W Ctr_>=B,M,epsilon` belongs to the same identity
only for `M>=B`.

The common Gamma-normalized projector kernel is

```text
Delta_(M,epsilon)(x)
 =H_x+integral_0^2481422 f_x(u)[D_(M,epsilon)-G_(T,epsilon)]du
   +sum_(m in T)w_m(A_m+B_m).                         (RA2)
```

Delete the already-certified A transition before taking a norm:

```text
mathcal A_(M,epsilon)
 =sum_(m=39853)^39936 w_m
   [A_m+A_-m+(1-chi_T(m))P_m].                        (RA3)

mathfrak R_A,(M,epsilon)
 =Delta_(M,epsilon)-chi_W Btr_(M,epsilon)
   -O_(M,epsilon)-mathcal A_(M,epsilon).               (RA4)
```

Equations (RA2)--(RA4) are executable without assigning independent norms
to any cancellation partner.  Expanding the common Fourier kernel, pairing
`+m` and `-m`, and using `I_m=P_m+A_m+B_m` gives the exact coefficient mask

```text
mode range       post-A joined summand
1..621           P_m+P_-m+A_m+A_-m+B_m+B_-m
622..39852       P_-m+A_m+A_-m+B_m+B_-m
39853..39894     P_-m+B_m+B_-m
39895..39936     P_-m+B_m+B_-m
39937..M         P_m+P_-m+A_m+A_-m+B_m+B_-m.          (RA5)
```

Together with `H_x+I_0`, (RA5), `-chi_W Btr`, and `-O` is exactly (RA4),
not an approximation.  The two middle rows have the same surviving vector
for different reasons: the target projector had already removed `P_m` in
the first row, while the A transition removes it in the second.  No
`P_-m` or B endpoint atom is deleted.

With the inherited equation-(9) physical functional

```text
P_t[F]=2(pi/(32t))^(1/4)
       Re[e^(-i*pi/8) integral_0^1 W_t(x)F(x)dx],     (RA6)
```

the object being targeted is precisely

```text
R_after_A
 =lim_(epsilon down 0) lim_(M to infinity)
    P_t[mathfrak R_A,(M,epsilon)],                    (RA7)

R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A.    (RA8)
```

The `M` limit is taken first at fixed positive `epsilon`; only then does the
common Gaussian Abel regulator tend to zero.  The separately certified B
outer block keeps its dominated Abel passage.  The endpoint half-current,
zero mode, negative bulk, remaining endpoint atoms, and remote-positive
pairs remain joined throughout.

Pi provenance: every `pi` in (RA1)--(RA8) comes from the original Kummer
quadratic phase, integer Fourier character, Gaussian Abel regulator, and
the inherited physical normalization.  No circular construction, fitted
constant, or new geometric occurrence of `pi` is introduced.

Proof boundary: exact finite-regulator representation, coefficient-mask
equivalence, cutoff compatibility, and common physical limit order at
`t=10^10` only.  No numerical or analytic upper bound for `R_after_A`,
complete `R_Dir` or `Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
