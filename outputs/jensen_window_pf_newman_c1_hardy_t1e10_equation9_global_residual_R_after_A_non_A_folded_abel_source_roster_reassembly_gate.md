# Folded Abel source-roster reassembly

Date: 2026-08-26

Status: exact source/ownership reassembly and rational transition schedule;
the physical non-A bound remains open.

## Source reconstruction

Write

```text
f_x(n)=(A+2n)exp(i*pi*x*(A+2n)^2/4),       0<=n<=L,
A=159577, B=5122421, L=2481422.                                      (SR1)
```

The inherited endpoint half-current and the Abel boundary layer are

```text
H_x=[f_x(0)+f_x(L)]/2,
F_x(0)=sum_(n=0)^(L-1)f_x(n),
F_x(1)=sum_(n=1)^L f_x(n).                              (SR2)
```

Therefore, coefficient by coefficient,

```text
H_x+[F_x(0)+F_x(1)]/2=sum_(n=0)^L f_x(n).              (SR3)
```

This is the complete odd source roster `A,A+2,...,B`, containing
`2481423` labels.  At `x=0`, (SR3) is the exact integer

```text
6553435661577.             (SR4)
```

No split `1/x` formula is needed at that endpoint.

## Ownership reduction

The Abel-limit source contributes `-I_m` on the extended notch.  On
`622..39852`, the endpoint row is `A_m+B_m=I_m-P_m`; on
`39853..39936`, it is
`B_m-A_-m=I_m-P_m-A_m-A_-m`.  Combining each row before a norm gives

```text
H_x+AbelSource+EndpointRows
 =sum_(n=0)^L f_x(n)
  -sum_(m=622)^39936 P_m(x)
  -sum_(m=39853)^39936 [A_m(x)+A_-m(x)].              (SR5)
```

Thus the 39,315 raw finite Fresnel coefficients are a valid independent
cross-check, but they are not the primary physical quadrature: every one of
their owned contributions cancels algebraically to the already defined bulk
and A-face atoms.  The B-trace, B-outer, and translated A-face subtractions
remain unchanged in the common non-A limit.

The post-A positive notch has the exact route partition

```text
622..39694:   39073 ordinary-core modes,
39695..39852:   158 remaining fold-owned interior modes,
39853..39936:    84 extracted A-transition modes.      (SR6)
```

This partition selects evaluator coordinates; it does not license separate
absolute bounds.

## Transition schedule

For `phi_m(y)=x*y^2/4-m*y`, the finite coefficient has stationary point
`y_*=2m/x`.  It enters at `x_B(m)=2m/B` and exits at `x_A(m)=2m/A`.
On `0<=x<=1/2` there are `78588` distinct
rational events.  The closest pair is

```text
x_A(630) < x_B(20223),
x_B(20223)-x_A(630)=882/817420575917
                    =1.0790039130230523e-09... .             (SR7)
```

All modes `39895..39936` have `x_A(m)>1/2`, so their finite-`y` stationary
points persist through the midpoint.  Their joint physical saddles lie in
the reflected `x>1/2` half; exact Kummer reflection supplies the full source
without introducing a second half-domain stationary family.  After quadratic
completion, the common physical phase and saddle are

```text
Psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x),
x_m=2*pi*m^2/(t+2*pi*m^2),
alpha_m=2m+t/(pi*m),                                   (SR8)
```

with `x_B<=x_m<=x_A` exactly when `A<=alpha_m<=B`.  At the saved height this
reproduces 39,231 interior modes, 42 lower A-transition modes, and 42
reflected A-transition modes.

## Route decision

Use the exact double-Weber source reconstruction, the closed full-line bulk,
the remaining 158 fold modes, and the already owned A/B channels in one
signed physical assembly.  Use the endpoint-plus-39,315-Fresnel form only as
an altered-representation check at sparse points.  The next numerical gate is
a bounded source-minus-owned-carriers pilot split at the rational events in
(SR7), followed by interval `x` quadrature if the signed scale is viable.

## Pi provenance

Every `pi` in (SR1)--(SR8) is inherited from the Gaussian Abel/Fourier
kernel, the original Kummer quadratic phase, or the physical Kummer weight.
The rational event locations themselves contain no `pi`; the `pi` in the
physical saddle is forced by differentiating the declared phase.

## Boundary

The exact Abel-zero reassembly of the inherited half-current, folded endpoint boundary layer, and two post-A endpoint rows into the complete 2481423-label odd source minus the contiguous P block and 84 paired A atoms; the exact x=0 source value; and the complete rational finite-phase transition schedule on 0<=x<=1/2 only. The route partition is not permission for independent norms. No numerical signed assembly, fast-evaluator error theorem, physical x quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
