# Symmetric endpoint-tail reassembly after bulk extraction

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the quantitative residual bound

Write the odd source roster as `alpha=A+2u`, where

```text
A=159577, B=5122421, L=2481422.
```

The portcullis coefficient is exactly the ordinary Fourier coefficient:

```text
J_m=(-1)^m integral_0^L f_x(u)exp[-i*pi*m(A+2u)]du
   =integral_0^L f_x(u)exp(-2*pi*i*m*u)du=I_m.       (SR1)
```

The last equality uses only that `A` is odd.  For every finite symmetric
cutoff `M>=39894`, define

```text
D_M(u)=sum_(m=-M)^M exp(-2*pi*i*m*u)
      =sin((2M+1)pi*u)/sin(pi*u),

G_T(u)=sum_(m=622)^39894 exp(-2*pi*i*m*u)
      =exp(-40516*pi*i*u)sin(39273*pi*u)/sin(pi*u),

C_M(u)=D_M(u)-G_T(u).                                (SR2)
```

The quotients use their removable values at integer `u`.  In particular,
`D_M(k)=2M+1`, `G_T(k)=39273`, and `C_M(k)=2M+1-39273`.

Let `H_x=[f_x(0)+f_x(L)]/2`, let `B_m` be the exact full-line bulk current,
and let `Q_m=J_m-B_m=P_A,m+P_B,m` be the grouped finite-endpoint tail.  Pure
finite algebra gives the cancellation-preserving residual

```text
S_M-sum_(m=622)^39894 B_m
 =H_x+integral_0^L f_x(u)C_M(u)du
     +sum_(m=622)^39894 Q_m,                          (SR3)

S_M=H_x+integral_0^L f_x(u)D_M(u)du.
```

Equation (SR3) is the exact object left after Section 11.349 closes the
Gamma bulk.  The previously proved symmetric-Poisson theorem gives the
limit of the `D_M` integral.  Since `G_T` is a fixed finite trigonometric
polynomial, it also proves the symmetric limit of the right side of (SR3).
No auxiliary Abel regulator is needed for this fixed finite roster.

The mode ledger identifies the cancellation precisely.  Negative modes
`-39894..-622` pair with the target positive modes; `-621..-1` pair with the
low positive completion; and `-M..-39895` pair with the outer positive
completion.  For every pair,

```text
I_m+I_-m=[-2 Delta f_x'+R_m+R_-m]/(2*pi*i*m)^2,       (SR4)
```

so the apparent `Delta f_x/(2*pi*i*m)` current vanishes before absolute
values.  The zero mode has no `1/m` current.  The endpoint half-current is
the finite-roster midpoint correction, not its canceller.

This also sharpens the guard from Section 11.350.  The completed-square
bare endpoint exponentials do not cancel directly: after odd parity they
are mode-independent.  Their compensating terms live in the Fresnel tails.
Only the grouped `J_m`, followed by symmetric `+m/-m` pairing, exposes the
true cancellation.  Bounding `P_A`, `P_B`, the negative modes, or the
endpoint half-current separately would destroy (SR3)--(SR4).

The next quantitative object is therefore one joint residual, not several
endpoint-error sums:

```text
R_end=lim_(M->infinity) [H + integral f C_M
                         +sum_target(P_A+P_B)].        (SR5)
```

The upper `B` part should first be enclosed nonstationarily in this grouped
form; the lower `A` characteristic part can then be spliced to the certified
fold atlas on the same current.

Pi provenance: every `pi` in (SR1)--(SR5) comes from the equation-(9)
Kummer quadratic phase or the integer Fourier-Poisson character.  No circle,
polygon, fitted period, or inserted geometric constant is used.

Proof boundary: exact mode partition, Dirichlet-kernel representation,
symmetric residual limit, and pair-current cancellation only.  No numerical
bound for (SR5), complete endpoint-tail theorem, ordinary/fold splice,
complete `T_upper`, height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
