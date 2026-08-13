# Midpoint discrete-parity interpolation guard

Date: 2026-08-13

Status: countermodel gate; not a proof of the half-domain completion

For each integer roster point `u`, `alpha=A+2u` is odd, so the direct
equation-(9) values at `x=1/2` do simplify:

```text
exp(i*pi*alpha^2/8)=exp(i*pi/8),
sum_(u=0)^L f_(1/2)(u)=exp(i*pi/8)(L+1)(A+B)/2.        (MG1)
```

That identity is valid and remains useful for checking the direct source
sum.  It does **not** imply that the canonical finite-Poisson interpolation
is linear.  Between integer points,

```text
f_(1/2)(u)=(A+2u)exp(i*pi(A+2u)^2/8),

f_(1/2)'(u)=exp(i*pi(A+2u)^2/8)
 [2+i*pi(A+2u)^2/2].                                   (MG2)
```

The quadratic chirp remains present for real noninteger `u`.  Consequently,
for nonzero integer `m`, the actual coefficient is

```text
I_m(1/2)=4m(-1)^m
 [F((B-4m)/2)-F((A-4m)/2)],                            (MG3)
```

where `F'(q)=exp(i*pi*q^2/2)`.  The explicit endpoint exponential in the
integration-by-parts formula vanishes because `A` and `B` have the same odd
character, but the Fresnel difference in (MG3) does not.

Actual-source interval witnesses disprove `I_m+I_-m=0`.  At `m=1`,

```text
|I_1+I_-1|=[8.007755091183675000469777162977408872614733174509926519687432367986794667849496888084563413760052182e-10 +/- 5.69e-89],
```

and at the fold edge `m=39894`,

```text
|I_39894+I_-39894|=[69468.70954509621284556606383375950775555445238184493883641109445892764822646920440910543057457938021 +/- 1.77e-84]>69000. (MG4)
```

The large midpoint value in (MG4) is not itself an error bound; the outer
`x` integral remains oscillatory.  Its role is to block an invalid argument:
discrete odd-square parity cannot be pushed through the continuous Poisson
interpolation coefficient.  Any half-domain completion must retain the
canonical chirp and use its phase, not replace it by a linear trapezoidal
interpolant.

Pi provenance: all `pi` factors come from the equation-(9) Kummer phase and
Fourier character.  No fitted constant is used.

Proof boundary: direct midpoint roster identity and explicit counterexamples
to coefficient-wise pair cancellation only.  No interior completion bound or
nonlinear B crossing theorem is proved.  No theorem establishing `T_upper`,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion follows here.
