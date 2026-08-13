# Newman Arbitrary-Multiplicity Tertiary Universal Layer Gate

Date: 2026-08-03

Status: exact universal `D^(-1/3)` correction and multiplicity-three collision refinement; full first Xi-jet order open; not a proof of RH.

## Extended Scaling

Set epsilon=D^(-1/6), lambda=lambda_*+rho*epsilon^6/8+(rho*q/4)*epsilon^8+(rho*r/4)*epsilon^10, and evaluate at (rho+rho*y*epsilon^4)/D.

After normalization by rho^m U(rho) epsilon^(4m), the finite Jensen cluster is P_m(y,q)+epsilon^2 Q_(m,n)(y,q,r)+O(epsilon^4) locally, where P_m=exp(q*partial_y^2+partial_y^3/12)y^m.

The complete correction is

```text
Q_(m,n)=((r-y/2)*partial_y^2+(2n+1)*partial_y/4)P_m.
```

All apparent fourth-order finite-Jensen and second heat-commutator terms cancel or recombine into the displayed first-order transport of P_m; no independent partial_y^4 term survives at epsilon^2.

The exact audit checks this identity for 14 symbolic multiplicities and 12 finite radial-heat/Jensen models at shifts zero and three.

## Shift And Appell Structure

The fixed Jensen shift first enters at epsilon^2 through (2n+1)partial_y P_m/4.

The universal family obeys partial_q P_m=partial_y^2 P_m and partial_y P_(m+1)=(m+1)P_m.

## Contact Correction

At any nondegenerate threshold double root P_m=P_m'=0, P_m''!=0, the next parameter and root coefficients are r_m=y_m/2 and s_m=(1-2n)/4.

For multiplicity three,

```text
r_3=2^(-5/3)
s_3=(1-2n)/4
lambda_D=lambda_*+rho/(8D)-rho*2^(-13/3)D^(-4/3)+rho*2^(-11/3)D^(-5/3)+o(D^(-5/3))
w_D=rho/D+rho*2^(-2/3)D^(-5/3)+rho*(1-2n)D^(-2)/4+o(D^(-2))
```

At the Newman shift `n=0`,

```text
For rho=-c^2 and n=0: t_D=Lambda-c^2/(8D)+c^2*2^(-13/3)D^(-4/3)-c^2*2^(-11/3)D^(-5/3)+o(D^(-5/3)).
```

The root-center coefficient rho*(1-2n)/4 agrees exactly with the existing double-zero finite-Jensen center correction, although the preceding multiplicity-three scales are different.

## First Xi-Specific Local Term

The local unit is absent from Q_(m,n). Its first linear contribution occurs at epsilon^4 and equals rho*[U'(rho)/U(rho)]*P_(m+1)(y,q).

This coefficient is independently visible in every exact model, but it is not the complete `epsilon^4` correction.

At epsilon^4 the Xi local unit appears together with a still-uncomputed universal background and the next heat/root parameters. Its isolated coefficient alone has no sign implication.

## Live Handoff

Compute the full `epsilon^4=D^(-2/3)` operator, including the next heat and root parameters, the universal background, and `rho U'(rho)/U(rho) P_(m+1)`. Then test the actual Xi logarithmic unit jet against exact generic-unit countermodels and obtain a degree-uniform analytic remainder.

## Pi Provenance

No pi enters this local operator calculation. All constants are rational, powers of two from the exact cubic discriminant, or the algebraic local data `rho`, `q_m`, and `y_m`.

## Proof Boundary

The local expansion remains noncontradictory and requires a degree-uniform analytic remainder before any collision exclusion can be claimed. This gate proves no Xi-specific collision exclusion, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
