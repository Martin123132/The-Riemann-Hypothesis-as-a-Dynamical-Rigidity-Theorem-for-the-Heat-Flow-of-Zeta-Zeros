# Newman Arbitrary-Multiplicity First Xi-Jet Layer Gate

Date: 2026-08-03

Status: complete fixed-m relative `epsilon^4` operator and exact multiplicity-three Xi-jet displacement; actual Xi field bound and degree-uniform exclusion open; not a proof of RH.

## Complete Expansion

Set epsilon=D^(-1/6), lambda=lambda_*+rho*epsilon^6/8+(rho*q/4)*epsilon^8+(rho*r/4)*epsilon^10+(rho*v/4)*epsilon^12, and evaluate at (rho+rho*y*epsilon^4)/D.

After normalization by rho^m U(rho) epsilon^(4m), the cluster is P_m+epsilon^2 A_(n,r)P_m+epsilon^4{R_(n,q,r,v)P_m+ell P_(m+1)}+O(epsilon^6), where ell=rho U'(rho)/U(rho).

The tertiary operator is

```text
A_(n,r)=(r-y/2)*partial_y^2+(2n+1)*partial_y/4.
```

The complete universal background is

```text
R P=v*P''+(A^2 P)/2+q*(2n+1)P'/2+(q*y+(2n+1)/8)P''+(q^2+r/2)P'''+q*P''''/4+P'''''/80.
```

Conjugating the remainder by exp(-q*partial_y^2-partial_y^3/12) gives [q*(2n+1)partial_y/2+(q*y+(2n+1)/8)partial_y^2+(r/2-q^2)partial_y^3+partial_y^5/80]y^m; conjugating y back adds 2q partial_y+partial_y^2/4 and yields R.

## First Local Xi Jet

The complete source-specific contribution at relative epsilon^4 is ell P_(m+1), with ell=rho U'(rho)/U(rho); U'' first enters later.

For multiplicity three, put `ell=rho U'(rho)/U(rho)`. The next contact coefficients are

```text
v_3=-(2*ell+3*n+2)/8
t_3=2^(-7/3)*(2*ell+n+4)
lambda_D=lambda_*+rho/(8D)-rho*2^(-13/3)D^(-4/3)+rho*2^(-11/3)D^(-5/3)-rho*(2*ell+3*n+2)D^(-2)/32+o(D^(-2))
w_D=rho/D+rho*2^(-2/3)D^(-5/3)+rho*(1-2*n)D^(-2)/4+rho*2^(-7/3)*(2*ell+n+4)D^(-7/3)+o(D^(-7/3))
```

At `n=0` and `rho=-c^2`,

```text
For rho=-c^2 and n=0: t_D=Lambda-c^2/(8D)+c^2*2^(-13/3)D^(-4/3)-c^2*2^(-11/3)D^(-5/3)+c^2*(ell+1)D^(-2)/16+o(D^(-2)).
```

## Canonical-Product Field

U'(rho)/U(rho)=F^(m+1)(rho)/((m+1)F^(m)(rho))

If F(s)=F(0)prod_j(1-s/rho_j) with multiplicities mu_j, then ell_k=rho_k*sum_(j!=k)mu_j/(rho_k-rho_j).

Negative-rootedness and positive coefficients do not fix the sign of ell at an interior multiple root. The two exact witnesses have `ell=-1/3` and `ell=4/3` while retaining only negative roots and positive coefficients.

Therefore the new local term moves the collision but does not contradict it without additional Xi information.

## Exact Audit

The builder checks 28 symbolic shift/operator instances over 14 multiplicities and 12 exact finite radial-heat/Jensen models with free `rho,q,r,v,u1`. The independent checker specializes every finite model to separate rational data and reconstructs the operator without importing this builder.

## Live Handoff

Bound the actual regular field `ell` at every hypothetical positive-boundary multiple zero using the Xi canonical product, zero dynamics, or the Fourier kernel. In parallel, derive a local-uniform remainder with constants uniform along the selected degree sequence. The independent first-jet winding theorem remains the determinant-free alternative.

## Pi Provenance

No pi enters this local operator calculation. The rational `1/80` is the fifth-order finite-Jensen/heat normal-order coefficient. Pi enters only in separate Fourier normalizations of the Xi kernel.

## Proof Boundary

The fixed-m local expansion remains noncontradictory. Collision exclusion still requires a remainder uniform in the selected Jensen degree or a separate global first-jet winding theorem. This gate proves no Xi regular-field sign, degree-uniform collision exclusion, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
