# Beta^-4 canonical fold: exact Airy-derivative reduction

Date: 2026-08-13
Status: exact Abel/Fourier reduction; not a proof of the ordinary-carrier bound

Put `X=lambda+y`, `epsilon=beta^-2`, and

```text
A(X)=Ai(-X),                    A_X=dA/dX,
D_M(y)=sum_(m in M) exp(-i*d_m*y).
```

For the corrected event normal form use

```text
M_4(z,y)=1+epsilon(a_1+i*r_1)
 +epsilon^2(a_2+i*r_2+i*a_1*r_1-r_1^2/2),            (AR1)

a_1=y/2-3z^2/4,
r_1=-2z^5/15+z^3y/3-zy^2/4,
a_2=13z^4/32-3z^2y/8,
r_2=17z^7/315-2z^5y/15+z^3y^2/12.
```

The real-line `z` integrals are symmetric Abel limits.  The Fourier
normalization of the Airy function gives, for `0<=k<=10`,

```text
integral_R^Abel z^k exp(i[z^3/3-Xz])dz
   =2*pi*i^k d^k A(X)/dX^k.                            (AR2)
```

This is the provenance of `pi`: it is exactly the `2*pi` in the Airy Fourier
transform, not a fitted constant or an arbitrary circle.  Since `A_XX=-X*A`,
write `A^(k)=P_k(X)A+Q_k(X)A_X`, where

```text
P_0=1, Q_0=0,
P_(k+1)=P_k'-X*Q_k,       Q_(k+1)=P_k+Q_k'.            (AR3)
```

Substitution through degree ten collapses the complete corrected two-variable
model to

```text
I_M^[4](beta,lambda,Y)
 =2*pi integral_0^Y exp(i*y^2/(4beta))*D_M(y)
    *[U(lambda,y,beta)A(lambda+y)
      +V(lambda,y,beta)A_X(lambda+y)]dy,               (AR4)

U=1+epsilon*U_1+epsilon^2*U_2,
V=  epsilon*V_1+epsilon^2*V_2,

U_1=-(13*lambda + 3*y)/60,
V_1=(8*lambda**2 - 4*lambda*y + 3*y**2)/60,
U_2=-(448*lambda**5 + 280*lambda**2*y**3 + 4565*lambda**2 - 105*lambda*y**4 + 30*lambda*y + 63*y**5 - 405*y**2)/50400,
V_2=-(-40*lambda**3 + 20*lambda**2*y - lambda*y**2 + 9*y**3 - 27)/1680.                                        (AR5)
```

Thus the beta^-4 fold model has only two scalar carrier channels, not an
unresolved two-dimensional oscillatory integral.  The leading formula in
Section 11.315 is recovered exactly by setting `epsilon=0`.

The remaining wall is now explicit.  One must prove an Airy branch expansion
with remainder in the ordinary corridors, match both phase branches and their
normalizations to the classical/Gamma carrier, retain the exact Airy form near
`X=0`, and prove overlap ownership without omission or duplication.  This gate
does not perform that asymptotic match and proves no bound for `Q_K-T`, no
complete `T_upper`, no height-uniform theorem, no `Lambda<=0`, no PF-infinity,
no RH, and no prize-level conclusion.
