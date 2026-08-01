# Jensen-Window PF Phi Pick-Kernel Target

Date: 2026-07-23

Status: exact endpoint-equivalent Pick-kernel target, exact first-wall
identity, and positive-mixture guard. This is not a proof of PF-infinity,
Jensen hyperbolicity for zeta, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_phi_pick_kernel_target`.

```text
work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json
python work/rh_compute/scripts/jensen_window_pf_phi_pick_kernel_target.py
python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py
```

## Exact Pick Reformulation

Set

```text
F(z)=integral_R Phi(u)*cosh(u*sqrt(z))*du, F(x)>0 for x>=0
C_u(z)=cosh(u*sqrt(z)), D_u(z)=partial_z C_u(z)=u*sinh(u*sqrt(z))/(2*sqrt(z))
R(z)=F'(z)/F(z)=integral Phi(u)D_u(z)du / integral Phi(u)C_u(z)du
```

For `Im(z)>0`, define

```text
P_Phi(z):=-Im(F'(z)*conj(F(z)))=|F(z)|^2*(-Im R(z))
K_z(u,v):=-(1/2)*Im(D_u(z)*conj(C_v(z))+D_v(z)*conj(C_u(z)))
P_Phi(z)=double_integral Phi(u)*Phi(v)*K_z(u,v)du dv
```

Christian Berg's Theorem 3.2 gives the Stieltjes/Pick characterization

```text
R is Stieltjes <=> R(x)>=0 for x>0, R is holomorphic on C\(-infinity,0], and Im R(z)<=0 for Im(z)>0
```

while Sokal's logarithmic-derivative criterion identifies `R` being
Stieltjes with `F` being LP+. Since the nonexponential Phi transform is
positive on the nonnegative real axis, these facts give the exact target

```text
For the nonexponential Phi transform: F in LP+ <=> P_Phi(z)>0 for every Im(z)>0
```

Strict positivity excludes an upper-half-plane zero because such a zero
would make `P_Phi=0`; reflection excludes lower-half-plane zeros, and
`F(x)>0` excludes nonnegative real zeros.

Under `z=w^2`, `w=a+i*b`, this becomes

```text
For z=w^2 and w=a+i*b with a,b>0: 2*|w|^2*P_Phi(w^2)=-Im(conj(w)*M'(w)*conj(M(w))), M(w)=integral Phi(u)cosh(u*w)du
```

This is a real two-parameter inequality in the known Phi kernel. It is
endpoint-equivalent, not a completed estimate.

## First Real-Axis Wall

Push the normalized Phi scale measure to `y=u^2` and tilt it by
`C_x(y)=cosh(sqrt(x*y))`. Then

```text
For y=u^2, C_x(y)=cosh(sqrt(x*y)), r_x(y)=partial_x log C_x(y), and dPi_x=C_x*dmu/F(x): R'(x)=E_Pi[r_x']+Var_Pi(r_x)
At x=0: R'(0)=(E[y^2]-3*E[y]^2)/12, so the first real-axis Stieltjes wall is E[y^2]<=3*E[y]^2
```

Thus even the first boundary wall is a concentration inequality: mixing
variance must not exceed the average fixed-scale negative curvature.
Passing this wall for every x is still only necessary for the global
upper-half-plane Pick condition.

## Exact Positive-Mixture Guard

For the positive two-scale law

```text
(9/10)*delta_(1/4)+(1/10)*delta_(5/2) in y=u^2
E[y]=19/40
E[y^2]=109/160
E[y^2]-3E[y]^2=7/1600
R'(0)=7/19200>0
```

so the Stieltjes derivative sign and the Pick numerator fail near
`z=i*epsilon`. This proves that the polarized cross-scale kernel cannot
be pointwise nonnegative in general. The witness is not the Xi density.

## Surviving Structural Routes

A direct global route is

```text
Prove P_Phi(z)>0 on the whole upper half-plane, or construct the positive resolvent representation, using verified Xi/Phi structure rather than zero reality
```

An operator route is

```text
A direct identity R(z)=gamma+<v,(I+z*T)^(-1)v> with T positive self-adjoint would give the required Stieltjes representation by the spectral theorem
```

The spectral theorem would then supply the positive Stieltjes measure;
entireness and Sokal's criterion would recover the discrete Edrei form.
Constructing such an operator from Phi is open and cannot use the zero
reality it is meant to prove.

## Sources

- Christian Berg, `Stieltjes-Pick-Bernstein-Schoenberg and their connection to complete monotonicity`, Theorem 3.2: https://web.math.ku.dk/~berg/manus/castellon.pdf
- Kalmykov and Karp, `When does a hypergeometric function belong to LP+?`, Proposition 6: https://doi.org/10.1016/j.jmaa.2022.126432
- `outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md`

## Proof Boundary

This artifact proves an exact Pick-kernel reformulation of the common LP+/coefficient-PF/Jensen/Stieltjes endpoint, a tilted first-wall identity, and an exact positive-mixture countermodel. It does not prove the global Phi-kernel sign, construct a positive resolvent, prove PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0.
