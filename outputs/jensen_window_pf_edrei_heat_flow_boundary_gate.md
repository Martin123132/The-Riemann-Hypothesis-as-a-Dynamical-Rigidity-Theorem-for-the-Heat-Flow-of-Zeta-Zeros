# Jensen-Window PF Edrei Heat-Flow Boundary Gate

Date: 2026-07-23

Status: exact heat-flow lemma, exact backward-invariance countermodel,
and open Xi/Phi handoff. This is not a proof of PF-infinity, Jensen
hyperbolicity for zeta, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_edrei_heat_flow_boundary_gate`.

```text
work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json
python work/rh_compute/scripts/jensen_window_pf_edrei_heat_flow_boundary_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py
```

## Exact Heat Hierarchy

For the unnormalized squared-variable Phi transform,

```text
partial_lambda F_lambda=(4*z*partial_z^2+2*partial_z)F_lambda
```

and normalization does not change `R=F'/F`. Direct differentiation gives

```text
For R=partial_z log(F): partial_lambda R=4*z*R_zz+8*z*R*R_z+6*R_z+4*R^2
For R(z)=sum_(n>=0)(-1)^n*a_n*z^n: a_n'= -2*(n+1)*(2*n+3)*a_(n+1)+4*(n+1)*sum_(k=0)^n a_k*a_(n-k)
For A(t)=sum_(n>=0)a_n*t^n=R(-t): partial_lambda A=-4*t*A_tt+(8*t*A-6)*A_t+4*A^2
```

This is the heat hierarchy for the exact endpoint moments
`a_n=p_(n+1)`, not for the original linear signed-Hankel sequence.

## Rank-One Boundary Orientation

At a repeated type-I factor,

```text
At H(z)=(1+beta*z)^m, a_n=m*beta^(n+1) and every shifted 2x2 Stieltjes minor Delta_s=a_s*a_(s+2)-a_(s+1)^2 vanishes
At that boundary, partial_lambda Delta_s=8*m^2*(m-1)*beta^(2*s+5)
```

Thus a multiplicity `m>=2` points strictly into every shifted `2x2`
Stieltjes wall under forward heat and strictly out under backward heat.

## Exact Double-Zero Countermodel

For `L=4z partial_z^2+2 partial_z`, the heat series terminates:

```text
exp(lambda*L)(1+beta*z)^2=beta^2*z^2+(2*beta+12*lambda*beta^2)*z+1+4*lambda*beta+12*lambda^2*beta^2
Disc_z=32*lambda*beta^3*(1+3*lambda*beta)
Delta_s(lambda)=32*lambda*beta^(2*s+5)*(1+3*lambda*beta)/(1+4*lambda*beta+12*lambda^2*beta^2)^(s+3)
```

The denominator in the shifted-minor formula is positive for real
`lambda`. Hence, for `-1/(3*beta)<lambda<0`, every shifted minor is
negative and the two zeros are nonreal. At `lambda=0` they collide at
`z=-1/beta`; for `lambda>0` they are distinct and negative.

Exact beta=1 witnesses:

```text
lambda=1/10:  Disc=104/25, Delta_0=8125/6859
lambda=-1/10: Disc=-56/25, Delta_0=-4375/729
```

## Consequence

For -1/(3*beta)<lambda<0 every Delta_s(lambda)<0 and the quadratic zeros are nonreal; generic backward Stieltjes-cone invariance is false.

Therefore the known real-zero regime at a positive de Bruijn time cannot
be propagated to lambda zero by a generic backward-invariance assertion
for the Stieltjes cone. The surviving target is:

```text
A backward argument from a de Bruijn real-zero time must use Xi/Phi-specific all-order rigidity or a uniform no-escape theorem, not generic invariance of the Stieltjes moment cone
```

The countermodel is a finite LP+ polynomial and not the Xi function. It
rejects only the generic shortcut; it leaves genuinely Xi/Phi-specific
rigidity open.

## Proof Boundary

This artifact proves the radial log-derivative and Edrei-moment heat hierarchies and gives an exact LP+ polynomial countermodel to generic backward Stieltjes-cone invariance. It does not prove an Xi/Phi-specific invariant, all-order Stieltjes positivity at lambda zero, PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0.
