# Newman Arbitrary-Multiplicity Jensen Boundary-Layer Gate

Date: 2026-08-03

Status: exact arbitrary-multiplicity fixed-shift cofinal Jensen boundary layer; one degree-uniform Xi handoff open; not a proof of RH.

## Positive-Boundary Coordinate

F_(0,t)(s)=2*H_t(i*sqrt(s)); a finite real multiplicity-m zero c of H_t maps to the finite negative multiplicity-m zero rho=-c^2 of F_(0,t).

Positive-boundary attainment therefore removes spatial escape for the contradiction `Lambda>0`. The relevant Jensen sequence can be taken at the single shift `n=0`; only its degree must diverge.

## Exact Jensen Operator

For T_DF(z)=J_(D,n,lambda)(z/D), T_DF(z)=[(1+(z/D)*partial_x)^D F(x)]_(x=0). After translation to x=z, its logarithm is -z^2*partial_x^2/(2D)+sum_(r>=3)(-1)^(r+1)z^r*partial_x^r/(rD^(r-1)).

Under z=rho+eta/sqrt(D), the quadratic operator tends to exp(-(rho^2/2)*partial_eta^2); every r>=3 logarithmic term is O(D^(1-r/2)) and vanishes on compact eta sets.

## Heat Scaling

For lambda=lambda_*+tau/D, the coefficient PDE tends under the same scaling to exp(4*rho*tau*partial_eta^2); the drift term is O(D^-1/2).

## Universal Multiplicity Layer

If F_*(z)=(z-rho)^m U(z), U(rho)!=0, then D^(m/2)J_(D,n,lambda_*+tau/D)((rho+eta/sqrt(D))/D)/U(rho) tends locally uniformly to K_(m,a)(eta)=exp(a*partial_eta^2)eta^m, a=4*rho*tau-rho^2/2.

Equivalently,

```text
K_(m,a)(eta)=m!*sum_(q=0)^floor(m/2) a^q*eta^(m-2q)/(q!*(m-2q)!).
```

The machine audit checks every coefficient, Hermite specialization, imaginary specialization, and parameter orientation for multiplicities `2` through `16`: 79 coefficient checks and 15 checks on each side. It also computes 14 exact finite-Jensen heat models at shifts zero and three and recovers the same limit independently.

## Root Geometry

For a=-sigma^2/2<0, K_(m,a)(eta)=sigma^m*He_m(eta/sigma), so all m roots are real and simple.

For a=beta^2/2>0, K_(m,a)(eta)=(i*beta)^m*He_m(eta/(i*beta)); for m>=2 it has a nonreal conjugate pair (and only the root zero is real when m is odd).

The layer changes type only at a=0, hence tau=rho/8. Root continuity and a diagonal compactness choice give a sequence D_j->infinity of finite Jensen collisions with lambda_j=lambda_*+rho/(8D_j)+o(D_j^-1) and polynomial root w_j=rho/D_j+o(D_j^(-3/2)).

For `m=2`, `K_(2,a)=eta^2+2a=eta^2+8rho*tau-rho^2`, exactly recovering the existing double-zero gate. The present theorem removes the unproved assumption that a hypothetical Newman boundary zero is double.

## Newman Consequence

If Lambda>0, positive-boundary attainment supplies finite c and m>=2. At fixed shift n=0, rho=-c^2, so a sequence D_j->infinity has collisions t_j=Lambda-c^2/(8D_j)+o(D_j^-1). No unbounded-shift sequence is needed.

Thus the prize-relevant Jensen escape is fixed shift and unbounded degree. A separate large-shift compactness theorem is not needed for this contradiction route.

## Quantifier Guard

The certified Section 11.195 tail is on [-100,0], not the positive Newman interval. Even if extended separately at every fixed order, a finite double collision of degree D uses determinant order k=D; more generally bounded collision multiplicity gives k=D-m_D+2 tending to infinity. Fixed-order eventual thresholds therefore do not control this layer without a new order-uniform theorem.

The all-shift degree-361 sector theorem is bounded-degree evidence. The fixed-order eventual-Hankel theorem excludes attained finite collisions when its upper-shift hypotheses hold, but supplies no determinant threshold uniform in the growing order used here.

## Live Handoff

Prove an Xi-specific estimate uniform in degree throughout the cofinal rho/(8D) layer, or use the independent first-jet boundary-winding route. Bounded-degree certificates and fixed-order Hankel tails cannot be promoted to this target.

## Pi Provenance

The universal layer introduces no pi. The value `rho=-c^2` comes from the algebraic map `s=-z^2`; pi appears only in separate Xi kernel normalizations and estimates.

## Proof Boundary

This gate proves the arbitrary-multiplicity local limit and the conditional fixed-shift cofinal collision sequence under `Lambda>0`. It proves no degree-uniform exclusion, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH, or prize-level conclusion.
