# Newman Phase-Cylinder Jacobian Transport Guard

Date: 2026-07-27

Status: exact Jacobian and winding-transport reduction with
route guards; not a proof of `Lambda<=0`, RH, PF-infinity,
or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_phase_cylinder_jacobian_transport_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_phase_cylinder_jacobian_transport_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_phase_cylinder_jacobian_transport_guard.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16 comes from the completed zeta normalization and Riemann-Siegel saddle; 4*pi(2N+1) is the cutoff-cell difference and 2*pi is the period of exp(i theta). The phase cylinder and its Jacobians do not define pi.

## Phase Cylinder

```text
Retain P_N(z,xi,mu)=sum_(k=0)^K z^k*sum_(m<=M_k,m odd)a_(k,m)exp(i xi log m)*(1+mu d_(2^k m))/(1+mu d_1), with z=exp(i theta). The synchronized base is (xi,mu)=(0,0); the actual finite prefix is (theta,xi,mu)=(omega_*log 2,omega_*,1).
P_theta=i*sum_alpha k_alpha w_alpha exp(i phi_alpha), P_xi=i*sum_alpha l_alpha w_alpha exp(i phi_alpha)
```

## Pair-Kernel Jacobian

```text
J_(theta,xi)=Im(conj(P_theta)*P_xi). It is the oriented real Jacobian det partial_(theta,xi)(Re P_N,Im P_N).
sum_(alpha<beta) w_alpha*w_beta*(k_alpha*l_beta-k_beta*l_alpha)*sin(phi_beta-phi_alpha), with phi_alpha=k_alpha*theta+l_alpha*xi, l_alpha=log(m_alpha). For complex correction coefficients, replace w_alpha*w_beta*sin(phi_beta-phi_alpha) by Im(conj(c_alpha)c_beta exp(i(phi_beta-phi_alpha))).
```

This Jacobian is the local degree density for unit-circle
zeros crossing as the odd Mellin phase changes.

## Physical Linked Direction

```text
On theta=h*xi, dP/dxi=h*P_theta+P_xi=i*mathcal D_hP=i*H_1.
Where H_0=P is nonzero, partial_xi arg P_link=Re(conj(H_0)H_1)/|H_0|^2.
At the actual point, H_(0,x)=-s_*'H_1+D_0. Thus the linked phase current iH_1 is only the phase component; the sigma_* amplitude drift and full d_n current D_0 must be restored before an x-transversality claim.
```

The phase-only linked direction is exact, but physical x
motion also changes amplitudes and correction factors.

## Winding Transport

```text
Assume P_N is nonzero on the two boundary circles xi=xi_0,xi_1 and has only isolated regular zeros in S^1_theta times (xi_0,xi_1). With J_(theta,xi)=Im(conj(P_theta)P_xi), the degree theorem gives wind_theta P(.,xi_1)-wind_theta P(.,xi_0)=-sum_(P=0)sign J_(theta,xi). For P=exp(i theta)-r(xi) with r increasing through one, J>0 at the crossing while wind_theta drops from one to zero.
```

This is the surviving role of the Jacobian. It counts signed
boundary crossings; it is not globally positive.

## Sign Guards

### Genuine dyadic/odd pair

```json
{
  "interpretation": "Even genuine dyadic/odd support with positive amplitudes has an oscillating pair current. The full Jacobian may still have structure, but no termwise-positive proof is available.",
  "linked_kernel": "J_23=w_2*w_3*log(3)*sin(xi*log(3/2))",
  "negative_point": "xi*log(3/2)=3*pi/2 modulo 2*pi",
  "positive_amplitudes": "w_2>0, w_3>0",
  "positive_point": "xi*log(3/2)=pi/2 modulo 2*pi",
  "support": "n=2=(k,m)=(1,1), n=3=(k,m)=(0,3)"
}
```

### Zero-free but nonmonotone linked block

```json
{
  "first_moment": "H_1=log(2)*r*exp(i*phi)",
  "interpretation": "A zero-free complete dyadic block need not have monotone argument along the linked physical phase.",
  "linked_argument_numerator": "Re(conj(P)*H_1)=log(2)*(r^2+r*cos(phi))",
  "negative_value": "at phi=pi it equals log(2)*(1/2-1/sqrt(2))<0",
  "nonvanishing": "|P(phi)|>=1-r>0",
  "prefix": "P(phi)=1+r*exp(i*phi), r=2^(-1/2)"
}
```

## Correction Cylinder

```text
P_mu=sum_(k,m) z^k a_(k,m)exp(i xi log m)*(d_(2^k m)-d_1)/(1+mu d_1)^2. Define J_(theta,mu)=Im(conj(P_theta)P_mu). If the boundary circles are nonzero and all interior zeros are regular, the same orientation gives wind_theta P(.,mu_1)-wind_theta P(.,mu_0)=-sum_(P=0)sign J_(theta,mu).
```

## Endpoint Cylinder

```text
At fixed physical x, after the finite prefix is restored, define F_N(z,tau)=P_N(z,omega_*,1)+tau*r_0 for 0<=tau<=1. Then F_theta=P_theta, F_tau=r_0, and J_(theta,tau)=Im(conj(P_theta)r_0). The same cylinder degree formula transports the finite-prefix winding to the endpoint-augmented value Z_0. This endpoint cylinder is indispensable.
```

The need for this final transport is exact:

```json
{
  "augmented_zero": "P(1)+e=0",
  "endpoint": "e=-5/3",
  "interpretation": "Unit-disk zero-freeness of the finite prefix does not survive an arbitrary endpoint addition. The actual Riemann-Siegel endpoint needs its own transport current.",
  "prefix": "P(z)=1+(2/3)z",
  "prefix_zero": "z=-3/2, outside the closed unit disk"
}
```

## Three-Cylinder Programme

```text
A sufficient strong route starts from the proved P_N(z,0,0), transports odd phase xi:0->omega_*, then correction mu:0->1, then endpoint tau:0->1. Each stage has an exact local Jacobian and signed winding-change formula. The final physical theorem may use a weaker linked-point argument, but no stage may be skipped merely because the preceding circle was zero-free.
The final physical derivatives remain Z_(0,x)=r_(0,x)-s_*'H_1+D_0 and Z_(A,x)=r_(A,x)-s_*''P_1-s_*'P_(1,x), where P_(1,x)=-s_*'(H_2-log(a)H_1)+D_1-log(a)D_0-H_0/(8*pi*a^2). After eta projection, these enter Gamma_ell=mathsf_X+i mathsf_A/ell and its exact argument flux. The adjacent recurrence is retained at cutoff boundaries.
```

## Route Decision

```text
Reject a global or termwise sign theorem for J_(theta,xi): its genuine n=2,3 pair kernel oscillates. Reject monotone linked argument: even the zero-free 1+2^(-1/2)exp(i phi) block has negative angular current on part of the phase cycle. Retain the Jacobians as local degree currents. The live theorem is a signed crossing budget, no-zero cylinder theorem, or Xi-specific cancellation theorem across the odd-phase, correction, and endpoint cylinders.
```

## Live Theorem

```text
On q=2tL^2>=1, prove enough boundary nonvanishing and signed-degree control for the three cylinders to transport the synchronized base to the endpoint-augmented physical first jet. A sufficient result is zero interior degree together with nonvanishing final boundary; a more local result may control only the linked point and complete successor path. It must still imply |mathsf_X|<=delta_L => |mathcal C_N|>A_L+epsilon_term and 0<=kappa_j<1.
```

The separate small-`q` obligation is

```text
The q=2tL^2<1 layer still needs the separate multiplicity-compatible parabolic/Hermite first-jet chart. The cylinder identities remain exact, but no endpoint simplicity or global Jacobian sign is imposed there.
```

## Boundary

The phase partials, pair-kernel Jacobian, linked derivative and argument current, cylinder winding-transport theorem, correction and endpoint cylinders, and three route guards are exact. No sign or nonzero lower bound for the full Xi Jacobian, no three-cylinder crossing budget, no joined Xi lower bound, no strict successor flux upper bound, no q<1 closure, no contact exclusion, no Lambda<=0 or PF-infinity theorem, no RH proof, and no Clay-prize conclusion is asserted.
