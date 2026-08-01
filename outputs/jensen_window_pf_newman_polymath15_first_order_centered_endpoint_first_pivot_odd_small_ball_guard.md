# Newman Endpoint First-Pivot Odd Small-Ball Guard

Date: 2026-07-27

Correction notice (2026-07-31): the sentence declaring `H_a` real is
superseded by the complex endpoint source-normalization gate and Formal
Core Section 11.146.  The complex small-ball reduction and phase formula
remain exact.

Status: exact first-pivot reduction and linked-phase route
countermodel; not a physical Xi pivot counterexample and not a
proof of `Lambda<=0`, PF-infinity, RH, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_first_pivot_odd_small_ball_guard.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_endpoint_first_pivot_odd_small_ball_guard.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_first_pivot_odd_small_ball_guard.py
```

## Pi Provenance

The pi in beta and T_0 comes from the completed-zeta normalization and Riemann-Siegel saddle. The 2*pi used for phase reduction is exactly the period of exp(i*theta). The dyadic split, common-denominator cancellation, Schur pivot, and triangle geometry introduce no new pi.

## Terminal Singleton

```text
Since 2^K<=N<2^(K+1), M_K=floor(N/2^K)=1. Hence the top dyadic layer contains only m=1 and B_(K,0)=C_K/(1+d_1), where C_K=exp[t*(K*h)^2/4-sigma_*K*h]*(1+d_(2^K)).
```

## Common Numerators

```text
Define O_N=sum_(m<=N,m odd)exp[t*log(m)^2/4-s_*log(m)]*(1+d_m). Then B_(0,0)=O_N/(1+d_1).
Define R_N=(-1)^N*beta*(T_0+i)*H_a/M_t(s). The branch-free phase-anchor formula gives r_0=R_N/(1+d_1), with beta>0 and H_a real.
Therefore b_0=r_0+B_(0,0)=(O_N+R_N)/(1+d_1) and b_K=B_(K,0)=C_K/(1+d_1). For nonzero D, |(O+R)/D|^2-|C/D|^2=(|O+R|^2-|C|^2)/|D|^2.
```

## Exact First Pivot

```text
Delta_K={|O_N+R_N|^2-rho_K^2}/|1+d_1|^2, where rho_K=exp[t*(K*h)^2/4-sigma_*K*h]*|1+d_(2^K)|. Thus the explicit normalization denominator affects neither the sign nor the zero set of the first pivot. The m=1 occurrence of d_1 remains inside O_N.
Delta_K>0 iff |O_N+R_N|>rho_K, equivalently O_N is outside the closed disk Dbar(-R_N,rho_K). The first Schur obligation is an endpoint-shifted odd Dirichlet-polynomial small-ball exclusion.
|O+R|^2=|O|^2+|R|^2+2Re(O*conj(R)). Hence Delta_K>0 iff |O_N|^2+|R_N|^2+2Re(O_N*conj(R_N))>rho_K^2.
```

## Endpoint Phase

```text
When H_a!=0, arg(R_N)=N*pi+arg(H_a)+arg(T_0+i)-arg(M_t(s)) modulo 2*pi; when H_a=0, R_N=0. This is branch-free phase bookkeeping, not a lower bound or an alignment theorem for O_N.
Any one of ||O_N|-|R_N||>rho_K or |Re(u*(O_N+R_N))|>rho_K for a proved unit complex u is sufficient. Equivalently one may prove the displayed correlation inequality. No certified current source supplies one of these uniformly.
```

## Linked-Phase Countermodel

```text
Set t=0, sigma_*=1/2, d_n=0, R_N=0, and N=5. Then K=2, rho_2=1/2, and O_5(omega)=1+3^(-1/2)exp(i*omega*log(3))+5^(-1/2)exp(i*omega*log(5)).
The continuous flow omega -> (omega*log(3),omega*log(5)) mod 2*pi is dense in the two-torus. By Kronecker's theorem it is enough that no nonzero integers u,v satisfy u*log(3)+v*log(5)=0; exponentiation would give 3^u*5^v=1, impossible by unique prime factorization.
The two rotating vectors have lengths a=1/sqrt(3) and b=1/sqrt(5), and their sums fill the closed annulus |a-b|<=|w|<=a+b. The exact inequalities |a-b|<1<a+b hold. a+b>1 because 30>7sqrt(15), as 900>735.
The annulus contains -1, so torus density gives inf_omega|O_5(omega)|=0. Therefore some linked logarithmic phases satisfy |O_5(omega)|<1/2=rho_2 and Delta_2<0. The diagnostic omega=62643/100 gives modulus 8.17177552525160484e-04; the exact conclusion uses density, not this decimal.
```

## Small-Endpoint Extension

```text
More generally, for fixed endpoint e, O_(5,e)=1+e+a*exp(i*omega*log(3))+b*exp(i*omega*log(5)) has infimum zero whenever |a-b|<|1+e|<a+b. In particular this holds for every |e|<delta_0=a+b-1=0.02456386468958370379098351424821270273573. Thus endpoint smallness alone does not restore the pivot.
```

The N=5 model preserves the actual odd logarithmic phase linkage and correction-free Dirichlet amplitudes, but omega is allowed to range independently of the physical cutoff and the recurrent Xi endpoint is replaced by zero or fixed e. It is an exact route countermodel, not a counterexample to an actual Xi pivot, Newman, or RH.

## Adjacent-Chart Boundary

The exact cutoff update Q_N=f_(N+1)+kappa_N*J_a controls the physical value jump. Existing adjacent certificates bound selected real projections, not the complex displacement of O_N relative to -R_N, so they do not imply disk avoidance.

## Route Decision

Retain the first pivot as a precise falsification coordinate, but stop treating decreasing amplitudes, logarithmic phase linkage, or endpoint smallness as candidate proofs. Full unit-disk stability now requires an Xi-specific correlation or small-ball theorem coupling O_N to R_N. Unless the source formulas reveal such a coupling, prioritize actual linked-point nonvanishing or the signed three-cylinder degree route, which asks less than control of every free z on |z|=1.

## Live Theorems

```text
Either prove uniformly on every q>=1 fixed-N Xi chart that O_N avoids Dbar(-R_N,rho_K), using the actual relation among omega_*, N, d_n, H_a, and M_t(s), or produce a certified physical chart where the first pivot fails and retire full disk stability.
For the weaker primary route, prove nonvanishing only at z_*=exp(i*omega_*log(2)), or prove the endpoint-complete signed crossing budget through the three exact transport cylinders while retaining H_2,D_0,D_1 and the recurrent endpoint first jet.
The q=2tL^2<1 multiplicity-compatible parabolic/Hermite chart remains separate. This first-pivot analysis supplies no small-q closure, finite connector, or chart join.
```

## Boundary

The terminal singleton, odd and endpoint numerators, common-denominator cancellation, pivot disk equivalence, correlation form, endpoint phase bookkeeping, linked-flow density guard, triangle-annulus pivot failure, and endpoint-neighborhood extension are exact. No uniform actual Xi small-ball exclusion, physical pivot failure, later Schur pivot, full endpoint disk stability, linked-point lower bound, signed crossing theorem, strict successor flux bound, q<1 closure, contact exclusion, Lambda<=0, PF-infinity, RH proof, or Clay-prize conclusion is asserted.
