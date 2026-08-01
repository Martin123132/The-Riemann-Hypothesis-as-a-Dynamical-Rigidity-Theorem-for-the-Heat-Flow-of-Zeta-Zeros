# Physical q=1 Saddle Phase And Variance Reduction

Date: 2026-08-01

Status: exact source reduction with `0 signed physical bounds` and `0 Phi_B
bounds`. This is not a proof of RH; the endpoint-composed phase estimate is
open.

## Physical Parameters

```text
On q=2tL^2=1 with L>=50, t=1/(2L^2), x=4*pi*exp(L), a^2=exp(L)+1/(32L^2), and T_0=2*pi*a^2=x/2+pi/(16L^2).

sigma=Re(s_*)=1/2+1/(8L)+log(1+x^(-2))/(16L^2)-1/[4L^2(1+x^2)].
1/2+1/(8L)-1/(4L^2x^2)<sigma<1/2+1/(8L)-1/(16L^2x^2).

Omega=-Im(s_*)=x/2+atan(x)/(8L^2)-3x/[4L^2(1+x^2)]=T_0-epsilon.
epsilon=atan(1/x)/(8L^2)+3x/[4L^2(1+x^2)].
3/(8L^2x)<epsilon<7/(8L^2x).
```

Here `Omega=-Im(s_*)>0`; the `tau` of Section 11.161 is `-Omega`.

## Amplitude And Anchor

```text
delta_a=log(a)-Re(alpha)=delta_0+(1/2)log(1+z), delta_0=1/(1+x^2)-(1/4)log(1+x^(-2)), z=pi/(8L^2x).
0<delta_a<1/x^2+pi/(16L^2x).
B_a=1/2-delta_a/(4L^2).
1/2-1/(4L^2x^2)-pi/(64L^4x)<B_a<1/2; in particular B_a>49/100.
For u_n=log(a/n), A_n=A_a*exp[B_a*u_n+u_n^2/(8L^2)], A_a=|eta/S_a|*exp[t*log(a)^2/4-sigma*log(a)]>0.
partial_u log A(u)=B_a+u/(4L^2)>49/100 for u>=0. Thus n<m<a implies A_n>A_m.

The retained endpoint scale is S_a=kappa*T_0 with kappa=(-1)^N*B_0, B_0>0, and T_0>0. Hence S_a is real and sign(S_a)=(-1)^N.
eta/S_a=[(-1)^N*eta]/(B_0*T_0), so the physical common unit phase is omega_c=(-1)^N*eta=(-1)^N*phi*(1+d_1)/|1+d_1|.
Since tau=Im(s_*)=-Omega, w_n=omega_c*A_n*exp(i*Omega*log n). With omega_a=omega_c*exp(i*Omega*log a), this is exactly w_n=omega_a*A_n*exp(-i*Omega*u_n).
```

Thus the physical common phase is fixed, while the relative moments depend
only on the positive amplitude profile and `Omega`.

## Saddle Transfer

```text
Define H_k(xi)=sum_(n<=B)u_n^k*A_n*exp(-i*xi*u_n) and M_k=omega_a*H_k(Omega), 0<=k<=4.
For 0<=k<=4, |H_k(Omega)-H_k(T_0)|<=epsilon*sum_(n<=B)u_n^(k+1)A_n.
```

Termwise use |exp(-i*Omega*u)-exp(-i*T_0*u)|<=|Omega-T_0|u=epsilon*u. No cancellation or denominator is used. The actual phase is exponentially close to the exact cutoff saddle phase in every unnormalized moment. Relative moments still require a lower bound on H_0 and are not inferred.

## Two Carriers

```text
For w_j=A_j(c_j+i s_j), c_j^2+s_j^2=1, the pure Turan numerator T_2=[f*f''-(f')^2]/4 satisfies 4T_2=-A_1^2u_1^2-A_2^2u_2^2-A_1A_2[(u_1^2+u_2^2)c_1c_2+2u_1u_2s_1s_2].
If A_1>A_2>0 and u_1>u_2>0, then 4T_2<=(A_2-A_1)(A_1u_1^2-A_2u_2^2)<0.
```

The cross bracket is the bilinear form diag(u_1^2+u_2^2,2u_1u_2) on two unit vectors, so it is at least -(u_1^2+u_2^2). The Section 11.161 witness has u_1=2log(2)>u_2=log(2) but A_1=1/2<A_2=1, opposite to the physical q=1 ordering. It is therefore excluded by the exact physical amplitude law.

This signs the two-carrier pure Turan numerator. The full P_bulk^(0) also contains u_(N,x)f*g/2; the old real guard has g=0, but no arbitrary-phase full-current theorem is claimed.

## Phase Variance

```text
On H_0(Omega)!=0 put mu_k=H_k(Omega)/H_0(Omega) and v=mu_2-mu_1^2. For W(y)=omega_a*H_0(Omega+y)=R(y)exp(i theta(y)), W'/W=-i*mu_1, (log R)''=-Re(v), theta'=-Re(mu_1), and theta''=-Im(v).

Where f=R cos(theta)!=0, 4P_bulk^(0)/f^2=-Re(v)-[Re(mu_1)]^2 sec(theta)^2+[2u_(N,x)+Im(v)]tan(theta).

Define D_phase=Re(v)+[Re(mu_1)]^2 sec(theta)^2-[2u_(N,x)+Im(v)]tan(theta). Then P_bulk^(0)=-f^2 D_phase/4.
```

The branch-free ordinary-fibre form is

```text
For W=M_0=f+ig, M_1, M_2, R2=|W|^2, A=Re(M_1*conj(W)), and C=(M_2W-M_1^2)conj(W)^2, put D_hat=f^2Re(C)+A^2R2-[2u_(N,x)R2^2+Im(C)]fg. On W!=0, 4P_bulk^(0)R2^2=-D_hat, so P_bulk^(0)<=0 iff D_hat>=0.
```

At W=0 the relative moments are undefined, but the primary identity gives P_bulk^(0)=-(Im M_1)^2/4<=0.

## Stronger Nonpromotion Guard

```text
Take a=16k and n=k,2k,8k, so u=(4d,3d,d), d=log(2). Use the exact q=1 correction-free profile A(u)=A_a exp(B_a u+c u^2), 49/100<B_a<1/2, 0<c=1/(8L^2)<=1/20000.
The artificial common logarithmic half-turn tau=pi/d with omega=-exp(i*tau*log k) gives weights (-A_1,+A_2,+A_3).
Writing r=exp(B_a d), gamma=exp(c d^2), the normalized amplitudes are (r^4 gamma^16,r^3 gamma^9,r gamma). Uniformly on the q=1 box, f/A_1>3/50 and f''/(A_1 d^2)>9, while f'=g=0. Hence P_bulk^(0)>27 A_1^2 d^2/200>0.
At B_a=1/2,c=0 the amplitudes are (4,2sqrt(2),sqrt(2)) and P_bulk^(0)=[(-185+134sqrt(2))/2]log(2)^2>0.
```

This guard obeys the exact finite-q=1 amplitude profile, strict amplitude ordering, dyadic integer support, and one common logarithmic phase. Its tau and common phase are deliberately not the physical Omega and omega_c. It proves that amplitude ordering and profile shape alone cannot close the current.

### Pi Provenance

The guard's pi is only the half-turn condition tau*log(2)=pi. The physical pi in T_0=2*pi*a^2 comes independently from the completed-zeta/Riemann-Siegel and Poisson normalization.

## Exact Remaining Target

```text
Retain the exact upstream decomposition h^2*Phi_B=E+P_bulk^(0)+R_corr. On W!=0 this is h^2*Phi_B=E-D_hat/(4|W|^4)+R_corr. On W=0 it is h^2*Phi_B=E-(Im M_1)^2/4+R_corr.

On W!=0 prove the physical source inequality E+R_corr-D_hat/(4|W|^4)<=h^2/400, preferably h^2/800, using Omega=T_0-epsilon, the exact A(u), and omega_a.

On W=0 prove E+R_corr<=h^2/400+(Im M_1)^2/4. No relative moment or phase division is allowed there.
```

First seek a division-free estimate for D_hat and the composed E+R_corr at the exact saddle phase T_0, then transfer the five unnormalized moments to Omega with the epsilon bounds. The endpoint and transpose phase-sum family must remain composed.

## Boundary

This proves exact physical q=1 parameter and common-phase laws, strict correction-free amplitude ordering, five saddle-phase transfer bounds, a two-carrier pure-Turan ordering theorem, the ordinary-fibre phase-variance identity, its branch-free discriminant, and one ordered-profile nonpromotion guard. It does not prove D_hat>=0 for the Xi source, bound E or R_corr, upper-bound Phi_B, establish a signed joint Type-I/II or Vaughan estimate, exclude contact, transfer the retained model to Xi, prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level result.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction.py
```
