# Newman Contact Signed-Transport Reduction

Date: 2026-07-28

Correction notice (2026-07-31): the historical endpoint atom and
`H_a!=0` division guard are superseded by the full complex endpoint atom
and zero-real-projection fibre in the source-normalization gate and Formal
Core Section 11.146.  The abstract signed-transport identities remain exact.

Status: exact endpoint-complete reduction and route guard;
not a proof of the Xi contact gap, `Lambda<=0`, RH,
PF-infinity, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contact_signed_transport_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contact_signed_transport_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contact_signed_transport_reduction.py
```

## Pi Provenance

The pi in a^2=x/(4*pi)+t/16, T_0, u_(N,x), and the cutoff recurrence is inherited from the completed-zeta normalization, Gaussian/theta Fourier normalization, and Riemann-Siegel saddle. The signed-mass partition, slope ordering, transport identity, and Gram audit introduce no new pi, fitted circle, or polygonal definition.

## Centered Components

```text
W_0=e+sum_n z_n=mathsf_X+i*mathsf_Y, W_A=g+s_*'*sum_n u_n*z_n, s_*'=c+i*b, mathsf_A=Re(W_A), alpha=c*u_N, and mathcal_C_N=mathsf_A-alpha*mathsf_X.
c_0=Re(e), d_0=Re(g)-alpha*Re(e); c_n=Re(z_n), d_n=c*(u_n-u_N)*Re(z_n)-b*u_n*Im(z_n); sum_(j=0)^N c_j=mathsf_X and sum_(j=0)^N d_j=mathcal_C_N, alpha=c*u_N..
e=kappa*H_a*(T_0+i), g=kappa*J_a*(T_0+i), with real kappa,H_a and generally complex J_a=H_(a,x)+mu_a*H_a.; c_0=kappa*T_0*H_a and d_0=kappa*[T_0*Re(J_a)-Im(J_a)-c*u_N*T_0*H_a].; If H_a!=0 and c_0!=0, the endpoint effective slope is h_0=d_0/c_0=Re(J_a)/H_a-Im(J_a)/(T_0*H_a)-c*u_N.; If H_a=0, the endpoint value vanishes but d_0=kappa*[T_0*Re(J_a)-Im(J_a)] remains in D_perp; no division by H_a is allowed..
```

## Contact Transport

```text
P={j:c_j>0}, N_-={j:c_j<0}, Z={j:c_j=0}; M_+=sum_P c_j, M_-=sum_N_(-c_j), h_+=sum_P d_j/M_+, h_-=sum_N_(-d_j)/M_-, D_perp=sum_Z d_j, and M=(M_++M_-)/2.
mathcal_C_N=M*(h_+-h_-)+(mathsf_X/2)*(h_++h_-)+D_perp.
On mathsf_X=0, M_+=M_-=M and mathcal_C_N=M*(h_+-h_-)+D_perp.
```

The contact equation is therefore a signed covariance, not
a geometry-only positive energy.

## Slope-Ordered Form

```text
Order the nonzero-real components so h_(1)<=...<=h_(m), where h_(j)=d_(j)/c_(j), and put P_k=sum_(j<=k)c_(j). mathcal_C_N=D_perp+h_(m)*mathsf_X-sum_(k=1)^(m-1)(h_(k+1)-h_k)*P_k. On mathsf_X=0, mathcal_C_N=D_perp-sum_(k=1)^(m-1)(h_(k+1)-h_k)*P_k.
```

The positive gaps are exact. Their cumulative real masses
are the open arithmetic input.

## Conditional Margin

```text
The exact band identity gives |mathcal_C_N|>=M*|h_+-h_-|-(|mathsf_X|/2)*|h_++h_-|-|D_perp|. Therefore on |mathsf_X|<=delta_L it is sufficient to prove M*|h_+-h_-|>A_L+epsilon_term+(delta_L/2)*|h_++h_-|+|D_perp|. Equivalently, the slope-ordered formula may be bounded directly before taking absolute values.
```

## Gram Rank Guard

For component amplitude variables a_j and real feature vectors g_j=(p_j,q_j/ell), mathsf_X^2+(mathcal_C_N/ell)^2=a^T*(p*p^T+q*q^T/ell^2)*a. The Gram matrix is positive semidefinite with rank at most two. On p.a=0 it reduces to (q.a)^2/ell^2, rank at most one; for m>=3 its contact kernel has dimension at least m-2. Squaring the scalar proves nonnegativity but supplies no positive lower bound.

## Endpoint-Shaped Null Guard

```text
1/16-i/2; e=(8+i)/100, g=e/8; z_1=omega, z_2=-(3/5)omega, z_3=-(2/5)omega, z_4=-e, omega=(-8+i)/sqrt(65); u_1,u_2,u_3,u_4=4,3,2,1; 1,3/5,2/5,sqrt(65)/100 are strictly decreasing. W_0=0; W_A=i*(7*sqrt(65)/80+13/320)!=0, so mathcal_C_N=Re(W_A)=0. The endpoint has the exact factored (T_0+i) shape with a real synthetic J_a/H_a, but the chosen H_a,J_a and carriers are synthetic and do not satisfy the full Xi recurrence or d_n phase chain.
```

## First Jet

```text
Let A=log(a), u_x=1/(8*pi*a^2), P_1=H_1-A*H_0, Z_0=r_0+H_0, Z_A=r_A-s_*'*P_1. Then Z_(0,x)=r_(0,x)-s_*'*H_1+D_0, P_(1,x)=-s_*'*(H_2-A*H_1)+D_1-A*D_0-u_x*H_0, and Z_(A,x)=r_(A,x)-s_*''*P_1-s_*'*P_(1,x). With eta_x=i*omega_eta*eta, mathsf_X_x=Pi_eta(Z_(0,x))-omega_eta*Im(eta*Z_0), mathsf_A_x=Pi_eta(Z_(A,x))-omega_eta*Im(eta*Z_A), and (mathcal_C_N)_x=mathsf_A_x-alpha_x*mathsf_X-alpha*mathsf_X_x, alpha_x=Re(s_*'')*u_N+c/(4*T_0). This retains H_0,H_1,H_2,D_0,D_1 together with r_0,r_A and all displayed endpoint derivatives.
```

## Cutoff Law

```text
For alpha_N=c*u_N and mathcal_C_N=mathsf_A_N-alpha_N*mathsf_X_N, mathcal_C_(N+1)-mathcal_C_N=Delta mathsf_A-alpha_(N+1)*Delta mathsf_X+(alpha_N-alpha_(N+1))*mathsf_X_N. On a contact of the N chart this becomes Delta mathcal_C=Delta mathsf_A-alpha_(N+1)*Delta mathsf_X. Delta Z_0=q_(N+1)+j_0 and Delta Z_A=-s_*'*log((N+1)/a)q_(N+1)+j_A; only their certified real projections are transferred.
```

The labels h_j=d_j/c_j are used only when c_j!=0. All c_j=0 contributions remain in the division-free D_perp. At a carrier cosine zero or H_a=0 the mass formula continues without a pole, while any h-order chart must end and restart.

## Legacy Symmetry Verdict

The legacy radial/transpose symmetry and graph-Laplacian positivity correspond only to a geometry-side norm or positive Gram construction. The contact restriction leaves one signed feature and a large nullspace. The missing input is arithmetic separation of the positive and negative effective-slope mass distributions, not another symmetric smoothing or skeletonization.

## Open Source Inequality

A proof may use either of two equivalent signed targets: (i) separate the contact-weighted means h_+ and h_- at the required scale while controlling D_perp, or (ii) order by effective slope and prove a one-sided cumulative-real-mass estimate for P_k strong enough that D_perp+h_m*mathsf_X-sum_k(h_(k+1)-h_k)P_k cannot enter [-(A_L+epsilon_term),A_L+epsilon_term]. Both targets include the endpoint atom and use the contact condition before absolute values.

## Route Decision

Retain the contact signed-transport formulation as a proof-facing diagnostic, not as a solved inequality. Partition a fixed-N physical chart only at real-component zeros and effective-slope order changes, keep the endpoint atom, and seek Xi-specific oscillatory control of the cumulative masses. Test every candidate first at q=1, W_0=0, H_a=0, the recurrent endpoint, and both sides of every cutoff recurrence. Do not infer coercivity from the positive Gram square.

The q<1 layer remains a separate multiplicity-compatible parabolic/Hermite or boundary-degree theorem. A uniform positive slope floor down to t=0 would impose endpoint simplicity and is not required by RH.

## Boundary

The terminal-centered component identity, endpoint atom, rate-free contact covariance, slope-ordered Abel transport, conditional margin, Gram-rank collapse, endpoint-shaped generic null guard, full five-current x jet, and adjacent centered-scalar jump are exact or reproducibly finite. The null guard is not an actual Xi counterexample. No Xi cumulative-mass estimate, Abel-scalar gap, successor inner degree theorem, q<1 closure, finite-height effectivity, contact exclusion, Lambda<=0, PF-infinity, RH proof, or Clay-prize conclusion is obtained.
