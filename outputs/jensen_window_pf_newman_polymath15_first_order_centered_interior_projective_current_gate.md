# Interior Carrier Projective-Current Gate

Date: 2026-07-28

Correction notice (2026-07-31): the endpoint effective-slope text is
superseded by the full complex endpoint projection in the
source-normalization gate and Formal Core Section 11.146.  The proved
interior projective-current theorem and exceptional edge packaging remain
exact.

Status: proved interior-carrier monotonicity and exact edge isolation;
not a proof of the Xi cumulative-mass gap, `Lambda<=0`, RH,
PF-infinity, or a Clay-prize result.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_interior_projective_current_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_interior_projective_current_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_interior_projective_current_gate.py
```

## Pi Provenance

The pi in x=4*pi*exp(L), a^2=x/(4*pi)+t/16, and u_x=1/(8*pi*a^2) comes from the completed-zeta factor pi^(-s/2)*Gamma(s/2)*zeta(s) and the Riemann-Siegel saddle a^2=T_0/(2*pi). The use of pi>3 in one rational budget is an elementary bound on that same constant. No circle fitted to the carriers and no prefix polygon introduces pi.

## Carrier Current

```text
For z_n=X_n+iY_n and z_(n,x)=(varrho_n+i*nu_n)z_n, X_(n,x)=varrho_nX_n-nu_nY_n and Y_(n,x)=nu_nX_n+varrho_nY_n.
Put c_s=Re(s_*'), b=Im(s_*'), k_n=u_n-u_N=log(N/n), so k_(n,x)=0. Then c_n=X_n and d_n=c_s*k_n*X_n-b*u_n*Y_n.
J_n:=c_n*d_(n,x)-d_n*c_(n,x)=c_(s,x)*k_n*X_n^2-(b_x*u_n+b*u_x)*X_n*Y_n-b*u_n*nu_n*(X_n^2+Y_n^2). The radial amplitude current varrho_n cancels identically from J_n.
When X_n!=0, h_n=d_n/c_n=c_s*log(N/n)-b*u_n*tan(theta_n), and h_(n,x)=c_(s,x)*log(N/n)-(b_x*u_n+b*u_x)*tan(theta_n)-b*u_n*nu_n*sec(theta_n)^2=J_n/X_n^2.
```

## Uniform Interior Bound

The source gives nu_n<=8449/x^2-u_n/2, -b>=1/2, |b|<1, and u_x=1/(8*pi*a^2). Also s_*''=-t*alpha''(s)/8. Since alpha''(s)=1/s^3+2/(s-1)^3-1/(2s^2), |s|=|s-1|=sqrt(1+x^2)/2, t<=1/2, and x>1, |alpha''|<=26/x^2 and hence |c_(s,x)|,|b_x|<=13/(8*x^2).

For n<=N-1, u_n>=log(N/(N-1))>1/N>=1/a and 0<=k_n=log(N/n)<L. Since a^2<2*exp(L), u_n^2>exp(-L)/2. The monotone boundary checks 13L<pi^2*exp(L) and 33796a<x^2 follow for L>=50 from 650<9*2^50 and 67592<144*2^75. Therefore 8449/x^2<u_n/4, |c_(s,x)|k_n<u_n^2/64, |b_x|u_n<u_n^2/64, and u_x<u_n^2/(8*pi).

```text
Put B_n=b_x*u_n+b*u_x and K_n=b*u_n*nu_n. Then nu_n<-u_n/4, K_n>=u_n^2/8, and |B_n|<u_n^2/16 because 1/64+1/(8*pi)<1/64+1/24=11/192<1/16. Moreover J_n=[X_n,Y_n]M_n[X_n,Y_n]^T with M_n=[[c_(s,x)k_n-K_n,-B_n/2],[-B_n/2,-K_n]].
Gershgorin gives -1/8+1/64+1/32=-5/64 and -1/8+1/32=-3/32. Hence every eigenvalue of M_n is at most -5*u_n^2/64, and J_n<=-(5/64)*u_n^2*|z_n|^2<0 for 1<=n<=N-1.
```

## Projective Join

The interior vector v_n=(c_n,d_n) never vanishes: if X_n=0, then Y_n!=0 and d_n=-b*u_n*Y_n!=0. Its projective angle psi_n=arg(c_n+i*d_n) therefore satisfies psi_(n,x)=J_n/(c_n^2+d_n^2)<0. On X_n!=0, h_(n,x)=J_n/X_n^2<0.

At X_n=0 the ratio h_n is not used, but J_n=-K_n*Y_n^2<=-(1/8)*u_n^2*Y_n^2<0. Thus the division-free projective orientation crosses every interior real-projection zero with the same strict sign.

## Exceptional Edge

For n=N, k_N=0 and J_N=-(b_x*u_N+b*u_x)X_NY_N-b*u_N*nu_N*(X_N^2+Y_N^2). At the left cutoff u_N=0, J_N=-b*u_x*X_NY_N, which has either sign under the currently proved bounds. At u_N=0 and X_N=0, the isolated terminal component (c_N,d_N)=(X_N,0) is the zero vector, so it has no projective label. The choices b=-1/2,u_x=1,(X_N,Y_N)=(1,1) and (1,-1) give J_N=1/2 and -1/2. This is a local algebraic nonpromotion guard, not two asserted Xi cutoff states.

```text
Define the exceptional edge block C_edge=c_0+c_N and D_edge=d_0+d_N. Then mathsf_X=C_edge+sum_(n=1)^(N-1)c_n and mathcal_C_N=D_edge+sum_(n=1)^(N-1)d_n exactly. If C_edge=0, retain D_edge division-free. This packages the only carrier not covered by the interior theorem with the recurrent endpoint instead of assigning either one an unsupported sign.
```

Distinguish two source objects that were both denoted J_a. Inside one chart, J_a^(der)=H_(a,x)+mu_aH_a is generally complex, and h_0=Re(J_a^(der))/H_a-Im(J_a^(der))/(T_0H_a)-c_s*u_N only when H_a!=0. At an adjacent cutoff, J_a^(adj)=H_a(p+2)+H_a(p) and Q_N=f_(N+1)+kappa_NJ_a^(adj)=E_[1],N+1-E_[1],N. Recompute the terminal centering and transport the aggregate by mathcal_C_(N+1)-mathcal_C_N=Delta mathsf_A-alpha_(N+1)*Delta mathsf_X+(alpha_N-alpha_(N+1))*mathsf_X_N, with Delta A_a=Re[Q_(N,x)-lambda_aQ_N-e_(N+1)d_(N+1,x)]. Never identify J_a^(der) with J_a^(adj), and infer no separate endpoint or terminal projective sign.

## Order-Swap Guard

h_1(x)=-x and h_2(x)=-x+sin(x)/2 both decrease strictly, because h_1'=-1 and -3/2<=h_2'<=-1/2. Nevertheless h_2-h_1=sin(x)/2 changes sign at every successive multiple of pi. This hypothetical exact pair is not an Xi carrier pair. It proves only that individual monotonicity, as a theorem type, cannot by itself bound slope-order changes.

## Remaining Theorem

Strict clockwise flow of every interior atom does not sign the cumulative real masses P_k in the slope-ordered Abel identity, does not order different h_n, and does not bound their order changes. The remaining source theorem is still an Xi-specific one-sided cumulative-mass estimate, now with N-1 monotone interior projective atoms and one exact exceptional edge block.

Retain the negative-definite interior current as a new analytic lemma. In the next stage, formulate the aggregate transport on oriented projective cells for n<=N-1, but carry (C_edge,D_edge) by the exact endpoint recurrence. Seek summation-by-parts or oscillatory control of the actual Xi cumulative masses; do not infer it from monotonicity alone or split the terminal and endpoint signs.

The q<1 layer remains a separate multiplicity-compatible parabolic/Hermite or boundary-degree theorem. The interior current calculation is not promoted to t=0 endpoint simplicity.

## Boundary

This proves the exact division-free carrier current, radial-current cancellation, quantitative negative definiteness and strict projective orientation for every interior carrier on the stated q>=1 domain, including q=1. It also proves only an algebraic terminal sign obstruction, an exact edge-block decomposition, and a generic order-swap nonpromotion guard. It does not prove an Xi cumulative-mass estimate, an edge-block sign, the Abel-scalar gap, successor count, q<1 closure, finite-height effectivity, contact exclusion, Lambda<=0, PF-infinity, RH, or a Clay-prize conclusion.
