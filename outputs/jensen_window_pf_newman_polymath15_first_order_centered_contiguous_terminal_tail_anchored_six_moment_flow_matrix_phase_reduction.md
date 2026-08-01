# Six-Moment Flow Matrix and Reciprocal Phase Reduction

Date: 2026-08-01

Status: exact matrix/phase reduction with `0 signed flow bounds`, `0
stationary-phase remainder bounds`, and `0 Phi_B bounds`; this is not a proof
of RH.

## Domain

Use the q=1, L>=50, fixed-N chart, B=N-M-1, and the exact auxiliary frequency segment T_0-epsilon<=xi<=T_0 from Section 11.164. All physical tail jets and coefficient rows are held fixed along xi.

## Six-Moment Shift

For 0<=j<=4, G_j'(xi)=-i*log(a)G_j(xi)+iG_(j+1)(xi). Thus the derivative of a degree-four observation uses exactly G_0,...,G_5 and no G_6.

For r=(r_0,...,r_4), pad(r)=(r_0,...,r_4,0), S(r)=(0,r_0,...,r_4), and delta_a(r)=i[S(r)-log(a)pad(r)]. Then d_xi Re(r dot G_[0:4])=Re(delta_a(r) dot G_[0:5]).

For y=(Re G_0,...,Re G_5,Im G_0,...,Im G_5)^T, let E select G_0,...,G_4, P select G_1,...,G_5, and J_5=[[0,-I_5],[I_5,0]]. Then x=Ey and x'=D_ay with D_a=J_5(P-log(a)E).

## Eight-Observation Flow

Retain the Section 11.159 complex coefficient rows v_0,n_0,a_0,q_0 and put dot(r)=delta_a(r). For y in R^12, V=tilde(v)^Ty, N=tilde(n)^Ty, A=tilde(a)^Ty, Q=tilde(q)^Ty and V_xi=dot(v)^Ty, N_xi=dot(n)^Ty, A_xi=dot(a)^Ty, Q_xi=dot(q)^Ty.

Psi'=V_xi*N+V*N_xi-A_xi*Q-A*Q_xi+X_T*N_xi+A_(T,x)*V_xi-A_T*Q_xi-X_(T,x)*A_xi. Equivalently, Psi'=V_xi[N+A_(T,x)]+[X_T+V]N_xi-A_xi[Q+X_(T,x)]-[A_T+A]Q_xi.

Let U_0=[tilde(v),tilde(n),tilde(a),tilde(q)], U_1=[dot(v),dot(n),dot(a),dot(q)], J=(1/2)diag_offdiag(+1,-1), and t_T=(A_(T,x),X_T,-X_(T,x),-A_T)^T. Then Psi'=y^T M_xi y+ell_xi^T y, where M_xi=U_0JU_1^T+U_1JU_0^T and ell_xi=U_1t_T.

With W=[U_0,U_1] and K=[[0,J],[J,0]], M_xi=WKW^T. Hence rank(M_xi)<=8, ell_xi lies in range(U_1), and Psi' depends on at most eight real observations.

Every y satisfying U_0^Ty=U_1^Ty=0 has Psi'(xi)=0. The common flow-invisible kernel has dimension at least four.

On X_T+V=A_T+A=0, Psi'=V_xi[N+A_(T,x)]-A_xi[Q+X_(T,x)]. This is division-free but has no algebraic sign.

## Rank and Inertia

An exact algebraic specialization log(N)=2, log(a)=3, rho_1=1+i, rho_2=2-i, rho_(1,x)=1-2i, rho_(2,x)=-1+i, s_*'=1-i/2, s_*''=2+i, chi_N=1+i, b=-1/2, b_x=1, and u_N=u_(N,x)=1 has rank(W)=8 and leading 8-by-8 minor 51095/16.

Since K=[[0,J],[J,0]] has inertia (4,4), the witness flow matrix has inertia (4,4,4). Therefore the rank-eight bound is generically sharp and the flow is not coefficient-blind semidefinite. The witness is algebraic, not a physical Xi state.

## Common Phase

Write G_j(xi)=c_xi*Ghat_j(xi), where c_xi=exp[i(Omega-xi)log(a)] and Ghat_j contains the full Mangoldt or balanced divisor sum but no c_xi. Then Ghat_j'(xi)=iGhat_(j+1)(xi), while c_xi'=-i log(a)c_xi.

The endpoint family carries c_xi, the Hermitian family carries c_xi*conj(c_xi)=1, and the transpose family carries c_xi^2. Hence log(a) cancels exactly from the Hermitian flow and remains only through c_xi' and (c_xi^2)' in the endpoint and transpose families. No common phase enters a divisor sum.

For J_5=[[0,-I],[I,0]], decompose a real symmetric current matrix M as M_H=(M-J_5MJ_5)/2 and M_T=(M+J_5MJ_5)/2. Then M_HJ_5=J_5M_H, M_TJ_5=-J_5M_T, and the common-rotation quadratic flow is log(a)(J_5M-MJ_5)=-2log(a)M_TJ_5. The Hermitian part contributes zero.

In logarithmic coordinates ell_n=log(n), the three exact multipliers are i(ell_n-log(a)), i(ell_n-ell_m), and i(ell_n+ell_m-2log(a)) for endpoint, Hermitian, and transpose respectively.

## Reciprocal Phases

Use e(t)=exp(2*pi*i*t). The T_0=2*pi*a^2 factor is inherited from the completed-zeta/Riemann-Siegel saddle. Dividing its radian phase by the standard Poisson 2*pi gives a^2. These are the two appearances of the same normalization; no circle or fitted pi is introduced.

At xi=T_0 and with Poisson character e(-kn), discard only n-independent unit phases when locating saddles. The endpoint cycle phase is a^2log(n)-kn. For a two-variable mode, use a^2log(n)-a^2log(m)-kn+rm in the Hermitian family and a^2log(n)+a^2log(m)-kn-rm in the transpose family, with k,r>0.

The endpoint mode has n_*=a^2/k, phase a^2[log(a^2/k)-1], curvature -k^2/a^2, and flow multiplier log(a/k). Its retained-support window 1<=n_*<=B is exactly a^2/B<=k<=a^2.

The Hermitian mode has (n_*,m_*)=(a^2/k,a^2/r), stationary phase a^2log(r/k), Hessian diag(-k^2/a^2,+r^2/a^2), and flow multiplier log(r/k). Thus the multiplier is exactly the stationary phase divided by a^2.

On the reciprocal diagonal k=r, both the Hermitian stationary phase and its differentiated-flow multiplier vanish exactly. The zero-phase reciprocal diagonal channel is therefore absent from Psi'_H, although near-diagonal modes still require an estimate.

The transpose mode has (n_*,m_*)=(a^2/k,a^2/r), phase a^2[log(a^4/(kr))-2], Hessian diag(-k^2/a^2,-r^2/a^2), and flow multiplier log(a^2/(kr)).

Because B<a, every interior retained saddle has k,r>a. Hence the endpoint multiplier log(a/k) and transpose multiplier log(a^2/(kr)) are negative there. They multiply complex coefficients and unit phases, so this is not a current-sign theorem. The Hermitian multiplier is antisymmetric under k<->r.

Since Omega=T_0-epsilon, c_(T_0)=exp[-i epsilon log(a)]. The endpoint Poisson family retains this unit phase, the Hermitian family has none, and the transpose family retains exp[-2i epsilon log(a)].

For each primal variable in [1,B], the only continuous stationary modes satisfy a^2/B<=mode<=a^2. Modes outside this window are nonstationary, but a discrete Poisson formula, hard-cutoff boundary terms, and uniform integration-by-parts constants remain to be proved.

## Handoff

Insert the balanced Mangoldt square-and-wing formula for each Ghat_0,...,Ghat_5 into the eight observations (V,N,A,Q,V_xi,N_xi,A_xi,Q_xi). Keep c_xi outside the sums. Apply one- and two-variable Poisson/B-process transforms before taking absolute values, preserving endpoint, Hermitian, and transpose terms in one signed expression.

The canonical zero phase-difference channel is the Hermitian reciprocal diagonal, and the differentiated kernel kills it exactly. The next estimate should exploit the resulting log(r/k) factor by dual antisymmetry or summation by parts. The endpoint and transpose saddles remain and must be composed with the certified terminal recurrence rather than bounded separately.

Prove Psi(T_0)-epsilon*integral_0^1 Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800, by an endpoint-complete signed reciprocal-Poisson estimate for the eight-observation normal form.

Any candidate estimate must retain the hard-cutoff boundary, adjacent saddle recurrence, c_xi, correction derivative D, M_0=0, X_T+V=0, double contact, W_0=0, p=+-1, removable p=+-1/2, prime edges, square/cube transports, n=64, and q=1.

## Proof Boundary

This proves an exact six-moment/12-real flow matrix, rank at most eight, generic inertia (4,4,4), common-phase separation, three reciprocal stationary-phase skeletons, and the Hermitian reciprocal-diagonal null. It proves no Poisson summation remainder bound, hard-cutoff boundary estimate, near-diagonal bilinear gain, signed flow estimate, saddle-proxy or Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_flow_matrix_phase_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_flow_matrix_phase_reduction.py
```
