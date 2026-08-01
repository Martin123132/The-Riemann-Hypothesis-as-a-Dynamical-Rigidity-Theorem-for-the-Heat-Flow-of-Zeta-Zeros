# Terminal-Tail Explicit Five-Moment Phi Quadratic Reduction

Date: 2026-07-31

Status: exact rank-four quadratic reduction with `0 Phi_B bounds`. This is
not a proof artifact; the signed arithmetic estimate and RH remain open.

## Five-Moment Coefficients

Fix the q=1, L>=50, fixed-N chart and admissible growing terminal block of the source reductions. Put B=N-M-1 and differentiate with N,M,B fixed. All G_r stop at B and use the same eta/S_a normalization as the tail.

```text
Let L=log(N), G=(G_0,...,G_4)^T. Then v_0=(1,rho_1,rho_2,0,0), v_1=(L,Lrho_1-1,Lrho_2-rho_1,-rho_2,0), v_2=(L^2,L^2rho_1-2L,L^2rho_2-2Lrho_1+1,rho_1-2Lrho_2,rho_2), e_0=(0,rho_(1,x),rho_(2,x),0,0), and e_1=(0,Lrho_(1,x),Lrho_(2,x)-rho_(1,x),-rho_(2,x),0).

q_0=chi_N*v_0+s_*'*v_1+e_0; q_1=chi_N*v_1+s_*'*v_2+e_1; a_0=s_*'*v_1+i*b*u_N*v_0; n_0=s_*''*v_1+s_*'*q_1+i*(b_x*u_N+b*u_(N,x))*v_0+i*b*u_N*q_0.
```

The five base vectors convert the correction-free moments to K_0,K_1,K_2,E_0,E_1. The four derived vectors then give the bulk value, value derivative, centered scalar, and centered-scalar derivative without storing a prefix path.

Define the four real observations by

```text
V=Re(v_0 dot G)=X_B, Q=Re(q_0 dot G)=X_(B,x), A=Re(a_0 dot G)=A_B, and N=Re(n_0 dot G)=A_(B,x).
```

## Exact Joined Form

```text
h^2*Phi_B=(X_T+V)N+V*A_(T,x)-(A_T+A)Q-A*X_(T,x).

For g=(Re G_0,...,Re G_4,Im G_0,...,Im G_4)^T and R(z)=(Re z_0,...,Re z_4,-Im z_0,...,-Im z_4)^T, put v=R(v_0), q=R(q_0), a=R(a_0), n=R(n_0). Then h^2*Phi_B=g^T M g+ell_T^T g, M=(v n^T+n v^T-a q^T-q a^T)/2 and ell_T=X_T n+A_(T,x)v-A_T q-X_(T,x)a.
```

Although five complex moments have ten real components, Phi_B depends on them only through four real observations. The exact quadratic matrix has rank at most four and its tail-linear term lies in the same four-row span.

The exact factorization is

```text
M=[v n a q] J [v n a q]^T with J=(1/2)*diag_offdiag(+1,-1); hence rank(M)<=4. The observation map g->(V,Q,A,N) has a real kernel of dimension at least 6, and Phi_B vanishes on that kernel.
```

## Indefiniteness And Fibres

An exact algebraic specialization L=2,rho_1=rho_2=rho_(1,x)=rho_(2,x)=0,s_*'=-i/2,s_*''=0,chi_N=1,u_N=u_(N,x)=1 has observation rank 4 and minor 5/16. By the displayed congruence its M has inertia (2,2,6). This is a generic nonpromotion guard, not a physical Xi state.

This guard forbids a generic positivity or negativity promotion. It does not
replace analysis of the actual Xi moment vector.

```text
On X_T+V=0, h^2*Phi_B=V*A_(T,x)-(A_T+A)Q-A*X_(T,x), with no division by the total value.

If also A_T+A=0, then Phi_B=-P_T>1/400. Therefore the desired bound must source-specifically exclude a simultaneous retained value/centered-scalar contact; it cannot assume it away.
```

The second identity explains why the zero fibre cannot be discarded: the
desired inequality must itself prevent the simultaneous contact.

## Open Arithmetic Target

```text
For the actual Xi moment vector g_B and certified tail jets, prove g_B^T M g_B+ell_T^T g_B<=h^2/400, preferably <=h^2/800.
```

Insert the exact balanced Mangoldt or Vaughan decomposition into the four observations V,Q,A,N, or equivalently expand the rank-four form as its original two-carrier kernel. Keep the endpoint, complete square, both wings, and all transpose and Hermitian cross terms before taking absolute values.

Any proposed signed criterion must survive X_T+V=0, W_0=0, zero observation rows, p=+-1, removable p=+-1/2, prime edges, square and cube transports, n=64, q=1, and adjacent cutoffs. The simultaneous value/scalar contact is an obstruction to be excluded by the theorem, not an admissible denominator.

## Boundary

No new `pi` occurs in the coefficient vectors or realification. Existing
`pi` factors retain their completed-zeta/Riemann--Siegel origin.

This proves an exact explicit rank-four real quadratic normal form, fibre identities, and a generic indefiniteness guard. It proves no upper bound for Phi_B, endpoint-coupled Type-I/II or Vaughan estimate, contact exclusion, retained aggregate sign, Xi residual transfer, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_explicit_phi_quadratic_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_explicit_phi_quadratic_reduction.py
```
