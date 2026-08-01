# Terminal-Tail Five-Moment Two-Carrier Kernel Reduction

Date: 2026-07-31

Status: exact endpoint/Hermitian/transpose reduction with `0 signed Type-I/II
bounds`. This is not a proof artifact; `Phi_B` and RH remain open.

## Carrier Polynomials

Use the fixed physical q=1, L>=50 terminal/bulk chart with B=N-M-1. Write ell_n=log(n), lambda_n=log(N/n), and use the correction-free normalized carrier w_n in G_r=sum ell_n^r w_n. The centered bulk rate is chi_N and is distinct from the established endpoint factor kappa_N.

```text
For ell=log(n), lambda=L-ell=log(N/n), put C(ell)=1+rho_1*ell+rho_2*ell^2, D(ell)=rho_(1,x)*ell+rho_(2,x)*ell^2, R(ell)=s_*'*lambda+i*b*u_N, Q(ell)=(chi_N+s_*'*lambda)C(ell)+D(ell), R_x(ell)=s_*''*lambda+i(b_x*u_N+b*u_(N,x)), and N(ell)=R(ell)Q(ell)+R_x(ell)C(ell).
```

```text
For w_n=(eta/S_a)exp[t(log n)^2/4-s_*log n], mathscr V=sum C_nw_n, mathscr Q=sum Q_nw_n, mathscr A=sum R_nC_nw_n, mathscr N=sum N_nw_n; the real observations are their real parts.
```

The exact radial-subtracted defect is

```text
R(ell)C(ell)-Q(ell)=[i*b*u_N-chi_N]C(ell)-D(ell).
```

## Endpoint And Pair Kernels

```text
L_n=X_T*N_n+A_(T,x)C_n-A_TQ_n-X_(T,x)R_nC_n=(X_TR_n-A_T)Q_n+[X_TR_(n,x)+A_(T,x)-X_(T,x)R_n]C_n.

H_(n,m)=C_n*conj(N_m)-R_nC_n*conj(Q_m)=C_n[(conj(R_m)-R_n)conj(Q_m)+conj(R_(m,x))conj(C_m)]; T_(n,m)=C_nN_m-R_nC_nQ_m=C_n[(R_m-R_n)Q_m+R_(m,x)C_m].

h^2*Phi_B=Re sum_(n<=B)L_nw_n+(1/2)Re sum_(n,m<=B)[H_(n,m)w_nconj(w_m)+T_(n,m)w_nw_m].
```

The two exact symmetrizations are

```text
H^+_(n,m)=(H_(n,m)+conj(H_(m,n)))/2=(1/2){[conj(R_m)-R_n][C_nconj(Q_m)-conj(C_m)Q_n]+C_nconj(C_m)[conj(R_(m,x))+R_(n,x)]}.

For Delta=lambda_m-lambda_n, T^+_(n,m)=(T_(n,m)+T_(m,n))/2=(1/2){s_*'^2Delta^2C_nC_m+s_*'Delta(C_nD_m-C_mD_n)+C_nC_m[R_(n,x)+R_(m,x)]}.
```

The Hermitian kernel multiplies w_n conjugate(w_m) and carries phase differences. The transpose kernel multiplies w_nw_m and carries phase sums. Both occur in the real projected current; controlling only one correlation family is incomplete.

T^+ is independent of chi_N. In H^+, C_nconj(Q_m)-conj(C_m)Q_n=C_nconj(C_m)[conj(chi_N)-chi_N+conj(s_*')lambda_m-s_*'lambda_n]+C_nconj(D_m)-conj(C_m)D_n, so Re(chi_N) cancels.

Expanding the rank-four moment quadratic directly produces one ordered two-carrier sum, not a four-variable product of two separately absolute-valued Mangoldt sums. The degree-four polynomials may be rejoined with the exact divisor identities while preserving both oscillatory kernels.

## Physical Anchor

```text
At n=N, hat(z_N)=C_Nw_N=-Q(p)exp(W_theta), Q_x/Q=-i*h*theta/2, and C_(N,x)/C_N=delta_N; hence chi_N=-i*h*theta/2+W_(theta,x)-delta_N.

Since delta_N=d_(N,x)/(1+d_N)-d_(1,x)/(1+d_1), the bounds |d_n|<1/2 and |d_(n,x)|<4223/x^2 give |delta_N|<16892/x^2. Together with |W_(theta,x)|<2h^2, x>a^2=h^-2, and 16892h^2<1, this gives |chi_N+i*h*theta/2|<3h^2.

For u_N=-log(1-h*theta), |s_*'+i/2|<h^2/6000 and 0<=u_N-h*theta<h^2 imply |i*b*u_N-chi_N|<4h^2. The full coefficient defect also contains -D(ell), which must remain in the bulk estimate.
```

The terminal phase anchor makes chi_N=O(h), with an O(h^2) error around -i*h*theta/2, and makes i*b*u_N-chi_N=O(h^2). Thus the real observations A and Q have a source-specific near-coincidence. This is genuine structure absent from the generic inertia witness, but no bound on the summed D(ell) defect or the two oscillatory kernels is yet proved.

## Leading Kernel And Guard

```text
At C=1,D=0,s_*'=-i/2,s_*''=0,b=-1/2,b_x=0,chi_N=i*b*u_N, one has R_n=Q_n=-i*u_n/2 with u_n=lambda_n+u_N, H^+_(n,m)=-(u_n+u_m)^2/8, and T^+_(n,m)=-(lambda_n-lambda_m)^2/8-i*u_(N,x)/2.
```

```text
u_1=h,w_1=1 gives P_bulk=-h^2/4
u_1=h,u_2=2h,w_1=1,w_2=-1/2 gives P_bulk=h^2/8
```

In the ideal q=1 leading kernel every single positive-distance carrier is clockwise, but two real opposite-phase carriers with distances h and 2h and amplitudes 1 and 1/2 give positive bulk current h^2/8. This is an algebraic kernel guard, not an attained Xi configuration or aggregate counterexample.

## Open Arithmetic Target

```text
Prove that the endpoint-linear sum plus half the real part of the complete Hermitian and transpose two-carrier sums is at most h^2/400, preferably h^2/800, for the actual Xi carriers.
```

Apply balanced Type-I/II or Vaughan decomposition jointly to the phase-difference and phase-sum kernels. Retain the complete square, both wings, correction derivative D(ell), terminal anchor, and adjacent-cutoff transport. A theorem for only the transpose-symmetric Mangoldt form cannot close Phi_B.

## Boundary

No new `pi` occurs in the kernel decomposition. The `pi` in
`Q_x/Q=-i*h*theta/2` comes from the inherited Riemann--Siegel terminal phase.

This proves an exact endpoint-linear plus Hermitian/transpose two-carrier decomposition, common-rate cancellations, the physical chi_N anchor bound, leading kernels, and a pair-sign guard. It proves no signed Type-I/II or Vaughan estimate, upper bound on Phi_B, contact exclusion, retained aggregate sign, Xi residual transfer, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_two_carrier_kernel_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_two_carrier_kernel_reduction.py
```
