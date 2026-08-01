# Contiguous Terminal-Tail Anchored Abel Cross-Current Reduction

Date: 2026-07-31

Status: exact cancellation-preserving reduction to one signed nonterminal
correlation, with `0 bulk closures`. This is not a proof artifact; the
`Phi_B` bound and RH remain open.

## Fixed-Chart Partition

Fix a physical q=1, L>=50, fixed-N chart and an admissible M>=1 with h*(M+1)^18<=1. Put B=N-M-1>=1 and keep M,B fixed while differentiating inside the chart.

```text
T=e+sum_(n=B+1)^N z_n, G_T=g+s_*'*sum_(n=B+1)^N u_n*z_n, F_B=sum_(n=1)^B z_n, U_B=sum_(n=1)^B u_n*z_n; W_0=T+F_B and W_A=G_T+s_*'*U_B.
U_B=u_N*F_B+R_B, with R_B=sum_(n=1)^B log(N/n)z_n.
```

The terminal block is exactly the endpoint-terminal edge together with the
carriers `N-M,...,N`; the bulk is `1,...,B`.

## Truncated Abel Identity

```text
F_k=sum_(n=1)^k z_n, h_k=log((k+1)/k), kappa_B=u_B-u_N=log(N/B); R_B=sum_(n=1)^B(u_n-u_N)z_n=kappa_B*F_B+sum_(k=1)^(B-1)h_k*F_k
Because N and B are fixed inside the chart, partial_x(u_n-u_N)=0; hence R_(B,x)=kappa_B*F_(B,x)+sum_(k=1)^(B-1)h_k*F_(k,x)
```

The differences `u_n-u_N=log(N/n)`, `kappa_B`, and every `h_k` are constant
under `x` differentiation while `N` and `B` are fixed. Thus no derivative of
an Abel weight appears.

## Tail-Anchored Centering

Use the same real positive source scale `S_a` and the same `u_N` shear as the
growing-prefix theorem:

```text
hat(T)=T/S_a, hat(F_k)=F_k/S_a, hat(R_B)=R_B/S_a; hat(W_0)=hat(T)+hat(F_B)
A_T=Re(hat(G_T))-c*u_N*Re(hat(T))
A_B=Re(s_*'*hat(R_B))-b*u_N*Im(hat(F_B))=c*Re(hat(R_B))-b*Im(hat(R_B))-b*u_N*Im(hat(F_B))
A=A_T+A_B=Re(hat(G_T)-s_*'*u_N*hat(T))-b*u_N*Im(hat(W_0))+Re(s_*'*hat(R_B))
```

In particular, the certified terminal block replaces the old isolated
endpoint without recentering the problem at `u_B`.

The differentiated bulk scalar is

```text
A_(B,x)=c_x*R_R+c*R_(R,x)-b_x*R_I-b*R_(I,x)-(b_x*u_N+b*u_(N,x))*F_I-b*u_N*F_(I,x), where F=hat(F_B), R=hat(R_B)
```

## Four Coordinates and Current

```text
For a block with normalized real value X and normalized centered scalar A, (C,D,E,N)=(X,A/h,X_x/h,A_x/h^2) and P=CN-DE
If X=C_0/S_a and A=A_0/S_a, then P=[C_0*A_(0,x)-A_0*C_(0,x)]/(h^2*S_a^2)=a^2*J/S_a^2; every S_(a,x) term cancels
```

Apply this map separately to the tail and bulk. Exact polarization gives

```text
P_ret=P_T+P_B+P_cross, P_B=C_B*N_B-D_B*E_B, P_cross=C_T*N_B+C_B*N_T-D_T*E_B-D_B*E_T
```

Equivalently, all nonterminal uncertainty is

```text
Phi_B=P_B+P_cross=h^-2*[(X_T+X_B)A_(B,x)+X_B*A_(T,x)-(A_T+A_B)X_(B,x)-A_B*X_(T,x)]
```

Compute hat(F_k) and hat(F_(k,x)) in index order, accumulate hat(R_B)=kappa_B*hat(F_B)+sum h_k*hat(F_k) and its derivative, then evaluate Phi_B in constant additional work. This is O(B) time and O(1) streaming state, rather than an O(B^2) pair expansion.

## Sharp Remaining Target

```text
The imported growing-prefix theorem supplies P_T=a^2*J_(tail,M)/S_a^2<-1/400.
The signed inequality Phi_B<=1/400 implies P_ret<0. The stronger Phi_B<=1/800 preserves P_ret<-1/800.
```

This is a signed upper-bound problem. Replacing `Phi_B` by separate absolute
pair bounds discards the prefix cancellation retained by the Abel form.

The reason separate clockwise theorems are insufficient is exact:

```text
At one point v_T=(1,0), v_(T,x)=(0,-1), v_B=(-1/2,0), v_(B,x)=(0,2): P_T=P_B=-1, P_cross=5/2, but P_ret=1/2
```

The certified terminal block now replaces the isolated endpoint in the Abel identity without changing the u_N shear. All nonterminal uncertainty is concentrated in one signed cumulative correlation Phi_B.

## Boundary

Do not recenter at u_B: u_N is the centering used by the certified current. An x-dependent change of shear contributes an extra shear-derivative current unless it is tracked.

The only pi occurs through the inherited completed-zeta/Riemann-Siegel normalization and h=1/a. The Abel weights are logarithmic arithmetic ratios and introduce no fitted geometry.

Prove the Xi-specific signed upper bound Phi_B<=1/400, preferably Phi_B<=1/800, from the actual multiplicative carrier phases and cumulative prefixes before taking absolute values.

This is an exact cancellation-preserving reduction, not the required bound. It proves no bulk aggregate closure, complete retained current sign, Xi-level current theorem, Abel gap, winding cap, contact exclusion, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.

## Reproduce

```text
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_abel_cross_current_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_abel_cross_current_reduction.py
```
