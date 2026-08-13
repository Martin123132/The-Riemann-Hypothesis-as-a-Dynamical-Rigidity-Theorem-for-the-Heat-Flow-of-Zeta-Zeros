# Physical P_lin Ideal-Cubic Fresnel Reciprocal Gate

Date: 2026-08-02

Status: exact weighted off-saddle transport, active-roster Fubini reindexing, reciprocal joint-phase self-duality, second-order Abel identity, and one absolute-family barrier. This is not a proof of a grouped h^2 estimate, signed flow, or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_gate.py
```

## Weighted Off-Saddle Jet

```text
Z(x)=exp(S_tilde(x))Y(x), with Z'=[Dhat+g_tilde(x)I_8]Z.
A_z=zJ^2/v, B_z=zJJ'-J, and K_z(epsilon)=A_z[Dhat+(-1/2+epsilon)I_8]+B_zI_8.
H_r(y)=[A_z(Dhat+g_vI_8)+B_zI_8]Z(x_r)+Z(x_(0,r)); c'_(bar P,r)(y)=exp(S(ell))wH_r(y)/(r sqrt(alpha_P)z^2), where w=(alpha^T,i beta^T).
```

## Exact Roster Transport

```text
P_d(epsilon)=exp[-d(1/2+epsilon)+td^2/4]E(d), so r/(r+1) times the shifted weighted jet equals P_d(epsilon)Z.
P_(d_2)(epsilon-td_1/2)P_(d_1)(epsilon)=P_(d_1+d_2)(epsilon) exactly.
L_d=[K_z(epsilon_v-td/2)P_d(epsilon_v),P_d(epsilon_0)], L_0=[K_z(epsilon_v),I_8].
Delta c'_r(y)=exp(S(ell))w[L_(d_r)-L_0]V_r/[r sqrt(alpha_P)z^2] at fixed y, V_r=(Z(x_r),Z(x_(0,r)))^T.
Delta^2c'_r(y)=exp(S(ell))w[L_(d_r+d_(r+1))-2L_(d_r)+L_0]V_r/[r sqrt(alpha_P)z^2] at fixed y.
At z=0 replace L_d by S_d=M(epsilon_0-td/2)P_d(epsilon_0); the exact first and second operators are S_d-S_0 and S_(d_1+d_2)-2S_(d_1)+S_0.
```

These are fixed-y identities. They retain the moving exponential, the 1/r denominator, the epsilon shift, the field term, and the saddle anchor in one block operator.

## Active Roster

```text
For fixed y and v=v(y), the active modes are T_y=T intersect {ceil(alpha_P v/B),...,floor(alpha_P v)}.
sum_(r in T)e(phi_r)integral_(y_(1,r))^(y_(B,r))f_r(y)dy=integral_R sum_(r in T_y)e(phi_r)f_r(y)dy exactly.
The moving integration limits become the two integer jumps of one contiguous active roster; no endpoint strip is discarded.
```

## Reciprocal Self-Duality

```text
For integer q,r, e(alpha_P log u-ru)=e(Psi_q(r,u)), Psi_q=alpha_P log u-r(u-q).
The unique positive continuous critical point is (r_*,u_*)=(alpha_P/q,q).
Psi_q(r_*,u_*)=alpha_P log q.
Hess_(r,u)Psi_q=[[0,-1],[-1,-alpha_P/q^2]], with determinant -1 and inertia (1,1).
The sequential u and r curvatures are -alpha_P/q^2 and q^2/alpha_P; their product is -1 and their Fresnel signatures cancel.
```

The reciprocal r-saddle returns the original carrier phase e(alpha_P log q), so the Morse roster is self-dual rather than an unrelated error lattice.

The determinant -1 is exact. The continuous r saddle need not be an integer; it is the legitimate stationary coordinate for the lattice phase, just as the continuous u saddle is used in the Morse chart.

## Reciprocal Scale

```text
At r_*=alpha_P/q, x_q=log(q/a), the ideal saddle coefficient is q exp(S(log q))wM(epsilon_q)Y(x_q)/(2alpha_P^(3/2)).
The reciprocal r-curvature width is sqrt(alpha_P)/q.
Their exact scalar product is exp(S(log q))/(2alpha_P) before the joined row and saddle matrix.
```

## Second-Order Abel And Barrier

```text
For S_j=sum_(k=m)^jz_k and T_j=sum_(l=m)^jS_l, sum_(r=m)^n w_rz_r=w_nS_n-Delta w_(n-1)T_(n-1)+sum_(r=m)^(n-2)Delta^2w_rT_r.
For bar P=1, the saddle scalar is m_q=epsilon_q^2+t/2-1/12 and the physical box gives |m_q|>1/13.
Since t>=0 and sigma<201/400, exp(S(log q))>=q^(-201/400), with strict inequality for q>1.
If A_B is the nonnegative sum of the reciprocal saddle coefficient scale times its curvature width for bar P=1, then A_B>200(B^(199/400)-1)/(2587alpha_P).
```

The factor B^(199/400)-1 is unbounded. Therefore summing the reciprocal q-family saddle scales by absolute values cannot prove a uniform C/alpha_P grouped estimate. This is a lower bound only on that chosen nonnegative scale budget, not on the signed Morse remainder.

## Handoff

Use the exact fixed-y L_d transport inside local r-cells, but rejoin the reciprocal stationary contributions with mathcal F, the physical q-carriers, both endpoint packages, and the terminal recurrence before summing over q. A second cancellation across the q-family is compulsory.

The fixed-y O(alpha_P^-2) field transport remains genuine, but it is only the inner half of a two-level problem. The reciprocal r saddles carry the original q phases. The next theorem must identify the cancellation between the transformed Morse carrier, the c-prime correction, and the endpoint-complete physical carrier before estimating the q sum.

## Proof Boundary

This gate proves exact weighted off-saddle and saddle mode transports, finite active-roster Fubini reindexing, reciprocal joint-phase geometry, one second-order Abel identity, and one nonnegative absolute-family barrier. It proves no grouped ideal-cubic or endpoint-composed h2 estimate, reciprocal q-family cancellation theorem, quadratic residual bound, signed flow estimate, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
