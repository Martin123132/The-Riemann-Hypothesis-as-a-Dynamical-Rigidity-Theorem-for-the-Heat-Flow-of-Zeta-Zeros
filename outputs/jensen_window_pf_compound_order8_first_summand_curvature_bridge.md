# Jensen-Window PF Compound Order-Eight First-Summand Curvature Bridge

Date: 2026-07-13

Status: exact first/full-kernel transfer with one open continuous
first-summand theorem. This is not a proof of order-eight entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_first_summand_curvature_bridge.py
```

## Fifth Stable Layer

```text
U(t)=6*B(t)-r(t-1)+2*r(t)-r(t+1); s(t)=2*r(t)-p(t)+log(1-exp(-U(t)))
s''=2*r''-p''+phi(U)*U''-chi(U)*(U')^2; phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2
W_k=s(k-1)-2*s(k)+s(k+1)=integral_[-1,1](1-|v|)*s''(k+v) dv
```

The completed order-seven first/full bounds prove

```text
min(U_j,U_j^(1))>=3/(2*j), j>=1249
first shifted numerator: 6*m**2 + 12581*m + 6352461>0,
full shifted numerator: 756*m**2 + 1456613*m + 639733131>0.
```

## Full-Kernel Transfer

The rebalanced 2/k^9 moment-tail theorem propagates as

```text
a_j=2*((j-1)^(-9)+2*j^(-9)+(j+1)^(-9)); |B_j-B_j^(1)|<=a_j
L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j
U1_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); |J_j-J_j^(1)|<=U1_j
V_j=2*L_j+8*j*U1_j; |h_j-h_j^(1)|<=V_j
W1_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); |R_j-R_j^(1)|<=W1_j
E_j=2*V_j+(5*j/7)*W1_j+L_j; |q_j-q_j^(1)|<=E_j
Z_j=4*a_j+E_(j-1)+2*E_j+E_(j+1); |S_j-S_j^(1)|<=Z_j
Y_j=2*E_j+V_j+(2*j/3)*Z_j; |p_j-p_j^(1)|<=Y_j
O_j=5*a_j+Y_(j-1)+2*Y_j+Y_(j+1); |T_j-T_j^(1)|<=O_j
N_j=2*Y_j+E_j+(2*j/3)*O_j; |r_j-r_j^(1)|<=N_j
P_j=6*a_j+N_(j-1)+2*N_j+N_(j+1); |U_j-U_j^(1)|<=P_j
C_j=2*N_j+Y_j+(2*j/3)*P_j; |s_j-s_j^(1)|<=C_j
|W_k-W_k^(1)|<=C_(k-1)+2*C_k+C_(k+1)<190/k^2, k>=1250
```

At the splice,

```text
1250^2*transfer_error=177.991946106896177163915096431
reserve below 190=12.0080538931038228360849035688
```

The degree-133 shifted numerator has 134 strictly positive
coefficients, so the transfer bound holds on the complete half-line.

## Remaining Continuous Theorem

```text
s_1''(t)<=4000/t^2 for every real t>=999
s_1''(t)<=4000/t^2 => W_k^(1)<=4000*[-log(1-1/k^2)]<4001/k^2, k>=1250
s_1''(t)<=4000/t^2 on t>=999 => W_k<4001/k^2+190/k^2=4191/k^2<4300/k^2, k>=1250
s_1''(t)<=4000/t^2 on t>=999 => Q_(8,n)(-100)>0 for every n>=1243
```

Thus the continuous theorem would contribute fewer than 4001 inverse
squares, the full-kernel transfer fewer than 190, and the total 4191
fits strictly inside the endpoint target 4300. The continuous theorem
itself remains open.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.md
outputs/formal_core.md
```
