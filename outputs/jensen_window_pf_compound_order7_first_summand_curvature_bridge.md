# Jensen-Window PF Compound Order-Seven First-Summand Curvature Bridge

Date: 2026-07-13

Status: exact first/full-kernel transfer with one open continuous
first-summand theorem. This is not a proof of order-seven entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order7_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order7_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_first_summand_curvature_bridge.py
```

## Fourth Stable Layer

```text
T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1); r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))
r''=2*p''-q''+phi(T)*T''-chi(T)*(T')^2; phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2
R_k=r(k-1)-2*r(k)+r(k+1)=integral_[-1,1](1-|s|)*r''(k+s) ds
```

For `j>=322`, the completed first and full order-six curvature bounds
and the first/full defect floors give

```text
min(T_j,T_j^(1))>=3/(2*j), j>=320
first shifted numerator: 4*m^2+1769*m+154480>0,
full shifted numerator: 101*m^2+34869*m+740684>0.
```

The two missing floor indices are rigorous finite compositions:

```text
j=320: full margin=0.00307952557376919261108880603829, first margin=0.00307653251939570966771523643548,
j=321: full margin=0.00307619327493849500430945284960, first margin=0.00307323740019425827044763704280.
```

Here the full gap is `T_j=log(1+K_(j-5))`, bounded below by
`K/(1+K)`, and the first gap is within `O_j` of it.

## Full-Kernel Transfer

The rebalanced `2/k^8` moment-tail theorem propagates as

```text
a_j=2*((j-1)^(-8)+2*j^(-8)+(j+1)^(-8)); |B_j-B_j^(1)|<=a_j
L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j
U_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); |J_j-J_j^(1)|<=U_j
V_j=2*L_j+8*j*U_j; |h_j-h_j^(1)|<=V_j
W_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); |R_j-R_j^(1)|<=W_j
E_j=2*V_j+(5*j/7)*W_j+L_j; |q_j-q_j^(1)|<=E_j
Z_j=4*a_j+E_(j-1)+2*E_j+E_(j+1); |S_j-S_j^(1)|<=Z_j
Y_j=2*E_j+V_j+(2*j/3)*Z_j; |p_j-p_j^(1)|<=Y_j
O_j=5*a_j+Y_(j-1)+2*Y_j+Y_(j+1); |T_j-T_j^(1)|<=O_j
N_j=2*Y_j+E_j+(2*j/3)*O_j; |r_j-r_j^(1)|<=N_j
|R_k-R_k^(1)|<=N_(k-1)+2*N_k+N_(k+1)<262/k^2, k>=321
```

At the splice,

```text
321^2*transfer_error=261.334434117711045553258857298
reserve below 262=0.665565882288954446741142702026
```

The degree-102 shifted numerator has 103 strictly positive
coefficients, so the transfer bound holds on the entire half-line.

## Remaining Continuous Theorem

```text
r_1''(t)<=600/t^2 for every real t>=320
r_1''(t)<=600/t^2 => R_k^(1)<=600*[-log(1-1/k^2)]<601/k^2, k>=321
r_1''(t)<=600/t^2 on t>=320 => R_k<863/k^2<900/k^2 for k>=321
r_1''(t)<=600/t^2 on t>=320 => Q_(7,n)(-100)>0 for every n>=315
```

Thus the continuous theorem would contribute fewer than 601 inverse
squares, the full-kernel transfer fewer than 262, and the total 863
fits strictly inside the endpoint target 900. The continuous theorem
is now proved through `t=V'(2)`; its finite and asymptotic mode rays remain
open.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md
outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md
outputs/formal_core.md
```
