# Jensen-Window PF Compound Order-Six First-Summand Curvature Bridge

Date: 2026-07-13

Status: exact first/full-kernel transfer with one open continuous
first-summand theorem. This is not a proof of order-six entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_first_summand_curvature_bridge.py
```

## Third Stable Layer

```text
S(t)=4*B(t)-q(t-1)+2*q(t)-q(t+1); p(t)=2*q(t)-h(t)+log(1-exp(-S(t)))
p''=2*q''-h''+phi(S)*S''-chi(S)*(S')^2; phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2
P_k=p(k-1)-2*p(k)+p(k+1)=integral_[-1,1](1-|s|)*p''(k+s) ds
```

The completed order-five curvature and defect bounds imply

```text
min(S_j,S_j^(1))>=3/(2*j), j>=321
first floor numerator after j=321+m: 2*m^2+1029*m+124101>0,
full floor numerator after j=321+m: 254*m^2+112693*m+9977039>0.
```

## Full-Kernel Transfer

The strengthened `2/k^7` moment-tail theorem propagates as

```text
a_j=2*((j-1)^(-7)+2*j^(-7)+(j+1)^(-7)); |B_j-B_j^(1)|<=a_j
L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j
U_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); |J_j-J_j^(1)|<=U_j
V_j=2*L_j+8*j*U_j; |h_j-h_j^(1)|<=V_j
W_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); |R_j-R_j^(1)|<=W_j
E_j=2*V_j+(5*j/7)*W_j+L_j; |q_j-q_j^(1)|<=E_j
Z_j=4*a_j+E_(j-1)+2*E_j+E_(j+1); |S_j-S_j^(1)|<=Z_j
Y_j=2*E_j+V_j+(2*j/3)*Z_j; |p_j-p_j^(1)|<=Y_j
|P_k-P_k^(1)|<=Y_(k-1)+2*Y_k+Y_(k+1)<100/k^2, k>=322
```

At the splice,

```text
322^2*transfer_error=97.4494017121864439638173083660
reserve below 100=2.55059828781355603618269163403
```

The degree-75 shifted numerator has 76 strictly positive
coefficients, so the transfer bound holds on the entire half-line.

## Remaining Continuous Theorem

```text
p_1''(t)<=200/t^2 for every real t>=321
p_1''(t)<=200/t^2 => P_k^(1)<=200*[-log(1-1/k^2)]<201/k^2, k>=322
p_1''(t)<=200/t^2 on t>=321 => P_k<301/k^2<320/k^2 for k>=322
```

Thus the continuous theorem supplies fewer than 201 inverse squares,
the full-kernel transfer supplies fewer than 100, and the total 301
fits strictly inside the endpoint target 320.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.md
outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md
outputs/formal_core.md
```
