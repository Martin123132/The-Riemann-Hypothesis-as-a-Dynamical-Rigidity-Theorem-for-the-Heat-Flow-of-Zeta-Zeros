# Jensen-Window PF Compound Order-Nine First-Summand Curvature Bridge

Date: 2026-07-13

Status: exact first/full-kernel transfer with one open continuous
first-summand theorem and a two-index finite splice. This is not a
proof of order-nine entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order9_first_summand_curvature_bridge.json
python work/rh_compute/scripts/jensen_window_pf_compound_order9_first_summand_curvature_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_first_summand_curvature_bridge.py
```

## Sixth Stable Layer

```text
V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1); w(t)=2*s(t)-r(t)+log(1-exp(-V(t)))
w''=2*s''-r''+phi(V)*V''-chi(V)*(V')^2; phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2
Y_k=w(k-1)-2*w(k)+w(k+1)=integral_[-1,1](1-|v|)*w''(k+v) dv
min(V_j,V_j^(1))>=4/(3*j), j>=1250
phi(min(V_j,V_j^(1)))<=1/min(V_j,V_j^(1))<=3*j/4
first shifted numerator: 13*m**2 + 17490*m + 1542497>0,
full shifted numerator: 3271*m**2 + 4140000*m + 62044250>0.
```

## Full-Kernel Transfer

The rebalanced `2/k^10` moment-tail theorem propagates through all
six stable logarithms. The final two stages are

```text
D_j=7*a_j+C_(j-1)+2*C_j+C_(j+1); |V_j-V_j^(1)|<=D_j
F_j=2*C_j+N_j+(3*j/4)*D_j; |w_j-w_j^(1)|<=F_j
|Y_k-Y_k^(1)|<=F_(k-1)+2*F_k+F_(k+1)<550/k^2, k>=1251
1251^2*transfer_error=533.841639334312505775394634804
reserve below 550=16.1583606656874942246053651960
```

The degree-168 shifted numerator has 169 strictly positive
coefficients, so the transfer bound holds on the complete half-line.

## Remaining Targets

```text
w_1''(t)<=4200/t^2 for every real t>=1250
w_1''(t)<=4200/t^2 => Y_k^(1)<=4200*[-log(1-1/k^2)]<4201/k^2, k>=1251
w_1''(t)<=4200/t^2 on t>=1250 => Y_k<4201/k^2+550/k^2=4751/k^2<4900/k^2, k>=1251
w_1''(t)<=4200/t^2 on t>=1250 => Q_(9,n)(-100)>0 for every n>=1243
prove Q_(9,n)(-100)>0 for n=1241,1242
```

Thus the continuous theorem would contribute fewer than 4201 inverse
squares and the full-kernel transfer fewer than 550. Their total 4751
fits strictly inside 4900. The continuous theorem and the two finite
signs remain open.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.md
outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md
```
