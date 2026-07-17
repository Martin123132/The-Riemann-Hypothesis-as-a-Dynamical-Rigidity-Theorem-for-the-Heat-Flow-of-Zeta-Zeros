# Jensen-Window PF Rebalanced Power-Twelve First-Summand Dominance

Date: 2026-07-16

Status: rigorous inverse-twelfth-power first-summand dominance at
`lambda=-100` for `k>=320`. This is not a proof of order ten,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power12_rebalanced_dominance_extension.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power12_rebalanced_dominance_extension.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power12_rebalanced_dominance_extension.py
```

## Reused Split

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
```

At `k=320`, all fourteen Arb endpoint and derivative gates are strict.
The same exact monotonicity identities then cover every larger index.

```text
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-12), k>=320
epsilon(0)*P_k(u<a(k))<k^(-12), k>=320
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^12 for every integer k>=320
0<=log(1+delta_k)<2/k^12, k>=320
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-12)+2*j^(-12)+(j+1)^(-12)), j>=321
```

This is the perturbation power used by the order-ten seventh-gap
full-kernel transfer. That transfer is a separate theorem.
