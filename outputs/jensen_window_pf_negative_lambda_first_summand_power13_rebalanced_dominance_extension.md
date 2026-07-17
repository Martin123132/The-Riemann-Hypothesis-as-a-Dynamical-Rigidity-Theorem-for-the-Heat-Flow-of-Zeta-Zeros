# Jensen-Window PF Rebalanced Power-Thirteen First-Summand Dominance

Date: 2026-07-16

Status: rigorous inverse-thirteenth-power first-summand dominance at
`lambda=-100` for `k>=340`. This is not a proof of order eleven,
PF-infinity, RH, or `Lambda<=0`.

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-13), k>=340
epsilon(0)*P_k(u<a(k))<k^(-13), k>=340
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^13 for every integer k>=340
0<=log(1+delta_k)<2/k^13, k>=340
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-13)+2*j^(-13)+(j+1)^(-13)), j>=341
```

At `k=340`, all fourteen Arb endpoint and derivative gates are strict.
Exact monotonicity then covers every larger integer. This is the
natural perturbation power for an eighth stable-log transfer; that
order-eleven transfer remains a separate theorem.
