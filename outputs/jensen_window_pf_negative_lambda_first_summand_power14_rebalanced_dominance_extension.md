# Jensen-Window PF Rebalanced Power-Fourteen First-Summand Dominance

Date: 2026-07-22

Status: rigorous inverse-fourteenth-power first-summand dominance at
`lambda=-100` for `k>=380`. This is not a proof of order twelve,
PF-infinity, RH, or `Lambda<=0`.

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-14), k>=380
epsilon(0)*P_k(u<a(k))<k^(-14), k>=380
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^14 for every integer k>=380
0<=log(1+delta_k)<2/k^14, k>=380
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-14)+2*j^(-14)+(j+1)^(-14)), j>=381
```

At `k=380`, all fourteen Arb endpoint and derivative gates are strict.
Exact monotonicity then covers every larger integer. This is the
natural perturbation power for a ninth stable-log transfer; that
order-twelve transfer remains a separate theorem.
