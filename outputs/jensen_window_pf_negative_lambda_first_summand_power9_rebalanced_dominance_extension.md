# Jensen-Window PF Rebalanced Power-Nine First-Summand Dominance

Date: 2026-07-13

Status: rigorous inverse-ninth-power first-summand dominance at
`lambda=-100`. This is not a proof of order eight, PF-infinity, RH,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.py
```

## Reused Split

The exact summand-ratio monotonicity, epsilon(0) cap, and split remain

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
```

The previous endpoint inequalities retain more than one additional
log(k) of reserve at k=300.

## Endpoint Gates

At k=300, outward-rounded Arb arithmetic gives

```text
high log comparison=[-38.11530344306961937941322954918953557539736335793766074672687222084081 +/- 4.04E-69]<0,
low log comparison=[-16.75045421515097096127501822929304465099741704966740075102521837149116 +/- 1.67E-69]<0,
high derivative margin=2.79130236236126579978175269934735796758531232744984256672944E+1>0,
low power-derivative margin=2.33664016042713732356706005808783621003054381058540423959103E-1>0.
```

All fourteen endpoint and derivative gates remain strictly positive,
and the same six exact ratio-log derivatives cover the full half-line.

## Strengthened Tail

```text
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-9), k>=300
epsilon(0)*P_k(u<a(k))<k^(-9), k>=300
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^9 for every integer k>=300
0<=log(1+delta_k)<2/k^9, k>=300
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-9)+2*j^(-9)+(j+1)^(-9)), j>=301
```

This is the power needed to pass through the fifth stable logarithm
while retaining an inverse-square order-eight curvature transfer.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.md
outputs/formal_core.md
```
