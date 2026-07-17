# Jensen-Window PF Rebalanced Power-Ten First-Summand Dominance

Date: 2026-07-13

Status: rigorous inverse-tenth-power first-summand dominance at
`lambda=-100`. This is not a proof of order nine, PF-infinity, RH,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.py
```

## Reused Split

The exact summand-ratio monotonicity, epsilon(0) cap, and split remain

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
```

## Endpoint Gates

At k=300, outward-rounded Arb arithmetic gives

```text
high log comparison=[-32.41152096841341831998200140289828145554766982285736534292552208526817 +/- 3.55E-69]<0,
low log comparison=[-11.04667174049476990184379008300179053114772351458710534722386823591852 +/- 1.18E-69]<0,
high derivative margin=2.69130236236126579978175269934735796758531232744984256672944E+1>0,
low power-derivative margin=2.33664016042713732356706005808783621003054381058540423959103E-1>0.
```

All fourteen endpoint and derivative gates remain strictly positive,
and the same six exact ratio-log derivatives cover the full half-line.

## Strengthened Tail

```text
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-10), k>=300
epsilon(0)*P_k(u<a(k))<k^(-10), k>=300
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^10 for every integer k>=300
0<=log(1+delta_k)<2/k^10, k>=300
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-10)+2*j^(-10)+(j+1)^(-10)), j>=301
```

This is the perturbation power needed to pass through the sixth stable
logarithm while retaining an inverse-square order-nine curvature transfer.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.md
outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md
outputs/formal_core.md
```
