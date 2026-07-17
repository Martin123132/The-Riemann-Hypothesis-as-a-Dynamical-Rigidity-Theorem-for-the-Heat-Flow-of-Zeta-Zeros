# Jensen-Window PF Rebalanced Power-Eight First-Summand Dominance

Date: 2026-07-13

Status: rigorous inverse-eighth-power first-summand dominance at
`lambda=-100`. This is not a proof of order seven, PF-infinity, RH,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py
```

## Rebalanced Split

The exact summand-ratio monotonicity and the rigorous
`epsilon(0)<0.0022` theorem are unchanged. Replace the former split by

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
q(a(k))=pi*k^(2/5).
```

The high-region exponential reserve remains large while the tilted
low-region probability gains the extra inverse power.

## Endpoint Gates

At `k=300`, outward-rounded Arb arithmetic gives

```text
high log comparison=[-43.81908591772582043884445769548078969524705689301795615052822235641345 +/- 4.52E-69]<0,
low log comparison=[-22.45423668980717202070624637558429877084711058474769615482656850706380 +/- 2.16E-69]<0,
high prefactor reserve=1.79791512439484206187063618727459739964492156552968028304521E-1>0,
minimum low power-derivative margin=2.33664016042713732356706005808783621003054381058540423959103E-1>0.
```

All fourteen endpoint and derivative gates are strictly positive.
The six displayed ratio-log derivatives prove that every comparison
strengthens for the complete half-line `k>=300`.

## Strengthened Tail

```text
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-8), k>=300
epsilon(0)*P_k(u<a(k))<k^(-8), k>=300
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^8 for every integer k>=300
0<=log(1+delta_k)<2/k^8, k>=300
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-8)+2*j^(-8)+(j+1)^(-8)), j>=301
```

The inverse-eighth-power wall error is the input needed by the fourth
stable logarithm in the order-seven complete-to-first transfer.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md
outputs/formal_core.md
```
