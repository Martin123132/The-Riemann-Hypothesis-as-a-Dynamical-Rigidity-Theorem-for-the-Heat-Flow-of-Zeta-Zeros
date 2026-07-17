# Jensen-Window PF Rebalanced Power-Eleven First-Summand Dominance

Date: 2026-07-16

Status: rigorous inverse-eleventh-power first-summand dominance at
`lambda=-100`. This is not a proof of order ten, PF-infinity, RH,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.py
```

## Reused Split

```text
a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.
```

At `k=300`, all fourteen Arb endpoint and derivative gates remain
strictly positive, and the six exact ratio-log derivatives cover the
complete half-line.

```text
epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-11), k>=300
epsilon(0)*P_k(u<a(k))<k^(-11), k>=300
0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^11 for every integer k>=300
0<=log(1+delta_k)<2/k^11, k>=300
|B_j-B_j^(1)|<=a_j=2*((j-1)^(-11)+2*j^(-11)+(j+1)^(-11)), j>=301
```

The next power requires a slightly later endpoint gate and is certified
separately rather than inferred from this artifact.
