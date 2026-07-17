# Jensen-Window PF Uniform Superpolynomial First-Summand Dominance

Date: 2026-07-13

Status: uniform compact-heat superpolynomial first-summand dominance
theorem. This is not a proof of uniform order-four positivity, PF-infinity,
RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json
python work/rh_compute/scripts/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py
python work/rh_compute/scripts/check_jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py
```

## Low And High Regions

The proved lambda=-100 split uses `a(k)=log(k)/8`. Its low-region
exponent satisfies

```text
B_low(k)=-8*k/(5*(log(k)+4/5))+(5/2)*log(k)+1+pi*(exp(2/5)-1)*sqrt(k)
lim_(k->infinity) (B_low(k)+p*log(k))*log(k)/k=-8/5.
```

The high-region estimate satisfies

```text
epsilon(a(k))<=17*exp(-3*pi*sqrt(k))
lim_(k->infinity) (log(17)-3*pi*sqrt(k)+p*log(k))/sqrt(k)=-3*pi.
```

Both limits are negative for every fixed `p>0`, so the two regions give

```text
for every p>0 there exists K_p such that 0<=delta_k(100)<=2*k^(-p), k>=K_p.
```

## Uniform Heat Transfer

The covariance monotonicity theorem gives `delta_k(T)<=delta_k(100)`.
Consequently

```text
for every p>0 there exists K_p such that sup_(0<=T<=100) delta_k(T)<=2*k^(-p), k>=K_p.
```

For `e_k(T)=log(1+delta_k(T))`, the elementary stencil bound

```text
|nabla^m e_k(T)|<=2^m*max_(0<=j<=m)e_(k+j)(T)
```

proves every fixed local log difference is uniformly superpolynomial.
Thus higher theta summands are negligible to every fixed order in the
uniform order-four ratio expansion. The remaining asymptotic problem is
only the first-summand heat tilt.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_lambda0_first_summand_dominance_transfer.md
outputs/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.md
outputs/signed_hankel_jensen_dependency_graph.md
```
