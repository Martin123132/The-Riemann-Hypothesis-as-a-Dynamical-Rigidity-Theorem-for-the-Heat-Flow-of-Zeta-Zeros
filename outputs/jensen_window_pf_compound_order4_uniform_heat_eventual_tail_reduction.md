# Jensen-Window PF Compound Order-Four Uniform-Heat Eventual-Tail Reduction

Date: 2026-07-13

Status: exact uniform-asymptotic order-four reduction with one open
heat-tilt theorem. This is not a proof of unconditional order-four
invariance, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py
```

## Uniform Ratio Contract

For `I=[-100,0]`, `M=n+6`, and backward shifts `0<=j<=6`, it is
sufficient to prove uniformly in `lambda` that

```text
log(A_(M-j)(lambda)/A_M(lambda))=-G_1(lambda,M)*j-sum_(m=2)^7 G_m(lambda,M)*h(lambda,M)^(m-1)*j^m+r_(lambda,M,j)
sup_(lambda in I) h(lambda,M)->0 and sup_(lambda in I)|G_2(lambda,M)-1|->0
sup_(lambda in I)|G_m(lambda,M)|=O(1), m=3,...,7
max_(lambda in I,0<=j<=6)|r_(lambda,M,j)|=o(h(lambda,M)^6) uniformly
```

## Universal Determinant

The exact 24-permutation calculation is parameter-blind:

```text
[h^0,...,h^6] det K=[0,0,0,0,0,0,768*G_2^6].
det K=768*G_2(lambda,M)^6*h(lambda,M)^6+o(h(lambda,M)^6) uniformly.
```

Since `G_2->1` uniformly, this would prove

```text
there exists N such that H_(4,n)(lambda)>0 for every n>=N and lambda in [-100,0].
```

The finite-confinement theorem would then give

```text
the uniform tail plus lambda=-100 entry and the cooperative flow imply H_(4,n)(lambda)>0 for all n>=0 and lambda in [-100,0].
```

## Xi-Specific Target

Relative to the published lambda-zero coefficient sequence, set

```text
R_T^(1)(k)=A_k^(1)(-T)/A_k^(1)(0), 0<=T<=100.
```

A sufficient compact-heat estimate is

```text
Delta_k^m log R_T^(1)(k)=O(log(k)/k^m) uniformly for 0<=T<=100, m=2,...,7.
```

O(log(k)/k^m)=o(k^(-(m-1))) for m>=2, so the bounded heat tilt does not change the universal G_2 limit or the graded ratio scales.
The superpolynomial dominance theorem now closes the complete-kernel
correction and all seven required local log differences. The only remaining
asymptotic input is the displayed first-summand heat-tilt estimate.

```text
outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md
outputs/jensen_window_pf_lambda0_first_summand_dominance_transfer.md
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.md
outputs/signed_hankel_jensen_dependency_graph.md
```
