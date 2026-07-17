# Jensen-Window PF Compound Order-Four Uniform-Heat Forward Invariance

Date: 2026-07-13

Status: all-shift contiguous order-four forward invariance through lambda
zero. This is not a proof of arbitrary-column order four, PF-infinity, RH,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py
```

## Uniform Ratio Expansion

The published lambda-zero Xi ratio theorem supplies the base graded
expansion. Two new uniform theorems control the heat deformation:

```text
Delta_k^m log(A_k^(1)(-T)/A_k^(1)(0))=O(log(k)/k^m), 2<=m<=7, uniformly for 0<=T<=100
all fixed local log differences of log(A_k(-T)/A_k^(1)(-T)) are O_p,m(k^-p) uniformly for every p>0
```

Therefore

```text
the complete A_k(lambda) family satisfies the uniform degree-seven ratio contract on -100<=lambda<=0 with the same G_2 limit 1.
```

On the seven determinant nodes, exact Newton interpolation is triangular:
the coefficient of `j^r` uses only differences of order at least `r`.
Thus the heat correction at degree `r` is `O(log(M)/M^r)`, and

```text
c_r(T,M)=O(log(M)/M^r) and h(M)~1/(2M) imply c_r(T,M)/h(M)^(r-1)=O(log(M)/M)=o(1), r=2,...,6.
```

## Uniform Positive Tail

The exact determinant cancellation is universal and gives

```text
H_(4,n)(lambda)=positive_scale(lambda,n)*(768*G_2(lambda,n)^6*h(lambda,n)^6+o(h(lambda,n)^6)) uniformly.
```

Since `G_2->1` uniformly, there is one finite, non-effective `N` with

```text
there exists N such that H_(4,n)(lambda)>0 for every n>=N and -100<=lambda<=0.
```

## Backward Cooperative Propagation

At `lambda=-100`, every shift is already positive. Throughout the interval
the order-three theorem gives

```text
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0.
```

Starting at the uniform positive tail boundary, variation of constants
propagates strict positivity backward through the finite remaining indices:

```text
Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds).
```

Hence

```text
H_(4,n)(lambda)>0 for every integer n>=0 and every lambda in [-100,0]
H_(4,n)(0)>0 for every integer n>=0.
```

This closes the previous non-effective lambda-zero finite splice without
making its threshold explicit. The next layer is not RH: noncontiguous
order-four minors and every compound order at least five remain open.

```text
outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md
outputs/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.md
outputs/signed_hankel_jensen_dependency_graph.md
```
