# Graded-Kernel Vandermonde Lemma At Every Fixed Order

Date: 2026-07-13

Status: exact all-fixed-order leading-determinant theorem and
compact-heat uniform eventual signed-tail theorem. This is not a proof of
all-shift sign regularity, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json
python work/rh_compute/scripts/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py
```

## Kernel And Gauge

Fix an integer `m>=1`, put `D=binom(m,2)`, `L=2(m-1)`,
`z_i=L-i`, and `y_j=j`. The normalized determinant kernel is

```text
K_(i,j)=exp(-sum_(r>=2)G_r*h^(r-1)*(L-i-j)^r)
```

The binomial identity

```text
(z-y)^r=z^r+(-y)^r
        +sum_(a,b>=1,a+b=r) binom(r,a)z^a(-y)^b
```

splits off row and column factors with constant term one. They cannot
alter the first nonzero determinant coefficient. The remaining mixed
kernel is

```text
M_h(z,y)=exp(sum_(a,b>=1)c_(a,b)*h^(a+b-1)*z^a*y^b)
c_(1,1)=2*G_2
```

## Valuation Lemma

Write `M_h(z,y)=sum c_(p,q)(h)z^p y^q`. A product of `ell` mixed
factors contributing bidegree `(p,q)` has exact `h`-degree

```text
p+q-ell>=max(p,q), because ell<=min(p,q).
```

Therefore

```text
ord_h c_(p,q)>=max(p,q),
[h^k]c_(k,k)=(2*G_2)^k/k!.
```

Equality in the diagonal formula requires `ell=k` and hence `k`
copies of the sole `(a,b)=(1,1)` quadratic interaction. No `G_r`
with `r>=3` can enter.

## Cauchy-Binet Floor

Expand the mixed kernel as a monomial matrix product. For increasing
nonnegative exponent tuples `P,Q` and a coefficient permutation `pi`,

```text
sum_i max(p_i,q_(pi(i)))
 >=max(sum_i p_i,sum_i q_i)>=D.
```

Equality forces `P=Q=(0,1,...,m-1)`. In that block

```text
sum_i max(i,pi(i))=D+(1/2)*sum_i|i-pi(i)|,
```

so only the identity permutation remains at degree `D`. The checker
stress-tests this identity on all `46233`
permutations through order eight.

## Universal First Term

The two monomial alternants are Vandermonde determinants. Since
`z_i` is decreasing and `y_j` increasing, their product leaves

```text
epsilon_m*2^D*prod_(j=1)^(m-1)j!*G_2^D
```

as the first raw coefficient. Thus `Q_(m,n)=epsilon_m H_(m,n)` has
the positive first coefficient

```text
2^D*prod_(j=1)^(m-1)j!*G_2^D>0
```

## Every Fixed Heat-Difference Order

O'Sullivan's Theorem 5.2 supplies the saddle expansion at arbitrary
fixed truncation order for suitable multipliers. Cauchy estimates in
the published complex disk make the compact Gaussian-log family
uniformly suitable at every fixed order. For `w=W(2k/pi)`, put

```text
F_0(w)=w^2,
F_(s+1)(w)=-s*F_s(w)+w/(1+w)*F_s'(w).
```

Induction using `w'=w/(k(1+w))` gives

```text
d^s(w^2)/dk^s=k^(-s)F_s(w), F_s(w)=O_s(w).
```

For one fixed difference order `s`, choose the published saddle
truncation `R>s`. The explicit correction terms gain `s` inverse
powers under the integral finite-difference formula. The difference
of the remainder is bounded directly, without differentiating it:

```text
2^s O_R(w^(3R)/k^R)=o(w/k^s).
```

Consequently

```text
Delta^s log R_T^(1)(k)=O_s(log(k)/k^s)
for every fixed s>=2, uniformly for 0<=T<=100.
```

The constants and the truncation may depend on `s`; this is not a
uniform-in-order estimate.

The all-order Xi ratio theorem, compact heat-tilt theorem, and
superpolynomial higher-theta suppression now imply

```text
for every fixed m, exists N_m such that Q_(m,n)(lambda)>0 for n>=N_m and -100<=lambda<=0
```

This is an `all fixed m` statement: the finite threshold may depend on
`m`. It neither reverses the quantifiers nor fills any finite prefix.

## Order Seven

At the new frontier `D=21` and `epsilon_7=-1`, so

```text
[h^0,...,h^21]det K=[0,...,0,-52183852646400*G_2^21],
52183852646400*G2**21*h**21 for Q_7=-H_7.
```

This closes the eventual order-seven tail. Endpoint entry at
`lambda=-100`, its analytic finite-prefix splice, and the all-order
sign-regularity bridge remain separate obligations.

```text
outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md
outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md
outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
