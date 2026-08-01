# Jensen-Window PF Mertens Spectral Anchor/High-Mode Reduction

Date: 2026-07-23

Status: exact spectral anchor separation and summable high-mode
reduction with one open low-mode gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.json
python work/rh_compute/scripts/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py
```

## Exact Current/Anchor Split

Fix `0<alpha<1`. Corollary 11.22Z.1 proves

```text
R_alpha<infinity
 iff
S_alpha
 :=sum_(K dyadic)K^(-alpha)
   sum_(r=1)^K |X_(K,r)|^2/(2r-1)^2
 <infinity.                                           (SAHR.1)
```

Write

```text
X_(K,r)=Y_(K,r)+beta_K*h_(K,r),

Y_(K,r)
 :=sum_(i=1)^(K-1)mu(K+i)phi_(K,r)(i),

h_(K,r):=phi_(K,r)(K).                               (SAHR.2)
```

Since

```text
K*theta_(K,r)
 =(2r-1)pi/2-theta_(K,r)/2,
```

the boundary amplitude is exactly

```text
h_(K,r)
 =2(-1)^(r-1)/sqrt(2K+1)
  *cos(theta_(K,r)/2),

|h_(K,r)|^2<=4/(2K+1)<=2/K.                         (SAHR.3)
```

The current block has full Parseval control:

```text
sum_(r=1)^K|Y_(K,r)|^2
 =sum_(i=1)^(K-1)mu(K+i)^2
 <=K.                                                (SAHR.4)
```

## High-Mode Inequality

For an integer `1<=R<=K`, the odd tail satisfies

```text
sum_(r>R)1/(2r-1)^2<=1/(2R+1).                      (SAHR.5)
```

Therefore

```text
sum_(r>R)|Y_(K,r)|^2/(2r-1)^2
 <=K/(2R+1)^2,

sum_(r>R)|beta_K*h_(K,r)|^2/(2r-1)^2
 <=2|beta_K|^2/[K(2R+1)].                           (SAHR.6)
```

Using `|a+b|^2<=2|a|^2+2|b|^2` gives the exact working bound

```text
sum_(r>R)|X_(K,r)|^2/(2r-1)^2
 <=2K/(2R+1)^2
   +4|beta_K|^2/[K(2R+1)].                          (SAHR.7)
```

The two terms have different origins. The current oscillatory block
decays like `R^(-2)`, while the affine future anchor decays only like
`R^(-1)`.

## Unconditional Summable Tail

Every coefficient in `beta_K` lies in `[0,1]`, so

```text
|beta_K|<=2K.                                        (SAHR.8)
```

Choose

```text
R_(alpha,K):=ceil(K^(1-alpha/2)).                    (SAHR.9)
```

Multiplying (SAHR.7) by `K^(-alpha)` gives

```text
K^(-alpha)
sum_(r>R_(alpha,K))|X_(K,r)|^2/(2r-1)^2
 <=2K^(-1)+8K^(-alpha/2).                           (SAHR.10)
```

Both terms are dyadically summable. Hence the high modes are closed
unconditionally, and

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-alpha)
 sum_(1<=r<=R_(alpha,K))
 |X_(K,r)|^2/(2r-1)^2
 <infinity.                                          (SAHR.11)
```

For every fixed positive `alpha`, only
`O(K^(1-alpha/2))` modes remain. This statement is pointwise in
`alpha`; it gives neither an `alpha=0` reduction nor uniformity as
`alpha->0`.

## Conditional Stronger Cutoff

Suppose, conditionally, that

```text
|beta_K|=O_epsilon(K^(1/2+epsilon)),
                  2epsilon<alpha.                   (SAHR.12)
```

For any `eta>0` satisfying `(1-alpha)/2+eta<1`, choose

```text
R=ceil(K^((1-alpha)/2+eta)).
```

Then the current contribution in (SAHR.7), after multiplication by
`K^(-alpha)`, is `O(K^(-2eta))`; the boundary contribution is also
dyadically summable because `2epsilon<alpha`. Under (SAHR.12), only

```text
O(K^((1-alpha)/2+eta))                               (SAHR.13)
```

low modes remain. This stronger cutoff is conditional.

## Guards

Davenport and partial summation give
`beta_K<<_A K*log^(-A)K`, but logarithmic savings do not change the
power threshold in (SAHR.7). A generic bounded-sequence argument cannot
do better: if the synthetic future coefficients equal one on
`[2K,3K]`, then every corresponding weight is at least `2/9`, so
`beta_K>=2K/9`.

Thus the unconditional low-mode count cannot be improved by coefficient
envelopes or generic orthogonality. A better cutoff requires genuine
Mobius power cancellation in the anchor, while closing the remaining
low modes requires the signed additive-twist estimate from Corollary
11.22Z.1.

## Open Gate

For every member of one fixed cofinal sequence `alpha_j->0`, prove the
low-mode series in (SAHR.11) without assuming RH. The high modes are now
an exact closed component. Low-mode square-root cancellation, anchor
power cancellation, the full Burnol bound, RH, PF-infinity, and
`Lambda<=0` all remain open.

Primary source:

- Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
