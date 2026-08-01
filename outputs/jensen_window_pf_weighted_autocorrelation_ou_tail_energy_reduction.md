# Jensen-Window PF Weighted Autocorrelation OU Tail-Energy Reduction

Date: 2026-07-23

Status: exact weighted far-ratio/Ornstein-Uhlenbeck reduction and cofinal
RH-equivalent reciprocal-Mobius tail energy, with one open comparison gate.
This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.json
python work/rh_compute/scripts/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py
```

## Quantitative Period Averaging

Fix `0<alpha<1` and put

```text
a=(1-alpha)/2,
beta=(1+alpha)/2.
```

The amplitude

```text
f_alpha(x)={x}x^(alpha-2)
```

is integrable at both endpoints. Its half-mean is

```text
L_alpha
 =(1/2)integral_(0,infinity)f_alpha(x)dx
 =-zeta(1-alpha)/(2(1-alpha)) > 0.                            (WAOTER.1)
```

Let `Psi(y)=B_2({y})/2`. Then `Psi` is continuous and periodic,
`Psi'={y}-1/2` almost everywhere, and `||Psi||_infinity=1/12`.
Splitting the centered autocorrelation integral at `h=1/lambda` and
integrating by parts on `[h,infinity)` gives, for `lambda>=1`,

```text
|A_alpha(lambda)-L_alpha|
 <=K_alpha*lambda^(-alpha),                                  (WAOTER.2)

K_alpha=1/(2alpha)+zeta(2-alpha)/6.
```

Indeed, the first interval costs at most
`lambda^(-alpha)/(2alpha)`. On `[h,infinity)`, the total variation of
`f_alpha` is at most

```text
h^(alpha-1)+2(zeta(2-alpha)-1).
```

The boundary term contributes one further `h^(alpha-1)`, and
`lambda^(-1)<=lambda^(-alpha)`. Thus (WAOTER.2) is fully explicit.

For the stationary kernel from the preceding Gram bridge,

```text
C_alpha(u)=exp(-a*u)A_alpha(exp(u)), u>=0,
```

evenness and `a+alpha=beta` now give

```text
C_alpha(u)
 =L_alpha*exp(-a|u|)+R_alpha(u),                              (WAOTER.3)

|R_alpha(u)|<=K_alpha*exp(-beta|u|).                          (WAOTER.4)
```

The leading term is the stationary Ornstein-Uhlenbeck covariance kernel.

## The Remainder Is Not Positive

The Fourier densities are

```text
C_alpha:
  |zeta(a+it)|^2/(a^2+t^2),

L_alpha*exp(-a|u|):
  2aL_alpha/(a^2+t^2),

R_alpha:
  [|zeta(a+it)|^2-2aL_alpha]/(a^2+t^2).                      (WAOTER.5)
```

At `alpha=1/2`, `a=1/4`, and `t=0`,

```text
|zeta(1/4)|^2-2*(1/4)*L_(1/2)
 =-0.0687554899394677... < 0.                                (WAOTER.6)
```

Therefore the natural remainder is not positive definite. The asymptotic
split cannot be promoted to a lower or upper comparison between its two
Gram forms.

## Exact OU Tail Energy

Set

```text
c_d=mu(d)d^(-beta),

G_(alpha,N)
 =sum_(d,e<=N)c_d*c_e*exp(-a|log(d/e)|).                     (WAOTER.7)
```

The elementary factorization

```text
exp(-a|u-v|)
 =2a*integral_(-infinity,min(u,v))
   exp(-a(u-x))exp(-a(v-x))dx                               (WAOTER.8)
```

and `beta+a=1` give

```text
G_(alpha,N)
 =(1-alpha)integral_0^N y^(-alpha)
   |sum_(y<=d<=N)mu(d)/d|^2dy                               (WAOTER.9)

 =sum_(k=1)^N [k^(1-alpha)-(k-1)^(1-alpha)]
   |sum_(d=k)^N mu(d)/d|^2.                                 (WAOTER.10)
```

This is a positive, exact reciprocal-Mobius tail energy.

## Cofinal RH Equivalence

The prime number theorem gives

```text
sum_(d>=1)mu(d)/d=0.
```

Hence, for each fixed `k`,

```text
sum_(d=k)^N mu(d)/d
 ->r_k:=sum_(d=k)^infinity mu(d)/d.
```

If `sup_N G_(alpha,N)<infinity`, Fatou and concavity give

```text
sum_(k>=1)k^(-alpha)|r_k|^2<infinity.                        (WAOTER.11)
```

For

```text
beta=(1+alpha)/2,

F_alpha(s)
 :=sum_(k>=2)r_k[k^(1-s)-(k-1)^(1-s)],                       (WAOTER.12)
```

Cauchy-Schwarz makes the series locally uniformly convergent in
`Re(s)>beta`. Summation by parts in `Re(s)>1` gives
`F_alpha(s)=1/zeta(s)`, so (WAOTER.12) analytically continues the
reciprocal zeta function and proves

```text
zeta(s)!=0 for Re(s)>(1+alpha)/2.                            (WAOTER.13)
```

Conversely, RH gives `M(x)=O_epsilon(x^(1/2+epsilon))`. Partial summation
then yields, uniformly in `N>=k`,

```text
|sum_(d=k)^N mu(d)/d|=O_epsilon(k^(-1/2+epsilon)).
```

Choosing `2epsilon<alpha` makes (WAOTER.10) uniformly summable. Therefore

```text
RH
 iff
sup_N G_(1/(j+1),N)<infinity for every integer j>=1.         (WAOTER.14)
```

This is an exact auxiliary RH-equivalent criterion. It does not establish
its antecedent.

## Why The Full Gram Is Still Open

The full Burnol Gram form is only

```text
Gamma_(alpha,N)
 =L_alpha*G_(alpha,N)+E_(alpha,N),                            (WAOTER.15)
```

where `E` uses `R_alpha(log(d/e))`. Equation (WAOTER.6) proves that no
positive-semidefinite ordering follows. Even (WAOTER.4), used absolutely,
costs `O(N^(1-alpha))`.

The Montgomery-Vaughan weighted Hilbert inequality also does not close the
gap. With frequencies `lambda_n=log n`, their local spacing is
`delta_n=log(1+1/n)` and its generic right side costs

```text
sum_(n<=N)|c_n|^2/delta_n
 asymptotic to a constant times sum_(n<=N)mu(n)^2n^(-alpha)
 asymptotic to a constant times N^(1-alpha).                 (WAOTER.16)
```

Moreover, that theorem controls the skew kernel
`1/(lambda_m-lambda_n)`, not the bounded ratio kernel here. The
multiplicative Hilbert matrix instead uses a product kernel. Both are
theorem-search anchors, not ready-to-apply solutions.

The open target is a Mobius-specific comparison, cancellation, or bilinear
estimate that controls the sign-indefinite remainder together with the
RH-equivalent OU tail energy on one cofinal sequence. The reduction itself
does not prove the full Burnol bound, RH, PF-infinity, or `Lambda <= 0`.

## Source Boundary

- [Báez-Duarte, Balazard, Landreau, and Saias](https://arxiv.org/abs/math/0306251)
  derive the `alpha=0` fractional-part autocorrelation, scaling, and Mellin
  formula. The weighted `0<alpha<1` period mean and (WAOTER.2) are derived
  here.
- [Yangjit](https://arxiv.org/abs/2203.14950) states the
  Montgomery-Vaughan weighted Hilbert spacing inequality used in the
  theorem-fit guard.
- [Brevig, Perfekt, Seip, Siskakis, and Vukotic](https://arxiv.org/abs/1411.7294)
  study the product-kernel multiplicative Hilbert matrix. It is not the
  present ratio-stationary matrix.
