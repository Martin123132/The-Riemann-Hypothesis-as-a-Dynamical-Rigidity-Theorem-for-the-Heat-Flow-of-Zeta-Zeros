# Jensen-Window PF Compound Order-Four Lambda-Zero Eventual Positivity

Date: 2026-07-13

Status: exact eventual lambda-zero order-four positivity theorem plus a
rigorous `0<=n<=500` prefix, with the effective splice still open. This is
not a proof of all-shift order four at lambda zero, PF-infinity, RH, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.py
```

## Exact Xi Normalization

The standard identities

```text
H_0(x)=xi((1+i*x)/2)/8
xi(1/2+z)=sum_(k>=0) gamma(k)*z^(2k)/k!
A_k(0)=2*k!/(2k)!*integral_0^infinity Phi(u)*u^(2k)du
```

give

```text
A_k(0)=gamma(k)/4^(k+1)
H_(4,n)[A(0)]=4^(-4*n-16)*H_(4,n)[gamma]
```

so the lambda-zero `A_k` and standard `gamma(k)` determinants have the
same sign.

## Published Ratio Input

Theorem 2.1 of Griffin, Ono, Rolen, Thorner, Tripp, and Wagner,
[Jensen Polynomials for the Riemann Xi Function](https://arxiv.org/abs/1910.01227),
gives, for `1<=j<M`,

```text
log(gamma(M-j)/gamma(M))=-sum_(m>=1) G_m(M)*Delta(M)^(2*m-2)*j^m
Delta(M)~1/sqrt(2*M) and G_m(M)->2^(m-1)/(m*(m-1))
|G_m(M)|<<_C(2*C)^m uniformly in m and M
```

This use is directly through the Xi coefficient asymptotics. No claim that
Jensen hyperbolicity by itself implies the signed Hankel sign is used.

## Universal Determinant Term

Set `M=n+6`, `h=Delta(M)^2`, and remove the affine `G_1` factor from
rows and columns. The remaining matrix has entries

```text
K_(i,j)(h)=exp(-sum_(m>=2) G_m*h^(m-1)*(6-i-j)^m)
```

Exact 24-permutation truncated-series algebra gives

```text
[h^0,...,h^6] det(K) =
['0', '0', '0', '0', '0', '0', '768*G2**6']
768*G2^6*h^6
```

All `G_3,...,G_7` contributions cancel from the first nonzero term. The
published tail bound controls `m>=8`, while `G_2(M)->1`. Consequently

```text
H_(4,n)[gamma]=gamma(M)^4*exp(-12*G_1(M))*(768*G_2(M)^6*Delta(M)^12+o(Delta(M)^12))
there exists N_H4 such that H_(4,n)(0)>0 for every n>=N_H4
```

This is an actual eventual positivity theorem for the zeta/Xi sequence,
not merely numerical evidence.

## Rigorous Prefix

A direct 24,576-bit Arb expansion of `xi(1/2+z)` through `A_506(0)`
give

```text
H_(4,n)(0)>0 for every 0<=n<=500
smallest recorded stable margin occurs at n=500
F_500=[5.481268356518127684426961803646365408360017160550209264948367373383654753206792239922810852195036203e-20 +/- 2.95e-120]
F_500/G_501^2=[0.003716298961440233294098765075005020047401588139987551949881424960342387775139012190166540704409788816 +/- 1.73e-103]
```

All 501 raw `H_4` balls and all 501 stable `F_n` balls are strictly
positive. The full interval rows and source hashes are stored in the JSON
artifact.

## Remaining Splice

The current result is

```text
finite theorem:   H_(4,n)(0)>0 for 0<=n<=500
eventual theorem: H_(4,n)(0)>0 for n>=N_H4 for some finite N_H4
open splice:      make N_H4 explicit and cover 501<=n<N_H4
```

The next analytic job is to extract explicit constants from the convergent
ratio expansion, bound the normalized determinant remainder against
`768*G_2^6`, and then run Arb only over the resulting finite gap.

```text
outputs/jensen_window_pf_compound_order4_forward_flow_reduction.md
outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md
outputs/arb_xi_lambda0_order4_prefix_certificate.md
outputs/formal_core.md
```
