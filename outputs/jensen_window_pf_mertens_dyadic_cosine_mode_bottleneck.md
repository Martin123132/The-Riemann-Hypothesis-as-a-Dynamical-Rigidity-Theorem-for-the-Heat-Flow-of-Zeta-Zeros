# Jensen-Window PF Mertens Dyadic Cosine-Mode Bottleneck

Date: 2026-07-23

Status: exact dyadic cosine-mode reduction with automatically summable
high modes and one open RH-equivalent low-mode gate. This is not a proof
of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json
python work/rh_compute/scripts/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
```

## Dyadic Tail Blocks

Fix `0<alpha<1` and put

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.                    (MDCOS.1)
```

For a dyadic integer `K`, define the length-`K` tail vector and its
weighted block energy by

```text
y_m^(K):=r_(K+m-1),             1<=m<=K,
B_K:=sum_(N=K)^(2K-1)N^alpha|r_N|^2.                  (MDCOS.2)
```

The dyadic intervals partition the positive integers, and

```text
R_alpha=sum_(K dyadic)B_K,

K^alpha||y^(K)||_2^2
 <=B_K
 <=2^alpha K^alpha||y^(K)||_2^2.                      (MDCOS.3)
```

The adjacent-tail identity becomes

```text
y_m^(K)-y_(m+1)^(K)=q_(K+m),       1<=m<K.            (MDCOS.4)
```

Notice what (MDCOS.4) does cleanly: the future tail `r_(2K)` and the
coefficient `q_(2K)` add the same constant to every entry of `y^(K)`.
They therefore affect only the mean mode below.

## Exact Neumann Cosine Diagonalization

Let `D_K` be the first-difference map in (MDCOS.4). The orthonormal
Neumann cosine basis is

```text
phi_0(m):=K^(-1/2),

phi_r(m):=sqrt(2/K)
          cos(pi*r*(m-1/2)/K),       1<=r<K,           (MDCOS.5)

D_K^*D_K phi_r=nu_r phi_r,
nu_r:=4sin^2(pi*r/(2K)).
```

For

```text
c_(K,r):=<y^(K),phi_r>,
```

Parseval and gradient Parseval give the exact identities

```text
||y^(K)||_2^2
 =sum_(r=0)^(K-1)|c_(K,r)|^2,                          (MDCOS.6)

sum_(m=1)^(K-1)|q_(K+m)|^2
 =sum_(r=1)^(K-1)nu_r|c_(K,r)|^2.                     (MDCOS.7)
```

The mean coordinate retains the endpoint information:

```text
c_(K,0)
 =K^(-1/2)sum_(m=1)^K r_(K+m-1)

 =sqrt(K)r_(2K)
  +K^(-1/2)sum_(ell=1)^K ell*q_(K+ell).               (MDCOS.8)
```

Every nonconstant coordinate is a smooth sine-weighted Mobius sum:

```text
c_(K,r)
 =[sqrt(2/K)/(2sin(pi*r/(2K)))]
   sum_(m=1)^(K-1)
     mu(K+m)(K+m)^(-1-alpha)sin(pi*r*m/K),

1<=r<K.                                                (MDCOS.9)
```

This follows from
`c_(K,r)=nu_r^(-1)<D_K y^(K),D_K phi_r>`. In particular,
the endpoint tail cancels before any estimate is made.

## High Modes Are Free

For `1<=R<K`, monotonicity of `nu_r` and
`sin x>=2x/pi` on `[0,pi/2]` give

```text
nu_R>=4R^2/K^2.
```

Since `|mu(n)|<=1`,

```text
sum_(m=1)^(K-1)|q_(K+m)|^2<=K^(-1-2alpha).
```

Consequently,

```text
K^alpha sum_(r=R)^(K-1)|c_(K,r)|^2
 <=K^(1-alpha)/(4R^2).                                 (MDCOS.10)
```

Thus any cutoff

```text
R_K:=ceil(K^beta),       beta>(1-alpha)/2,             (MDCOS.11)
```

makes the omitted high-mode blocks dyadically summable. The universal
choice `R_K=ceil(sqrt(K))` gives

```text
K^alpha sum_(r>=R_K)|c_(K,r)|^2
 <=(1/4)K^(-alpha).
```

Combining this with (MDCOS.3) and (MDCOS.6) yields the exact reduction

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^alpha
  sum_(0<=r<ceil(sqrt(K)))|c_(K,r)|^2
 <infinity.                                             (MDCOS.12)
```

For any fixed cofinal sequence `alpha_j->0`, the previously checked
Mertens/tail criterion therefore gives

```text
RH
 iff
(MDCOS.12) holds for every alpha_j.                     (MDCOS.13)
```

This removes all but `O(sqrt(K))` smooth coordinates at scale `K`
without using cancellation. It does not estimate the remaining
coordinates.

## What The Lowest Modes Ask For

The first oscillatory coordinate is

```text
c_(K,1)
 =[sqrt(2/K)/(2sin(pi/(2K)))]
   sum_(m=1)^(K-1)q_(K+m)sin(pi*m/K).                  (MDCOS.14)
```

The sine window is positive on the whole block and the prefactor is
asymptotic to `sqrt(2K)/pi`. This is a one-hump anchored Mobius average,
not a rapidly oscillating mode. The zero mode is the block mean of the
tails.

Davenport's theorem states, for every `A>0`, uniformly in real `theta`,

```text
|sum_(n<=X)mu(n)e(n theta)|
 <<_A X(log X)^(-A).                                   (MDCOS.15)
```

Partial summation, with the sine in (MDCOS.9) split into the phases
`theta=+/-r/(2K)`, gives

```text
|c_(K,0)|
 <<_(alpha,A)K^(1/2-alpha)(log K)^(-A),

|c_(K,r)|
 <<_(alpha,A)
 K^(1/2-alpha)/(r(log K)^A),       1<=r<K.             (MDCOS.16)
```

Hence Davenport gives only

```text
K^alpha
sum_(0<=r<sqrt(K))|c_(K,r)|^2
 <<_(alpha,A)K^(1-alpha)(log K)^(-2A).                 (MDCOS.17)
```

The positive power `K^(1-alpha)` defeats every fixed logarithmic saving.
This does not prove that the actual low modes are large. It proves that
classical arbitrary-log linear-phase orthogonality does not close the
required dyadic series.

More generally, a hypothetical uniform estimate with `sigma<1+alpha`,

```text
|sum_(n<=X)mu(n)e(n theta)|<<X^sigma
```

would produce the block scale

```text
K^(2sigma-1-alpha).                                    (MDCOS.18)
```

It closes the dyadic sum for fixed `alpha` precisely at the level

```text
sigma<(1+alpha)/2.
```

In particular, the stronger unproved hypothesis

```text
sup_theta |sum_(n<=X)mu(n)e(n theta)|
 <<_epsilon X^(1/2+epsilon),

epsilon<alpha/2,                                       (MDCOS.19)
```

would make each low-mode block
`O(K^(-alpha+2epsilon))` and prove `R_alpha<infinity`.
Equation (MDCOS.19) is only a calibration of the missing exponent. It
contains the `theta=0` RH-scale Mertens estimate and is not claimed here
as a theorem or as a consequence of RH alone.

## Guards And Open Gate

Translation-averaged cancellation and almost-all short-interval theorems
do not automatically control the fixed dyadic mean and one-hump modes.
An explicit de-averaging theorem would be required. Likewise, this
reduction concerns `R_alpha` and the logarithmically summed
stable-prefix component `P`; it neither proves nor obstructs a cofinal
theorem for the full Burnol energy `Q`.

The primary linear-phase sources checked were:

- H. Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
- B. Green and T. Tao, *Quadratic Uniformity of the Mobius Function*,
  especially Example 3, Section 5, and the summation-by-parts appendix:
  https://doi.org/10.5802/aif.2401

A focused search found later exponential-sum refinements but no
unconditional result removing the `theta=0` power deficit needed by
(MDCOS.18) on a cofinal alpha sequence. This is a search report, not a
novelty, priority, or exhaustiveness claim.

The sharpened open target is

```text
for every alpha in one fixed cofinal sequence alpha_j->0,

sum_(K dyadic)K^alpha
  sum_(0<=r<ceil(sqrt(K)))|c_(K,r)|^2
 <infinity.                                             (MDCOS.20)
```

This is RH-equivalent by (MDCOS.13). The exact spectral reduction proves
the high modes harmless and identifies the low smooth modes, but supplies
no Mobius-specific gain. RH, PF-infinity, and `Lambda <= 0` remain open.
