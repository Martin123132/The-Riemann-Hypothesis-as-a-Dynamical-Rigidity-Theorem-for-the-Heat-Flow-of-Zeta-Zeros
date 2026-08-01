# Jensen-Window PF OU/Mertens Mean-Square Reduction

Date: 2026-07-23

Status: exact OU-tail/Mertens/dyadic/anchored-correlation reduction with one
open RH-equivalent mean-square gate. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_ou_mertens_mean_square_reduction.json
python work/rh_compute/scripts/jensen_window_pf_ou_mertens_mean_square_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_ou_mertens_mean_square_reduction.py
```

## Tail And Mertens Coordinates

Fix `0<alpha<1`. The prime number theorem gives

```text
sum_(n>=1)mu(n)/n=0,
M(N)/N->0,
```

so define

```text
r_k:=sum_(n>=k)mu(n)/n,
M(k):=sum_(n<=k)mu(n),

B_alpha:=sum_(k>=1)k^(-alpha)|r_k|^2,
M_alpha:=sum_(k>=1)M(k)^2/k^(2+alpha).                       (OMMSR.1)
```

The difference and two Abel identities are

```text
r_k-r_(k+1)=mu(k)/k,                                         (OMMSR.2)

M(k)=sum_(n=1)^k r_n-k*r_(k+1),                              (OMMSR.3)

r_(k+1)
 =-M(k)/(k+1)+sum_(n>=k+1)M(n)/(n(n+1)).                    (OMMSR.4)
```

Let

```text
Hf(k)=(1/k)sum_(n<=k)f(n),
Cf(k)=sum_(n>=k)f(n)/n.
```

A direct Schur test on `l2(k^(-alpha))` gives the nonsharp explicit
constants

```text
||H||<=h_alpha:=sqrt(2(3+alpha))/(1+alpha),
||C||<=c_alpha:=sqrt(2(3-alpha))/(1-alpha).                  (OMMSR.5)
```

Equations (OMMSR.3)-(OMMSR.5), including the one-step shift, imply

```text
M_alpha^(1/2)
 <=(h_alpha+2^(alpha/2))*B_alpha^(1/2),

B_alpha^(1/2)
 <=(1+c_alpha)*M_alpha^(1/2).                                (OMMSR.6)
```

Thus the reciprocal-Mobius tail norm is finite exactly when the weighted
Mertens mean square is finite.

## OU Energy Is The Same Barrier

For

```text
w_(alpha,k):=k^(1-alpha)-(k-1)^(1-alpha),
```

concavity gives

```text
(1-alpha)k^(-alpha)
 <=w_(alpha,k)
 <=k^(-alpha).                                               (OMMSR.7)
```

Therefore the limiting OU energy from Lemma 11.22Q satisfies

```text
(1-alpha)B_alpha
 <=sum_(k>=1)w_(alpha,k)|r_k|^2
 <=B_alpha,                                                  (OMMSR.8)
```

and is norm-equivalent to `M_alpha`. The OU extraction supplies useful
positive structure, but it does not bypass the classical Mertens
mean-square barrier.

The continuous version is exact:

```text
I_alpha
 :=integral_1^infinity M(x)^2*x^(-2-alpha)dx

 =(1/(1+alpha))*sum_(k>=1)M(k)^2
   [k^(-1-alpha)-(k+1)^(-1-alpha)],                          (OMMSR.9)

2^(-2-alpha)M_alpha<=I_alpha<=M_alpha.                       (OMMSR.10)
```

The classical Mellin identity

```text
1/zeta(s)
 =s*integral_1^infinity M(x)x^(-s-1)dx                      (OMMSR.11)
```

starts in `Re(s)>1`. If `I_alpha<infinity`, Cauchy-Schwarz continues the
Mellin transform to `Re(s)>(1+alpha)/2`, and Mellin-Plancherel gives, with
`beta=(1+alpha)/2`,

```text
I_alpha
 =(1/(2*pi))*integral_R
   dt/[(beta^2+t^2)|zeta(beta+it)|^2].                       (OMMSR.12)
```

Finiteness excludes a reciprocal-zeta pole on that boundary. This is a
reciprocal-zeta Hardy energy, not an unconditional estimate for it.

## Exact Dyadic Criterion

Define

```text
D_j
 :=2^(-2j)sum_(2^j<=k<2^(j+1))M(k)^2.                       (OMMSR.13)
```

Block comparison gives

```text
2^(-2-alpha)sum_(j>=0)2^(-alpha*j)D_j
 <=M_alpha
 <=sum_(j>=0)2^(-alpha*j)D_j.                               (OMMSR.14)
```

The root test therefore gives

```text
M_alpha<infinity for every alpha>0
 iff
limsup_(j->infinity)D_j^(1/j)<=1.                            (OMMSR.15)
```

Composed with the cofinal OU criterion:

```text
RH
 iff
for every epsilon>0,
sum_(K<=k<2K)M(k)^2
 =O_epsilon(K^(2+epsilon))
on dyadic K.                                                 (OMMSR.16)
```

This is an exact RH-equivalent average-cancellation target. Equation
(OMMSR.16) has not been proved here.

## Anchored Correlation Form

Put

```text
S(X):=sum_(k<=X)M(k)^2,
C_h(Y):=sum_(m<=Y)mu(m)mu(m+h).
```

Expanding every prefix square and counting the containing prefixes gives

```text
S(X)
 =sum_(n<=X)mu(n)^2(X-n+1)

 +2sum_(h=1)^(X-1)sum_(m=1)^(X-h)
   mu(m)mu(m+h)(X-m-h+1)                                    (OMMSR.17)

 =sum_(n<=X)mu(n)^2(X-n+1)
  +2sum_(h=1)^(X-1)sum_(Y=1)^(X-h)C_h(Y).                   (OMMSR.18)
```

The squarefree diagonal is `O(X^2)`. Hence (OMMSR.16) is equivalently the
signed, origin-anchored integrated-correlation estimate

```text
sum_(h=1)^(X-1)sum_(Y=1)^(X-h)C_h(Y)
 =O_epsilon(X^(2+epsilon)).                                  (OMMSR.19)
```

The signs and the cumulative `Y` summation are essential. Taking absolute
values creates a much stronger target.

## Two Nonpromotion Guards

First, a generic weighted Hardy inequality applied to
`M(k)=sum_(n<=k)mu(n)` gives only

```text
sum_(k<=N)M(k)^2/k^(2+alpha)
 <=C_alpha*sum_(n<=N)mu(n)^2*n^(-alpha)
 =O_alpha(N^(1-alpha)).                                      (OMMSR.20)
```

This reproduces the divergent absolute Gram scale. Hardy/Copson machinery
transfers coordinates but supplies no Mobius cancellation.

Second, Matomaki, Radziwill, and Tao prove an averaged Chowla theorem of
the form

```text
sum_(h_1,...,h_k<=H)
 |sum_(n<=X)mu(n+h_1)...mu(n+h_k)|
 =o(H^k*X),                                                  (OMMSR.21)
```

with the Mobius extension stated in the paper. This averages all base
shifts. By contrast, `C_h(Y)` fixes the exceptional origin slice
`h_1=0` and then integrates every terminal length `Y`. Even a hypothetical
anchored terminal estimate

```text
sum_(h<=Y)|C_h(Y)|=o(Y^2)
```

used naively for every `Y<=X` gives cubic scale after summing over `Y`,
not the `X^(2+epsilon)` target. Almost-all short-interval cancellation
likewise leaves the single prefix `M(X)` uncontrolled. Therefore averaged
Chowla, logarithmically averaged Chowla, and almost-all short-interval
cancellation cannot be promoted into (OMMSR.19) without a new anchored
transfer theorem and an additional full power of cancellation.

## Open Handoff

The exact live target is any one of the equivalent statements

```text
D_j=2^(o(j)),

sum_(K<=k<2K)M(k)^2=O_epsilon(K^(2+epsilon)),

sum_(h<X)sum_(Y<=X-h)C_h(Y)=O_epsilon(X^(2+epsilon)),
```

proved without assuming a zero-free half-plane, RH, or the desired OU
bound. The reduction does not prove this target, the full Burnol norm,
RH, PF-infinity, or `Lambda <= 0`.

## Source Boundary

- [Lee and Leong](https://arxiv.org/abs/2208.06141) record modern explicit
  Mertens and reciprocal-zeta bounds and the standard Mellin connection.
  The exact tail/Mertens and dyadic reductions above are derived here.
- [Das and Manna](https://arxiv.org/abs/2508.00388) review and sharpen
  discrete Hardy/Copson inequalities. Only generic operator bounds are
  used here, and they do not supply arithmetic cancellation.
- [Matomaki, Radziwill, and Tao](https://arxiv.org/abs/1503.05121) prove
  averaged Chowla and explicitly note the Mobius extension. Their theorem
  does not state the anchored cumulative estimate (OMMSR.19).
- [Matomaki and Radziwill](https://annals.math.princeton.edu/2016/183-3/p06)
  prove cancellation for Mobius in almost all short intervals. The
  exceptional origin-anchored prefix remains outside that conclusion.
