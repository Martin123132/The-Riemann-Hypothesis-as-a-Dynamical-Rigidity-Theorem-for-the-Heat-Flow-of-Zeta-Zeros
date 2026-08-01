# Jensen-Window PF Mertens Ordinary-Vector Vaughan Handoff

Date: 2026-07-23

Status: exact ordinary-Mobius feature Gram and pre-square vector Vaughan
reduction with one open signed gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py
```

## Input From The Affine-Tent Lemma

Fix `0<alpha<1`, put `s=1+alpha`, and let `K` run over dyadic
integers. Formal Lemma 11.22Y proves

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)H_K<infinity,               (MOVH.1)
```

where

```text
beta_K
 :=mu(2K)
   +sum_(2K<n<4K)
     ((2K)/n)^s*(4K-n)/(2K)*mu(n),

H_K
 :=sum_(t=1)^K
   |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2.                (MOVH.2)
```

Write

```text
b_K(n):=((2K)/n)^s*(4K-n)/(2K),     2K<=n<4K.
```

Then `0<=b_K(n)<=1`.

## Bounded Feature Vector

For `K<n<4K` define the vector `a_K(n)` in `R^K` by

```text
a_(K,t)(K+i):=1_(t<=i),        1<=i<K,
a_(K,t)(n):=b_K(n),            2K<=n<4K,
                                      1<=t<=K.       (MOVH.3)
```

The anchored suffix vector is exactly

```text
B_K:=sum_(K<n<4K)mu(n)a_K(n),

B_(K,t)
 =beta_K+sum_(i=t)^(K-1)mu(K+i),

H_K=||B_K||_2^2.                                    (MOVH.4)
```

Thus (MOVH.1) becomes

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)||B_K||_2^2<infinity.      (MOVH.5)
```

All arithmetic coefficients in this formulation are ordinary Mobius
values multiplied by deterministic numbers in `[0,1]`.

## Explicit Gram Kernel

Let

```text
L_K(n,m):=<a_K(n),a_K(m)>.
```

For `1<=i,j<K` and `2K<=n,m<4K`,

```text
L_K(K+i,K+j)=min(i,j),
L_K(K+i,n)=i*b_K(n),
L_K(n,m)=K*b_K(n)b_K(m).                             (MOVH.6)
```

Consequently `L_K` is positive semidefinite and
`0<=L_K(n,m)<=K`. The exact expansion is

```text
H_K
 =sum_(n,m in (K,4K))mu(n)mu(m)L_K(n,m)
 =Delta_K+2C_(vec,K),                                (MOVH.7)
```

where

```text
Delta_K
 :=sum_(i=1)^(K-1)i*mu(K+i)^2
   +K*sum_(2K<=n<4K)b_K(n)^2mu(n)^2,

C_(vec,K)
 :=sum_(K<n<m<4K)mu(n)mu(m)L_K(n,m).                 (MOVH.8)
```

Since `b_K(n)<=(4K-n)/(2K)`,

```text
K*sum_(2K<=n<4K)b_K(n)^2
 <=K*sum_(r=1)^(2K)(r/(2K))^2
 =(2K+1)(4K+1)/12.
```

Therefore

```text
0<=Delta_K
 <=K(K-1)/2+(2K+1)(4K+1)/12
 =(14K^2+1)/12
 <=(5/4)K^2.                                        (MOVH.9)
```

The normalized diagonal is dyadically summable for every fixed
`alpha>0`. Hence

```text
R_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)
 K^(-2-alpha)C_(vec,K)<infinity.                    (MOVH.10)
```

This is still a signed criterion. The envelope
`||a_K(n)||_2<=sqrt(K)` and `O(K)` supported coefficients give only
`||B_K||_2=O(K^(3/2))`, hence `H_K=O(K^3)` and normalized scale
`O(K^(1-alpha))`. Absolute values therefore lose one full power.

## Vaughan Before Squaring

The vector identity admits a cleaner handoff than the scalar
off-diagonal expansion. Split the support into the standard intervals
`X<n<=2X` for `X=K,2K`. Put

```text
U_X=V_X=floor(sqrt(X)),

A_X(d)
 :=sum_(bc=d, b<=U_X, c<=V_X)mu(b)mu(c),

D_X(d)
 :=sum_(c|d, c>V_X)mu(c).                            (MOVH.11)
```

Define the finite vectors

```text
V_(I,K,X)
 :=sum_(d<=U_X*V_X)A_X(d)
   sum_(X/d<w<=2X/d)a_K(dw),

V_(II,K,X)
 :=sum_(V_X<d<=2X/U_X)D_X(d)
   sum_(max(U_X,X/d)<w<=2X/d)mu(w)a_K(dw).           (MOVH.12)
```

The inequalities in the inner sums mean integer ranges. Green and Tao's
finite Vaughan identity applies separately to every one of the `K`
coordinates, so no infinite-dimensional extension is being invoked:

```text
sum_(X<n<=2X)mu(n)a_K(n)
 =-V_(I,K,X)+V_(II,K,X).                             (MOVH.13)
```

Set

```text
V_(r,K):=V_(r,K,K)+V_(r,K,2K),       r in {I,II}.
```

The feature at `n=4K` is zero, and the two intervals cover all nonzero
features. Therefore

```text
B_K=-V_(I,K)+V_(II,K),                               (MOVH.14)
```

and the live theorem is exactly

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)
 ||-V_(I,K)+V_(II,K)||_2^2
 <infinity.                                          (MOVH.15)
```

This is Vaughan before squaring: the signed Type I/II cancellation
survives inside the norm. Separate bounds for the two vector norms may
be useful if they save the required power, but they are not equivalent
to (MOVH.15). Even a uniform coordinate estimate `B_(K,t)=o(K)` gives
only `H_K=o(K^3)`, still above the required essentially quadratic scale.

## All-Interval Logarithmic Guard

Matomaki, Shao, Tao, and Teravainen prove, in every interval of length
`H>=X^(5/8+epsilon)`, a Mobius bound of size
`H*log^(-A)X` for arbitrary fixed `A`. This is much stronger than an
unquantified `o(H)` statement, but it is still scalar. Applying it
coordinatewise to the `O(K)` long suffixes and squaring gives at best

```text
K^3*log^(-2A)K,                                     (MOVH.16)
```

up to harmless shorter-suffix and smooth-anchor terms. No logarithmic
choice of `A` turns this into `K^(2+epsilon)`. The cited theorem does
not give the vector correlation across nested suffix endpoints required
by (MOVH.15).

## Open Gate

For every member of one fixed cofinal sequence `alpha_j->0`, prove

```text
sum_(K dyadic)K^(-2-alpha_j)
 ||-V_(I,K)+V_(II,K)||_2^2
 <infinity                                            (MOVH.17)
```

without assuming RH. A viable estimate must save essentially one power
over coefficientwise control while preserving the signed Type I/II
difference and the vector correlation across suffix endpoints. No such
estimate, full Burnol bound, RH, PF-infinity, or `Lambda<=0` conclusion
is supplied here.

Primary source:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
