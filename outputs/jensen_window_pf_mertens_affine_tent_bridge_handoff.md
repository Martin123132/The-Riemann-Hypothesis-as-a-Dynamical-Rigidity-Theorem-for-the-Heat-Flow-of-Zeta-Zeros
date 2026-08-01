# Jensen-Window PF Mertens Affine-Tent/Bridge Handoff

Date: 2026-07-23

Status: exact affine scale-localization and unified local Vaughan
reduction with one open signed gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_affine_tent_bridge_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_affine_tent_bridge_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_affine_tent_bridge_handoff.py
```

## Affine Mean As A Scale Difference

Fix `0<alpha<1` and set

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.
```

For dyadic `K`, write

```text
u_K
 :=sum_(N=K)^(2K-1)r_N
  =sqrt(K)c_(K,0).                                    (MATBH.1)
```

The affine mean series from Lemma 11.22X is

```text
M_alpha
 =sum_(K dyadic)K^(alpha-1)|u_K|^2.                   (MATBH.2)
```

Counting the tails containing each coefficient gives

```text
u_K=sum_n q_n w_K(n),

w_K(n):=min((n-K)_+,K).                               (MATBH.3)
```

The plateau of `w_K` is the nonlocal part. Apply the forced scale
difference

```text
v_K:=u_K-(1/2)u_(2K).                                 (MATBH.4)
```

Define the compact tent

```text
W_K(n)
 :=
  n-K,          K<n<=2K,
  (4K-n)/2,     2K<n<4K,
  0,            otherwise.                            (MATBH.5)
```

Since `W_K=w_K-(1/2)w_(2K)`, the infinite tail cancels exactly:

```text
v_K
 =sum_(K<n<4K)mu(n)n^(-1-alpha)W_K(n).                (MATBH.6)
```

Put

```text
p:=(alpha-1)/2,
d_K:=K^p u_K,
e_K:=K^p v_K,
rho_alpha:=2^(-(1+alpha)/2).                          (MATBH.7)
```

Then

```text
e_K=d_K-rho_alpha*d_(2K).                             (MATBH.8)
```

The coefficient envelope gives

```text
|d_K|=O_alpha(K^((1-alpha)/2)),

rho_alpha^m d_(2^m K)=O_alpha(2^(-alpha*m))->0.       (MATBH.9)
```

Thus iteration of (MATBH.8) is legitimate without assuming the target:

```text
d_K=sum_(m>=0)rho_alpha^m e_(2^m K).                 (MATBH.10)
```

The unilateral shift and its geometric inverse give

```text
(1-rho_alpha)||d||_2
 <=||e||_2
 <=(1+rho_alpha)||d||_2.                              (MATBH.11)
```

Consequently

```text
M_alpha<infinity
 iff
E_aff:=sum_(K dyadic)|e_K|^2<infinity.                (MATBH.12)
```

This is an exact localization, not an estimate for `M_alpha`.

## Signed Affine Tent

Expanding the compact square,

```text
|e_K|^2=D_(aff,K)+2C_(aff,K),                         (MATBH.13)

D_(aff,K)
 :=K^(alpha-1)
   sum_(K<n<4K)
   mu(n)^2 n^(-2-2alpha)W_K(n)^2,

C_(aff,K)
 :=K^(alpha-1)
   sum_(K<n<m<4K)
   mu(n)mu(m)(nm)^(-1-alpha)W_K(n)W_K(m).             (MATBH.14)
```

Since there are fewer than `3K` terms and `W_K<=K`,

```text
0<=D_(aff,K)<=3K^(-alpha),

sum_(K dyadic)D_(aff,K)<infinity.                     (MATBH.15)
```

It follows that

```text
M_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)C_(aff,K)<infinity.          (MATBH.16)
```

The affine mean has therefore become a compact signed off-diagonal
problem on `(K,4K)`.

## Unified Local Kernel

Retain the Brownian-bridge energy from Lemma 11.22X:

```text
E_(br,K)
 :=K^(-2-alpha)||[S_K]||_K^2
  =D_(br,K)+2C_(br,K),                                (MATBH.17)

D_(br,K):=K^(-2-alpha)Delta_K,
C_(br,K):=K^(-2-alpha)C_K.
```

Its diagonal satisfies

```text
0<=D_(br,K)<=K^(-alpha)/6.                            (MATBH.18)
```

Define

```text
E_loc
 :=sum_(K dyadic)(|e_K|^2+E_(br,K)).                  (MATBH.19)
```

Lemmas 11.22X and (MATBH.12) give the exact equivalence

```text
R_alpha<infinity iff E_loc<infinity.                  (MATBH.20)
```

For `h>=1`, define

```text
G_(K,h)(n)
 :=
 K^(alpha-1)[n(n+h)]^(-1-alpha)
 W_K(n)W_K(n+h)
 1_(K<n<n+h<4K)

 +K^(-3-alpha)(n-K)(2K-n-h)
 1_(K<n<n+h<2K).                                     (MATBH.21)
```

The first line is the affine-tent kernel and the second is the
normalized Brownian-bridge kernel. They have the same scale:

```text
0<=G_(K,h)(n)<=(5/4)K^(-1-alpha),                     (MATBH.22)
```

and both vanish unless `K<n<n+h<4K`. For a finite dyadic cutoff,

```text
E_loc(J)=D_loc(J)+2O_loc(J),                          (MATBH.23)

O_loc(J)
 :=sum_(K<=2^J)
   sum_h sum_n
   mu(n)mu(n+h)G_(K,h)(n).
```

The diagonal obeys

```text
D_(loc,K)
 <=(3+1/6)K^(-alpha)
 =(19/6)K^(-alpha),                                  (MATBH.24)
```

so it is dyadically summable. Hence

```text
R_alpha<infinity
 iff
sup_J O_loc(J)<infinity.                              (MATBH.25)
```

This replaces the two open pieces of Lemma 11.22X by one compact signed
correlation criterion.

## Two-Interval Vaughan Handoff

Split the base variable into

```text
K<n<=2K,
2K<n<4K.                                              (MATBH.26)
```

For `a in {0,1}`, let

```text
g_(K,h,a)(n)
 :=mu(n+h)G_(K,h)(n)
   1_(2^a K<n<=2^(a+1)K).                             (MATBH.27)
```

Each piece now has the standard `(X,2X]` form. Applying the finite
Vaughan identity of Green and Tao separately to the two intervals and
then summing over `h` gives

```text
O_(loc,K)=-TI_(loc,K)+TII_(loc,K).                    (MATBH.28)
```

Therefore the single live criterion is

```text
sup_J sum_(K dyadic<=2^J)
 (-TI_(loc,K)+TII_(loc,K))
 <infinity.                                           (MATBH.29)
```

The signed difference is essential. Separate absolute values see
`O(K^2)` shift/base pairs of size `O(K^(-1-alpha))` and restore the
divergent scale `K^(1-alpha)`. The natural `o(K^2)` averaged-Chowla
pair scale likewise gives only `o(K^(1-alpha))`, one full power above
the summable diagonal calibration `K^(-alpha)`.

## Positive Compact Tail Square

There is a simpler Cholesky form before transferring the centered path
to ordinary Mertens increments. Retain

```text
T_(K,m)
 :=sum_(i=1)^m q_(K+i),                    0<=m<K,

Z_K
 :=K*q_(2K)
   +sum_(2K<n<4K)(4K-n)q_n/2.                         (MATBH.30)
```

The endpoint `q_(2K)` belongs to the affine tent, while the centered
length-`K` path uses only `q_(K+1),...,q_(2K-1)`. Thus

```text
v_K
 =sum_(i=1)^(K-1)i*q_(K+i)+Z_K.                       (MATBH.31)
```

The affine rank-one kernel exactly fills the constant direction removed
by the Brownian bridge:

```text
K^(alpha-1)i*j
 +K^alpha[min(i,j)-i*j/K]
 =K^alpha min(i,j),                 1<=i,j<K.          (MATBH.32)
```

Finite variance completion therefore gives the exact positive identity

```text
|e_K|^2+K^alpha||[T_K]||_K^2

 =K^alpha sum_(t=1)^K
   |Z_K/K+sum_(i=t)^(K-1)q_(K+i)|^2.                  (MATBH.33)
```

For `t=K` the suffix sum is empty. DCT Parseval, dyadic weight
comparison, and the invertible scale filter now yield the alternative
one-gate criterion

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^alpha
 sum_(t=1)^K
 |Z_K/K+sum_(i=t)^(K-1)q_(K+i)|^2
 <infinity.                                           (MATBH.34)
```

Every coefficient in (MATBH.34) lies in `(K,4K)`. This positive local
tail-lattice square is equivalent to the signed-kernel target, but it
does not estimate the Mobius suffixes.

## Ordinary-Mobius Anchored Form

The power weights on the current-block suffixes can be removed exactly.
Put `s=1+alpha`, `z_K=Z_K/K`, and

```text
beta_K
 :=(2K)^s*z_K
  =mu(2K)
   +sum_(2K<n<4K)
     ((2K)/n)^s*(4K-n)/(2K)*mu(n).                   (MATBH.35)
```

Every coefficient in `beta_K` lies in `[0,1]`, and its support is
`2K<=n<4K`. For `1<=i<K` set

```text
a_i:=(K+i)^(-s),       x_i:=mu(K+i),
a_K:=(2K)^(-s),        x_K:=beta_K,

A_(K,t):=sum_(i=t)^K a_i*x_i
        =z_K+sum_(i=t)^(K-1)(K+i)^(-s)mu(K+i),

B_(K,t):=sum_(i=t)^K x_i
        =beta_K+sum_(i=t)^(K-1)mu(K+i).              (MATBH.36)
```

Writing `A_K` and `B_K` for these length-`K` vectors, finite Abel
summation gives the mutually inverse triangular maps

```text
A_(K,t)
 =a_t*B_(K,t)
  +sum_(i=t+1)^K(a_i-a_(i-1))*B_(K,i),

B_(K,t)
 =a_t^(-1)*A_(K,t)
  +sum_(i=t+1)^K(a_i^(-1)-a_(i-1)^(-1))*A_(K,i).
                                                            (MATBH.37)
```

The first matrix has absolute row sums at most `2K^(-s)` and column
sums at most `3K^(-s)`. Its inverse has absolute row sums at most
`4K^s` and column sums at most `8K^s`. Schur's test therefore proves
the uniform two-sided estimate

```text
(1/(4sqrt(2)))*K^(-s)||B_K||_2
 <=||A_K||_2
 <=sqrt(6)*K^(-s)||B_K||_2.                          (MATBH.38)
```

Combining (MATBH.34) with (MATBH.38) yields the ordinary-Mobius
anchored criterion

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)
 sum_(t=1)^K
 |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2
 <infinity.                                           (MATBH.39)
```

Thus, if

```text
H_K
 :=sum_(t=1)^K
   |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2
 =O_epsilon(K^(2+epsilon)),       epsilon<alpha,      (MATBH.40)
```

then the criterion holds at that `alpha`. The coefficientwise bound
gives only `H_K=O(K^3)`, so cancellation sufficient to save almost one
full power is still required.

## Quantifier And Critical Guards

The strongest directly relevant short-interval theorems do not have the
required calibration. Matomaki and Teravainen prove `o(H)` cancellation
in every interval of length `H=x^theta` for `theta>0.55`; an unquantified
`o(H)` bound is still far above the square-root scale needed after
summing the suffix squares. Matomaki-Radziwill and the newer
higher-uniformity results give strong logarithmic cancellation outside
small exceptional sets of base points. The present family uses one fixed
dyadic base point at each scale, and no cited theorem de-averages those
exceptional sets along that sequence.

The parameter boundary is also strict:

```text
sum_(K dyadic)K^(-alpha)<infinity
```

holds for every fixed `alpha>0` and fails at `alpha=0`. The cofinal
criterion asks separately for each positive `alpha_j->0`; it does not
permit setting `alpha=0` or passing to a bound uniform in `alpha`.

## Open Gate

The scale filter retires the nonlocal affine tail and the quotient
retires the Abel/low-DCT commutator. Neither operation supplies the
missing arithmetic cancellation. The remaining task is:

```text
Prove (MATBH.29), equivalently (MATBH.34) and (MATBH.39), for every
member of one fixed cofinal sequence alpha_j->0
without assuming RH.                                  (MATBH.41)
```

A viable argument must preserve cancellation between the Type I and
Type II pieces and across dyadic scales. Rowwise absolute values,
post-collapse Vaughan, generic operator bounds, and direct promotion of
averaged Chowla remain insufficient. A signed Type I/II power gain, the
full Burnol bound, RH, PF-infinity, and `Lambda <= 0` all remain open.

Primary sources:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
- Matomaki and Teravainen, *On the Mobius Function in All Short
  Intervals*:
  https://arxiv.org/abs/1911.09076
- Matomaki and Radziwill, *Multiplicative Functions in Short
  Intervals*:
  https://arxiv.org/abs/1501.04585
- Matomaki, Radziwill, Shao, Tao, and Teravainen, *Higher Uniformity
  of Arithmetic Functions in Short Intervals II*:
  https://doi.org/10.1007/s00222-026-01408-6
