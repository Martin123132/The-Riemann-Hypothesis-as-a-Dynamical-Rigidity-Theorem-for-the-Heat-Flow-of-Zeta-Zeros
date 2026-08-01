# Jensen-Window PF Mertens Planar Abel Handoff

Date: 2026-07-24

Status: exact planar Abel identity and logarithmic mixed-kernel
variation with one open arithmetic energy gate.
This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_planar_abel_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_planar_abel_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_abel_handoff.py
```

## Exact Planar Abel Identity

Retain the notation of Corollary 11.22Z.4. Thus

```text
N=2K+1,
R=ceil(K^(1-alpha/2)),
rho_K=pi^2*K/N^2,
O_(alpha,K)=sum_(K<n<m<4K)mu(n)mu(m)G(n,m).
```

Use the rectangle

```text
a_K:=K+1<=n<=b_K:=4K-2,
1<=h<=H_K:=3K-2.
```

Extend both the arithmetic array and the kernel by zero:

```text
a_K(n,h):=
  mu(n)mu(n+h),  n+h<=4K-1,
  0,              n+h>4K-1,

Wtilde_K(n,h):=
  G(n,n+h),       n+h<=4K-1,
  0,              n+h>4K-1.                  (PAH.1)
```

Also set `Wtilde_K(b_K+1,h)=Wtilde_K(n,H_K+1)=0`. Define

```text
S_K(x,H)
 :=sum_(n=K+1)^x sum_(h=1)^H a_K(n,h),        (PAH.2)

Delta Wtilde_K(n,h)
 :=Wtilde_K(n,h)-Wtilde_K(n+1,h)
   -Wtilde_K(n,h+1)+Wtilde_K(n+1,h+1).        (PAH.3)
```

Two finite Abel summations give the exact identity

```text
O_(alpha,K)
 =sum_(n=K+1)^(4K-2)sum_(h=1)^(3K-2)
   S_K(n,h)Delta Wtilde_K(n,h).               (PAH.4)
```

No absolute value over the base point or shift has been taken in
(PAH.4).

## Current-Current Curvature

Put

```text
C_R(u)=sum_(r=1)^R cos((2r-1)u)/(2r-1)^2,

D_R(u):=sum_(r=1)^R cos((2r-1)u)
       =sin(2Ru)/(2sin u),

L_R:=integral_0^pi |D_R(u)|du.                (PAH.5)
```

On `(0,pi/2]`,

```text
|D_R(u)|<=min(R,pi/(4u)).
```

Symmetry and a split at `u=pi/(4R)` therefore give

```text
L_R<=(pi/2)(1+log(2R)).                       (PAH.6)
```

For `n=K+i`, `m=n+h=K+j`, and `s=i+j=2i+h`, the
current-current formula is

```text
Wtilde_K(n,h)
 =(2/N)[C_R(h*pi/N)-C_R(s*pi/N)].
```

The first term disappears under the mixed difference. Whenever
`i+h<=K-2`,

```text
Delta Wtilde_K(n,h)
 =(2/N)integral_0^(pi/N)integral_0^(2pi/N)
   D_R(s*pi/N+u+v)dvdu.                       (PAH.7)
```

For each `s` there are at most `K` pairs `(i,h)`. For fixed `v`, the
`pi/N` strips indexed by `s` are disjoint. Hence the total interior
current-current contribution to mixed variation is at most

```text
4*pi*K*L_R/N^2.                               (PAH.8)
```

The logarithm in (PAH.6) is the ordinary `L1` cost of the finite
Dirichlet kernel. No arithmetic estimate has entered.

## Transition And Future Curvature

For real `2K<=x<=4K`, write

```text
b_K(x)=(2K)^alpha(4K-x)x^(-1-alpha).
```

It is decreasing and strictly convex because

```text
b_K'(x)
 =(2K)^alpha[
   -4K(1+alpha)x^(-2-alpha)+alpha*x^(-1-alpha)
  ]<0,

b_K''(x)
 =(2K)^alpha(1+alpha)x^(-3-alpha)
  [4K(2+alpha)-alpha*x]>0.                    (PAH.9)
```

Set `b_x=b_K(x)` on `[2K,4K]`, extend it by zero for `x>=4K`, and
put `d_x=b_x-b_(x+1)`. Then

```text
d_x>=d_(x+1)>=0,
sum_(x>=2K)d_x=1,
d_(2K)<=(2+alpha)/(2K)<3/(2K).               (PAH.10)
```

Also retain `p_i=P_(K,R)(i,K)` and `c=P_(K,R)(K,K)`. Then

```text
0=p_0<p_1<...<p_K=c<=rho_K.                  (PAH.11)
```

The sine-integral boundary strips are disjoint, so

```text
sum_(i=1)^(K-2)|P(i,K)-P(i,K-1)|
 <pi^2/(2N).                                  (PAH.12)
```

All remaining mixed cells fall into four explicit families.

At the second-argument transition `m=2K-1`,

```text
Delta W
 =P(i,K-1)-p_i-p_(i+1)d_(2K),

sum_i|Delta W|
 <pi^2/(2N)+(3/2)rho_K.                      (PAH.13)
```

For current `n=K+i` and future `m>=2K`,

```text
Delta W=p_i d_m-p_(i+1)d_(m+1),

sum_(i,m)|Delta W|<=(5/2)rho_K.              (PAH.14)
```

At the first-argument transition `n=2K-1`,

```text
Delta W=p_(K-1)d_m-c*d_(m+1),

sum_m|Delta W|<=2rho_K.                       (PAH.15)
```

Finally, for `2K<=n<m`,

```text
Delta W
 =c[b_n(d_m-d_(m+1))+d_n d_(m+1)]>=0,

sum_(n,m)Delta W<=2rho_K.                     (PAH.16)
```

The joined treatment at `m=2K` is essential. Extending the
current-current and current-future blocks separately by zero would
manufacture a false `O(1)` interface variation that is absent from
the actual kernel.

## Logarithmic Mixed Variation

Define

```text
V_(alpha,K)
 :=sum_(n=K+1)^(4K-2)sum_(h=1)^(3K-2)
   |Delta Wtilde_K(n,h)|.
```

Adding (PAH.8) and (PAH.13)-(PAH.16) gives

```text
V_(alpha,K)
 <=4*pi*K*L_R/N^2+pi^2/(2N)+8rho_K.
```

Using (PAH.6), `N>2K`, and `rho_K<pi^2/(4K)` yields the convenient
uniform bound

```text
V_(alpha,K)
 <3*pi^2*(1+log(2R))/K.                       (PAH.17)
```

Thus, with

```text
M_K:=max_(K+1<=x<=4K-2,1<=H<=3K-2)|S_K(x,H)|,
```

(PAH.4) implies

```text
|O_(alpha,K)|
 <3*pi^2*(1+log(2R))*M_K/K.                  (PAH.18)
```

This is the first joint base/shift Abel interface: separate
square-root estimates are no longer summed over the other axis.

## Curvature-Weighted Energy Gate

The maximal condition in (PAH.18) is stronger than necessary. Define

```text
E_(alpha,K)
 :=sum_(n,h)|Delta Wtilde_K(n,h)||S_K(n,h)|^2. (PAH.19)
```

Weighted Cauchy-Schwarz applied only after (PAH.4) gives

```text
|O_(alpha,K)|^2
 <=V_(alpha,K)E_(alpha,K).                    (PAH.20)
```

Consequently, for one fixed `alpha>0`, the conditional estimate

```text
E_(alpha,K)=O_epsilon(K^(1+epsilon))          (PAH.21)
```

on dyadic `K`, for every `epsilon>0`, implies absolute summability of
`sum_K K^(-alpha)O_(alpha,K)`: choose `epsilon<2alpha` in
(PAH.17) and (PAH.20).

The stronger planar maximum

```text
M_K=O_epsilon(K^(1+epsilon))                  (PAH.22)
```

also suffices, now by choosing `epsilon<alpha` in (PAH.18). Both
(PAH.21) and (PAH.22) have the square-root scale for a planar region
with `O(K^2)` terms. Neither estimate is proved.

The energy gate (PAH.21) is the narrower target: it samples prefixes
only against the mixed-curvature measure of the actual kernel.

## Strength And Countermodel Guards

The terminal prefix is

```text
S_K(4K-2,3K-2)
 =1/2[A_K^2-Q_K],

A_K:=sum_(K<n<4K)mu(n),
Q_K:=sum_(K<n<4K)mu(n)^2.                    (PAH.23)
```

Therefore the maximal condition (PAH.22) already implies
`A_K=O_epsilon(K^(1/2+epsilon))`. It contains RH-scale block-Mertens
input and is not an unconditional shortcut. The curvature energy
(PAH.21) does not isolate that single corner in the same way.

Separate row and column square-root prefixes still do not imply a
planar square-root estimate. Let `K=L^2` and form a `K` by `K`
zero-one matrix from `L` diagonal all-one blocks of size `L` by `L`.
Every row and column prefix is at most `L=sqrt(K)`, while the full
rectangle sum is

```text
L^3=K^(3/2).                                  (PAH.24)
```

This abstract countermodel blocks promotion from the two
one-dimensional interfaces of Corollary 11.22Z.4.

Matomaki-Radziwill-Tao averaged Chowla has terminal scale `o(HX)`.
It neither provides the varying two-parameter prefixes in (PAH.2)
nor a fixed power saving from `K^2` to `K^(1+epsilon)`. The
all-interval Mobius/nilsequence theorem supplies logarithmic decay
for scalar linear tests, not the quadratic curvature energy
(PAH.19). Both are therefore power- or quantifier-short here.

## Vaughan Form And Open Gate

Let `T_(x,H)` be the finite triangular band operator whose quadratic
form is `S_K(x,H)`. Inserting the finite Vaughan identity
`mu=-u_I+u_II` before any square gives

```text
S_K(x,H)
 =<u_I,T u_I>+<u_II,T u_II>
  -<u_I,T u_II>-<u_II,T u_I>.                (PAH.25)
```

The operator need not be self-adjoint, so its two cross terms must
both remain. This is a rational-frequency, large-sieve, or
Hilbert-valued Vaughan interface; it is not permission to estimate
the four terms separately.

The exact surviving criterion is

```text
sup_J sum_(K dyadic<=2^J)K^(-alpha)
 sum_(n,h)S_K(n,h)Delta Wtilde_K(n,h)
 <infinity                                    (PAH.26)
```

for every member of one fixed cofinal positive-`alpha` sequence.
The concrete sufficient theorem target is (PAH.21), or a direct
signed estimate of (PAH.26) that is weaker.

The planar Abel identity, logarithmic mixed-variation estimate, and
Vaughan expansion are exact. The curvature-energy estimate, weighted
joint Mobius gain, full Burnol bound, RH, PF-infinity, and
`Lambda<=0` remain open.

The separate finite float64 diagnostic

```text
outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md
work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py
```

tests the actual-Mobius quantities through `K=1024`. It is a
falsification scout only and is not used in the proof of (PAH.17).

Primary sources:

- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
