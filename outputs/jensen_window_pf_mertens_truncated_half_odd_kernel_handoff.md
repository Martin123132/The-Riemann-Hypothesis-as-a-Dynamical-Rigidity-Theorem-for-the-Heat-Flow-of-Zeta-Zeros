# Jensen-Window PF Mertens Truncated Half-Odd Kernel Handoff

Date: 2026-07-24

Status: exact truncated half-odd kernel, uniformly summable diagonal,
and signed off-diagonal/Vaughan reduction with one open arithmetic
gate. This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py
```

## Exact Truncated Kernel

Fix `0<alpha<1`, let `K` be dyadic, and put

```text
R=R_(alpha,K):=ceil(K^(1-alpha/2)),
N:=2K+1,
q_r:=2r-1,
theta_r:=q_r*pi/N,

phi_r(i):=2/sqrt(N)*sin(i*theta_r),       1<=i<=K.
                                                        (THOK.1)
```

Corollary 11.22Z.2 proves that the reciprocal-tail criterion is
equivalent to summability of

```text
L_(alpha,K)
 :=sum_(r=1)^R |X_(K,r)|^2/q_r^2.          (THOK.2)
```

Define

```text
P_(K,R)(i,j)
 :=sum_(r=1)^R phi_r(i)phi_r(j)/q_r^2.      (THOK.3)
```

For every `x in R^K`,

```text
sum_(r=1)^R |<x,phi_r>|^2/q_r^2
 =x^T P_(K,R)x.                             (THOK.4)
```

Put

```text
C_R(u):=sum_(r=1)^R cos(q_r*u)/q_r^2,
S_R(u):=sum_(r=1)^R sin(q_r*u)/q_r.
```

Product-to-sum gives the exact finite formula

```text
P_(K,R)(i,j)
 =2/N[
    C_R((i-j)pi/N)-C_R((i+j)pi/N)
   ].                                       (THOK.5)
```

Moreover,

```text
C_R'(u)=-S_R(u),

S_R'(u)
 =sum_(r=1)^R cos(q_r*u)
 =sin(2R*u)/(2sin u),                       (THOK.6)
```

where the last quotient is interpreted continuously. Hence

```text
P_(K,R)(i,j)
 =2/N integral_(|i-j|pi/N)^((i+j)pi/N)
      S_R(u)du.                              (THOK.7)
```

The sine polynomial is strictly positive:

```text
S_R(u)>0,                         0<u<pi.    (THOK.8)
```

Indeed, `S_R(pi-u)=S_R(u)`. On `(0,pi/2]`, write

```text
S_R(u)=integral_0^u sin(2R*v)/(2sin v)dv.
```

Pair each negative half-wave with the preceding positive half-wave.
The amplitude `1/sin v` is decreasing there, so every completed pair
is nonnegative; a terminal partial half-wave has the same property.
The first positive half-wave makes the inequality strict. Symmetry
handles `(pi/2,pi)`.

It follows from (THOK.7) that

```text
P_(K,R)(i,j)>0                    (THOK.9)
```

entrywise. Independently, (THOK.3) proves that `P_(K,R)` is positive
semidefinite, has rank `R`, and is nested in Loewner order:

```text
P_(K,R+1)-P_(K,R)
 =phi_(R+1)phi_(R+1)^T/q_(R+1)^2 >=0.       (THOK.10)
```

## Completion To The Min Kernel

For `0<=u<=pi`, the full odd cosine series is

```text
sum_(r>=1)cos((2r-1)u)/(2r-1)^2
 =pi^2/8-pi*u/4.                             (THOK.11)
```

Inserting (THOK.11) into (THOK.5) gives the exact infinite completion

```text
P_(K,infinity)(i,j)
 =pi^2/N^2*min(i,j).                         (THOK.12)
```

If `C_K(i,j)=min(i,j)`, the omitted tail is an absolutely convergent
sum of rank-one positive semidefinite matrices. Therefore

```text
0<=P_(K,R)<=pi^2/N^2*C_K                     (THOK.13)
```

in Loewner order. At the full finite rank `R=K`, comparison in the
common sine eigenbasis gives

```text
4/N^2*C_K
 <=P_(K,K)
 <=pi^2/N^2*C_K.                             (THOK.14)
```

The lower bound uses `sin(theta/2)>=theta/pi`; the upper bound uses
`sin(theta/2)<=theta/2`. No lower comparison with `C_K` is possible
when `R<K`, because `P_(K,R)` then has rank `R`.

Orthonormality also gives

```text
tr P_(K,R)
 =sum_(r=1)^R 1/q_r^2
 <pi^2/8,

0<P_(K,R)(i,i)<=pi^2*i/N^2.                  (THOK.15)
```

## Ordinary-Mobius Pullback

Retain

```text
b_K(n)
 :=((2K)/n)^(1+alpha)*(4K-n)/(2K),
                                      2K<=n<4K.
```

For `K<n<4K`, define

```text
tau_K(n):=
  n-K,                              K<n<2K,
  K,                               2K<=n<4K,

w_K(n):=
  1,                                K<n<2K,
  b_K(n),                           2K<=n<4K.
```

Then the ordinary-Mobius feature from Corollary 11.22Z.1 is exactly

```text
g_(K,r)(n)=w_K(n)phi_r(tau_K(n)).             (THOK.16)
```

Pull the truncated kernel back to the ordinary support:

```text
G_(alpha,K,R)(n,m)
 :=w_K(n)w_K(m)
   P_(K,R)(tau_K(n),tau_K(m)).                 (THOK.17)
```

This kernel is positive semidefinite and entrywise positive. The low
energy has the exact finite expansion

```text
L_(alpha,K)
 =sum_(K<n,m<4K)
   mu(n)mu(m)G_(alpha,K,R)(n,m)

 =D_(alpha,K)+2O_(alpha,K),                    (THOK.18)

D_(alpha,K)
 :=sum_(K<n<4K)mu(n)^2G_(alpha,K,R)(n,n),

O_(alpha,K)
 :=sum_(K<n<m<4K)mu(n)mu(m)G_(alpha,K,R)(n,m).
```

The current diagonal is less than `pi^2/8` by (THOK.15). On the future
interval,

```text
sum_(2K<=n<4K)b_K(n)^2
 <=sum_(ell=1)^(2K)(ell/(2K))^2
 =(2K+1)(4K+1)/(12K).                         (THOK.19)
```

Using `P_(K,R)(K,K)<=pi^2*K/N^2`, the future diagonal is less than
`pi^2/6`. Thus

```text
0<=D_(alpha,K)<7pi^2/24                       (THOK.20)
```

uniformly in `alpha`, `K`, and `R`. Its normalized dyadic series is
automatically finite. Combining (THOK.18)-(THOK.20) with Corollary
11.22Z.2 proves

```text
R_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)
 K^(-alpha)O_(alpha,K)<infinity.              (THOK.21)
```

This is an exact low-mode, ordinary-Mobius, signed off-diagonal
criterion. Kernel positivity does not estimate that signed sum.

## Endpoint Abel Structure

Let

```text
p_i:=P_(K,R)(i,K),                    0<=i<=K,
p_0:=0.
```

The exact endpoint phase gives

```text
p_i
 =4/N sum_(r=1)^R
   (-1)^(r-1)sin(i*theta_r)cos(theta_r/2)/q_r^2.
                                                        (THOK.22)
```

The integral form is

```text
p_i
 =2/N integral_((K-i)pi/N)^((K+i)pi/N)
      S_R(u)du.                              (THOK.23)
```

Since `S_R>0`, these intervals are strictly nested:

```text
0=p_0<p_1<...<p_K.                            (THOK.24)
```

Put

```text
A_(K,t):=sum_(i=t)^(K-1)mu(K+i).
```

Finite suffix Abel summation now gives

```text
sum_(i=1)^(K-1)mu(K+i)p_i
 =sum_(t=1)^(K-1)
   A_(K,t)(p_t-p_(t-1)),                     (THOK.25)
```

and therefore

```text
|sum_(i=1)^(K-1)mu(K+i)p_i|
 <=pi^2*K/N^2*max_t|A_(K,t)|.                (THOK.26)
```

Thus the current-anchor cross term is an exact positive Abel average,
not an uncontrolled collection of modes. This geometry does not by
itself supply the missing power: Davenport-scale logarithmic bounds
remain one power short. If both `|beta_K|` and
`max_t|A_(K,t)|` were conditionally
`O_epsilon(K^(1/2+epsilon))`, the cross term would be
`O_epsilon(K^(2epsilon))`; both inputs are unproved and the other
off-diagonal blocks would still remain.

## Signed Vaughan Kernel Form

For `X in {K,2K}`, retain `U_X=V_X=floor(sqrt(X))` and the coefficients
`A_X(d),D_X(d)` from Lemma 11.22Z. Collect the finite Vaughan
incidences as

```text
u_(I,X)(n)
 :=1_(X<n<=2X)
   sum_(d|n,d<=U_X*V_X)A_X(d),

u_(II,X)(n)
 :=1_(X<n<=2X)
   sum_(d|n,V_X<d<=2X/U_X,n/d>U_X)
   D_X(d)mu(n/d).                              (THOK.27)
```

Sum the two disjoint intervals and write the resulting coefficients as
`u_(I,K),u_(II,K)`. The feature at `4K` is zero. Green and Tao's finite
identity gives

```text
mu(n)=-u_(I,K)(n)+u_(II,K)(n)                 (THOK.28)
```

on the feature support, and

```text
T_(q,K,r)
 =sum_(K<n<4K)u_(q,K)(n)g_(K,r)(n),
                                      q in {I,II}.
```

Consequently the entire truncated square is

```text
L_(alpha,K)
 =<u_I,G u_I>
  +<u_II,G u_II>
  -2<u_I,G u_II>.                             (THOK.29)
```

The cross term is part of the exact identity. Replacing (THOK.29) by
separate absolute Type I and Type II estimates is not equivalent and
may discard the only available cancellation.

## Open Gate

For every member of one fixed cofinal sequence `alpha_j->0`, prove the
weighted off-diagonal bound (THOK.21), equivalently the joint signed
Vaughan Gram estimate (THOK.29), without assuming RH. The diagonal,
high modes, kernel algebra, and endpoint Abel geometry are closed.
The Mobius-specific off-diagonal gain, full Burnol bound, RH,
PF-infinity, and `Lambda<=0` remain open.

Primary sources:

- Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
