# Jensen-Window PF Mertens Mixed-Boundary Sine Vaughan Handoff

Date: 2026-07-23

Status: exact mixed-boundary sine and modewise Vaughan reduction with
one open additive-twist gate. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py
```

## Mixed-Boundary Sine Diagonalization

Fix `0<alpha<1`. Lemma 11.22Z gives

```text
R_alpha<infinity
 iff
E_alpha
 :=sum_(K dyadic)K^(-2-alpha)||B_K||_2^2
 <infinity.                                           (MBSV.1)
```

Collapse the future anchor into the terminal coefficient:

```text
x_(K,i):=mu(K+i),       1<=i<K,
x_(K,K):=beta_K,
B_(K,t)=sum_(i=t)^K x_(K,i).                         (MBSV.2)
```

If `C_K(i,j)=min(i,j)`, then

```text
||B_K||_2^2=x_K^T C_K x_K.                           (MBSV.3)
```

The inverse `Q_K=C_K^(-1)` is tridiagonal, with off-diagonal `-1`,
diagonal `2` for `i<K`, and final diagonal `1`. It is the discrete
Laplacian with mixed boundary conditions

```text
v_0=0,                  v_(K+1)=v_K.                 (MBSV.4)
```

For `1<=r<=K`, put

```text
theta_(K,r):=(2r-1)pi/(2K+1),

phi_(K,r)(i)
 :=2/sqrt(2K+1)*sin(i*theta_(K,r)).                  (MBSV.5)
```

The finite sine identity

```text
sum_(i=1)^K
 sin(i*theta_(K,r))*sin(i*theta_(K,q))
 =(2K+1)delta_(r,q)/4
```

shows that the `phi_(K,r)` form an orthonormal basis. Direct substitution
in the interior recurrence and the terminal row gives

```text
C_K phi_(K,r)=lambda_(K,r)phi_(K,r),

lambda_(K,r)
 :=1/[4sin^2(theta_(K,r)/2)].                        (MBSV.6)
```

Define

```text
X_(K,r):=<x_K,phi_(K,r)>.
```

Then

```text
||B_K||_2^2
 =sum_(r=1)^K lambda_(K,r)|X_(K,r)|^2.               (MBSV.7)
```

## Ordinary-Mobius Spectral Features

Write

```text
b_K(n):=((2K)/n)^(1+alpha)*(4K-n)/(2K),
                                      2K<=n<4K.
```

Substituting the exact formula for `beta_K` gives

```text
X_(K,r)
 =2/sqrt(2K+1)*
  [sum_(i=1)^(K-1)mu(K+i)sin(i*theta_(K,r))
   +beta_K*sin(K*theta_(K,r))]

 =sum_(K<n<4K)mu(n)g_(K,r)(n),                       (MBSV.8)
```

where

```text
g_(K,r)(K+i):=phi_(K,r)(i),              1<=i<K,
g_(K,r)(n):=b_K(n)phi_(K,r)(K),          2K<=n<4K.
                                                               (MBSV.9)
```

Every test is compact and satisfies

```text
|g_(K,r)(n)|<=2/sqrt(2K+1).                          (MBSV.10)
```

Thus the spectral modes are normalized ordinary-Mobius additive sine
tests, with the future interval entering only through the terminal
mixed-Neumann amplitude.

## Odd-Frequency Square Function

For

```text
y_(K,r):=(2r-1)pi/(4K+2),
```

the elementary inequalities `2y/pi<=sin(y)<=y` give

```text
4/[pi^2(2r-1)^2]
 <=K^(-2)lambda_(K,r)
 <=9/[4(2r-1)^2].                                   (MBSV.11)
```

Consequently, if

```text
S_alpha
 :=sum_(K dyadic)K^(-alpha)
   sum_(r=1)^K |X_(K,r)|^2/(2r-1)^2,                (MBSV.12)
```

then

```text
(4/pi^2)S_alpha<=E_alpha<=(9/4)S_alpha.              (MBSV.13)
```

In particular,

```text
R_alpha<infinity iff S_alpha<infinity.               (MBSV.14)
```

This is the exact half-odd additive-twist square function.

## Sharp Calibration And Guards

The feature envelope gives `|X_(K,r)|=O(sqrt(K))`; since
`sum_r(2r-1)^(-2)<infinity`, it gives only

```text
K^(-alpha)sum_r |X_(K,r)|^2/(2r-1)^2
 =O(K^(1-alpha)).                                    (MBSV.15)
```

This is one full power short. Davenport's uniform estimate for
Mobius exponential sums, with partial summation for the smooth future
weight, improves this only to

```text
O_A(K^(1-alpha)log^(-2A)K).                          (MBSV.16)
```

No fixed logarithmic saving is dyadically summable here. The
all-interval `H*log^(-A)X` theorem has the same scalar limitation.
Generic Parseval or large-sieve bounds likewise do not create
Mobius-specific cancellation.

The sharp simple sufficient condition is transparent. If, uniformly in
`1<=r<=K`, the bracket in (MBSV.8) satisfies

```text
|sum_(i=1)^(K-1)mu(K+i)sin(i*theta_(K,r))
 +beta_K*sin(K*theta_(K,r))|
 =O_epsilon(K^(1/2+epsilon)),                        (MBSV.17)
```

then `|X_(K,r)|=O_epsilon(K^epsilon)`, and (MBSV.12)
converges for every `2epsilon<alpha`. This square-root-plus-epsilon
bound is a conditional calibration, not an unconditional theorem.

## Vaughan Mode By Mode

Split the support into `X<n<=2X` for `X=K,2K`, and put

```text
U_X=V_X=floor(sqrt(X)),

A_X(d)
 :=sum_(bc=d,b<=U_X,c<=V_X)mu(b)mu(c),

D_X(d)
 :=sum_(c|d,c>V_X)mu(c).                             (MBSV.18)
```

Define

```text
T_(I,K,r,X)
 :=sum_(d<=U_X*V_X)A_X(d)
   sum_(X/d<w<=2X/d)g_(K,r)(dw),

T_(II,K,r,X)
 :=sum_(V_X<d<=2X/U_X)D_X(d)
   sum_(max(U_X,X/d)<w<=2X/d)mu(w)g_(K,r)(dw).
                                                               (MBSV.19)
```

Green and Tao's finite Vaughan identity gives, after summing the two
intervals,

```text
X_(K,r)=-T_(I,K,r)+T_(II,K,r),                       (MBSV.20)
```

where `T_(q,K,r)` is the sum of the `X=K` and `X=2K` pieces. Therefore
the live criterion is exactly

```text
sum_(K dyadic)K^(-alpha)
 sum_(r=1)^K
 |-T_(I,K,r)+T_(II,K,r)|^2/(2r-1)^2
 <infinity.                                          (MBSV.21)
```

Vaughan is applied before the mode square. Replacing the signed
difference by separate absolute Type I and Type II bounds is not an
equivalent step.

## Open Gate

Prove (MBSV.21) for every member of one fixed cofinal sequence
`alpha_j->0` without assuming RH. A viable estimate must deliver
Mobius-specific square-root cancellation on average in the half-odd
mode measure, or an equally strong joint Type I/II bilinear bound.
The square-root estimate, signed power gain, full Burnol bound, RH,
PF-infinity, and `Lambda<=0` all remain open.

Primary sources:

- Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
