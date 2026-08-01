# Jensen-Window PF Mertens Weighted-Prefix/Affine-Defect Reduction

Date: 2026-07-23

Status: exact weighted-prefix, first-block affine-defect, and Vaughan
handoff reduction with two open arithmetic gates. This is not a proof of
RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json
python work/rh_compute/scripts/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
```

## Weighted Prefix Is The Same Hardy Energy

Fix `0<alpha<1` and set

```text
M(x):=sum_(n<=x)mu(n),
A_alpha(x):=sum_(n<=x)mu(n)n^(-alpha),

I_alpha:=integral_1^infinity M(x)^2*x^(-2-alpha)dx,
J_alpha:=integral_1^infinity A_alpha(x)^2*x^(alpha-2)dx,

beta:=(1+alpha)/2,
lambda:=(1-alpha)/2.                                      (MWPADR.1)
```

Stieltjes integration by parts in both directions gives

```text
A_alpha(x)
 =x^(-alpha)M(x)
  +alpha*integral_1^x M(u)u^(-alpha-1)du,                  (MWPADR.2)

M(x)
 =x^alpha*A_alpha(x)
  -alpha*integral_1^x A_alpha(u)u^(alpha-1)du.             (MWPADR.3)
```

With `x=e^t`, define

```text
m(t):=e^(-beta*t)M(e^t),
a(t):=e^(-lambda*t)A_alpha(e^t).
```

Then `||m||_2^2=I_alpha`, `||a||_2^2=J_alpha`, and

```text
a=m+alpha*V_lambda*m,
m=a-alpha*V_beta*a,

(V_c f)(t):=integral_0^t exp(-c(t-u))f(u)du.                (MWPADR.4)
```

The forward Fourier multiplier is

```text
(beta+it)/(lambda+it).
```

Its modulus lies between `1` and `beta/lambda`. Therefore

```text
I_alpha
 <=J_alpha
 <=((1+alpha)/(1-alpha))^2*I_alpha.                         (MWPADR.5)
```

This is an exact norm comparison, not a generic estimate for either
arithmetic function.

The Mellin transform makes the shared spectrum explicit:

```text
integral_1^infinity A_alpha(x)x^(-s-1)dx
 =1/[s*zeta(s+alpha)]                                      (MWPADR.6)
```

initially for `Re(s)>1-alpha`. If the energies are finite,
Mellin-Plancherel gives

```text
J_alpha
 =(1/(2*pi))*integral_R
   dt/[(lambda^2+t^2)|zeta(beta+it)|^2],

I_alpha
 =(1/(2*pi))*integral_R
   dt/[(beta^2+t^2)|zeta(beta+it)|^2].                     (MWPADR.7)
```

Thus both coordinates contain the same reciprocal-zeta boundary data.

## Discrete Energy

Since both summatory functions are constant on unit cells,

```text
J_alpha
 =(1/(1-alpha))*sum_(k>=1)A_alpha(k)^2
   [k^(alpha-1)-(k+1)^(alpha-1)].                          (MWPADR.8)
```

Put

```text
E_alpha:=sum_(k>=1)k^(alpha-2)A_alpha(k)^2,
M_alpha:=sum_(k>=1)M(k)^2/k^(2+alpha).
```

Unit-cell comparison and (MWPADR.5) give

```text
2^(-2-alpha)M_alpha
 <=E_alpha
 <=2^(2-alpha)*((1+alpha)/(1-alpha))^2*M_alpha.             (MWPADR.9)
```

Consequently, for any fixed cofinal sequence `alpha_j->0`,

```text
RH
 iff
E_(alpha_j)<infinity for every j.                          (MWPADR.10)
```

This is a new exact coordinate for the existing RH-equivalent energy, not
a proof that it is finite.

## First q-Block And Its Missing Mode

Set `alpha=2omega`. Then `A_alpha` is exactly the weighted Mobius prefix in
the Burnol tail-discrepancy reduction. For `N<k<=2N`,
every omitted divisor has `floor(k/d)=1`, so

```text
H_(omega,N)(k)=A_alpha(k)-A_alpha(N),

D_(omega,N)(k)
 =A_alpha(k)-A_alpha(N)-k*r_(omega,N).                     (MWPADR.11)
```

Discrete Abel summation also gives

```text
r_(omega,N)
 =-A_alpha(N)/(N+1)
  +sum_(n>N)A_alpha(n)/(n(n+1)).                           (MWPADR.12)
```

Define

```text
V_N:=sum_(N<k<=2N)|D_(omega,N)(k)|^2,
E_D(N):=sum_(N<k<=2N)k^(alpha-2)|D_(omega,N)(k)|^2,
F_alpha(N):=sum_(N<k<=2N)k^(alpha-2)|A_alpha(k)|^2.
```

Then

```text
2^(alpha-2)N^(alpha-2)V_N
 <=E_D(N)
 <=N^(alpha-2)V_N.                                        (MWPADR.13)
```

For

```text
W_q(N):=sum_(N<k<=2N)k^(alpha-2+q),
```

the missing affine line has the exact Gram

```text
||A_alpha(N)+k*r_N||_w^2
 =A_alpha(N)^2*W_0
  +2A_alpha(N)r_N*W_1
  +r_N^2*W_2,                                              (MWPADR.14)
```

where `W_0*W_2-W_1^2>0`. In particular,

```text
F_alpha(N)^(1/2)
 <=N^((alpha-2)/2)V_N^(1/2)
   +N^((alpha-1)/2)|A_alpha(N)|
   +2^(alpha/2)N^((1+alpha)/2)|r_N|.                       (MWPADR.15)
```

The first q-block is invariant under adding the same constant to
`A_alpha(N)` and every later prefix in the block. This is a genuine
unobserved anchor mode, not a notational artifact.

For a generic finitely supported coefficient sequence, take its prefix to
be a nonzero constant `C` from `N` onward. Then

```text
r_N=0,
D_N(k)=0 for every k>N,
F_alpha(N)=C^2*W_0(N)>0.                                   (MWPADR.16)
```

The exact tail identity (MWPADR.12) still holds. Hence no operator
inequality can bound the absolute prefix energy from `V_N` and the forced
reciprocal-tail rate alone. A separate arithmetic theorem must control the
constant anchor `N^((alpha-1)/2)A_alpha(N)`.

## Anchored Correlation Collapses, But Does Not Decouple

The signed two-dimensional correlation from the preceding reduction has
the exact one-dimensional form

```text
sum_(h=1)^(X-1)sum_(Y=1)^(X-h)C_h(Y)
 =sum_(n=2)^X (X-n+1)mu(n)M(n-1).                          (MWPADR.17)
```

This follows by using the larger index `n=m+h`. It is the discrete energy
increment identity

```text
M(n)^2-M(n-1)^2
 =2mu(n)M(n-1)+mu(n)^2.                                    (MWPADR.18)
```

Thus the collapse removes one summation variable but leaves the unknown
Mertens prefix inside the test function.

## Vaughan Audit

Green and Tao's Mobius Vaughan identity gives, for `1<=U,V<=N`,

```text
sum_(N<n<=2N)mu(n)f(n)
 =-sum_(d<=U*V)a_d
    sum_(N/d<w<=2N/d)f(d*w)

  +sum_(V<d<=2N/U)b_d
    sum_(max(U,N/d)<w<=2N/d)mu(w)f(d*w),                   (MWPADR.19)

a_d:=sum_(b*c=d, b<=U, c<=V)mu(b)mu(c),
b_d:=sum_(c|d, c>V)mu(c).
```

The identity is finite and valid for an arbitrary test vector `f`.
Applied after (MWPADR.17), however, the test is

```text
f_X(n):=(X-n+1)M(n-1).
```

It is endogenous, and

```text
sum_(n<=X)|f_X(n)|^2
 <=X^2*sum_(k<=X)M(k)^2.                                   (MWPADR.20)
```

Therefore the standard Cauchy-Schwarz elimination of the Type I/II
coefficients feeds the target energy back into the estimate. Vaughan's
identity alone supplies no absorbable small factor and does not prove RH.

There is one noncircular handoff worth retaining. Before collapsing the
shift sums, apply (MWPADR.19) to each dyadic `m`-block with

```text
f_(h,Y)(m)=mu(m+h)1_(m<=Y),
```

which is externally bounded by one. Let `N` run over powers of two below
`X`, choose `1<=U_N,V_N<=N`, and define

```text
a_(N,d):=sum_(b*c=d,b<=U_N,c<=V_N)mu(b)mu(c),
b_(N,d):=sum_(c|d,c>V_N)mu(c),

TI_X
 :=sum_N sum_(h<X) sum_(Y<=X-h)
   sum_(d<=U_N*V_N) a_(N,d)
   sum_(N/d<w<=min(2N,Y)/d) mu(d*w+h),

TII_X
 :=sum_N sum_(h<X) sum_(Y<=X-h)
   sum_(V_N<d<=2N/U_N) b_(N,d)
   sum_(max(U_N,N/d)<w<=min(2N,Y)/d)
      mu(w)mu(d*w+h).                                      (MWPADR.21)
```

Empty sums vanish. The dyadic blocks `(N,2N]` cover `2<=m<=X`; the
omitted `m=1` boundary is

```text
B_1(X):=sum_(n=2)^X mu(n)(X-n+1)=O(X^2).
```

Termwise use of (MWPADR.19) gives the exact aggregate identity

```text
sum_(h<X)sum_(Y<=X-h)C_h(Y)
 =B_1(X)-TI_X+TII_X.                                      (MWPADR.22)
```

The remaining obligation is therefore the concrete signed estimate

```text
-TI_X+TII_X
 =O_epsilon(X^(2+epsilon)).                                (MWPADR.23)
```

The summation must retain cancellation jointly in the shift and terminal
variables. Rowwise absolute values restore the previously identified cubic
loss.

## Open Handoff

Two gates remain:

```text
1. Control N^((alpha-1)/2)A_alpha(N) in a summable dyadic
   theorem, using the limiting Burnol energy or another arithmetic input.

2. Prove (MWPADR.23) for the pre-collapse Vaughan Type I/II
   aggregate without assuming RH.
```

Neither gate is proved here. The reduction proves no full Burnol bound,
RH, PF-infinity, or `Lambda <= 0`.

## Source Boundary

- [Green and Tao, Quadratic uniformity of the Mobius
  function](https://doi.org/10.5802/aif.2401), Lemma 4.1, supplies the
  finite Vaughan Type I/II identity. Their bounded-test estimates are not
  attributed to the endogenous Mertens test.
- [Burnol](https://arxiv.org/abs/math/0202166) supplies the published
  weighted natural-approximant setting. The first-block and anchor-defect
  reductions above are derived here.
- [Das and Manna](https://arxiv.org/abs/2508.00388) provide modern
  Hardy/Copson context. The sharp scalar multiplier comparison
  (MWPADR.5) is derived directly here.
