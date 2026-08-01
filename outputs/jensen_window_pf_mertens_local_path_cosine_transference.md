# Jensen-Window PF Mertens Local-Path Cosine Transference

Date: 2026-07-23

Status: exact local Mertens-path DCT transference with an isolated
generic first-mode barrier and one open two-component gate. This is not a proof
of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_local_path_cosine_transference.json
python work/rh_compute/scripts/jensen_window_pf_mertens_local_path_cosine_transference.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_local_path_cosine_transference.py
```

## Exact Local-Path Transform

Fix `0<alpha<1` and retain

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.                    (MLPATH.1)
```

For dyadic `K`, set

```text
y_j:=r_(K+j-1),                                      1<=j<=K,

T_(K,m):=sum_(ell=1)^m q_(K+ell)
        =r_K-r_(K+m),                                0<=m<=K,

z_j:=T_(K,j-1).                                      (MLPATH.2)
```

Thus `T_(K,0)=0` and, exactly,

```text
y_j=r_K-z_j.                                         (MLPATH.3)
```

Let `phi_r` be the orthonormal Neumann DCT basis from Lemma 11.22V:

```text
phi_0(j)=K^(-1/2),

phi_r(j)=sqrt(2/K)cos(pi*r*(j-1/2)/K),      1<=r<K.
```

Write

```text
c_(K,r):=<y,phi_r>,
tau_(K,r):=<z,phi_r>.
```

Because every nonconstant DCT vector has zero sum, (MLPATH.3) gives the
affine transference

```text
c_(K,0)=sqrt(K)r_K-tau_(K,0),
c_(K,r)=-tau_(K,r),                         1<=r<K.  (MLPATH.4)
```

Discrete summation by parts gives the equivalent nonconstant formula

```text
c_(K,r)
 =-sqrt(2/K) sum_(m=1)^(K-1)
    T_(K,m)cos(pi*r*(m+1/2)/K),             1<=r<K.  (MLPATH.5)
```

The factor `1/sin(pi*r/(2K))` in the sine-weight formula from
Lemma 11.22V has canceled exactly. The low modes are therefore the
ordinary cosine coefficients of the reciprocal-weighted local Mertens
path, not generic high-denominator oscillatory sums.

Parseval now reads

```text
||y||_2^2
 =|sqrt(K)r_K-tau_(K,0)|^2
  +sum_(r=1)^(K-1)|tau_(K,r)|^2,                     (MLPATH.6)

sum_(r=1)^(K-1)|tau_(K,r)|^2
 =sum_(m=0)^(K-1)|T_(K,m)|^2
  -K^(-1)|sum_(m=0)^(K-1)T_(K,m)|^2.                 (MLPATH.7)
```

So the nonconstant spectrum is exactly the centered path variance. If
the constant DCT row and the fixed column `T_(K,0)=0` are deleted, the
remaining square transform `U_K` satisfies

```text
U_K^*U_K=I_(K-1)-K^(-1)11^T.                         (MLPATH.8)
```

Its singular values are `1` with multiplicity `K-2` and `K^(-1/2)`
with multiplicity one. The weak direction is the near-constant local
path. It cannot be separated from the affine mean defect in (MLPATH.4)
without losing the exact cancellation.

## RH-Equivalent Two-Component Gate

Let `R_K=ceil(sqrt(K))`. Lemma 11.22V already proves that every mode
`r>=R_K` is dyadically summable from coefficient size alone. Equations
(MLPATH.4)-(MLPATH.6) sharpen its criterion to

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^alpha
 (
   |sqrt(K)r_K-tau_(K,0)|^2
   +sum_(1<=r<R_K)|tau_(K,r)|^2
 )
<infinity.                                            (MLPATH.9)
```

For every member of any one fixed cofinal sequence `alpha_j->0`,
(MLPATH.9) is equivalent to RH by the checked reciprocal-tail
criterion. It has two genuinely different pieces:

```text
affine mean defect:
  sqrt(K)r_K-tau_(K,0);

nonconstant local-path spectrum:
  tau_(K,r),                         1<=r<ceil(sqrt(K)).
```

The first contains the incoming global tail and local path mean. The
second contains the smooth shape of the anchored path. Controlling one
does not algebraically control the other.

## Transfer To Ordinary Local Mertens Increments

Define

```text
S_(K,m):=M(K+m)-M(K)
        =sum_(ell=1)^m mu(K+ell),

a_m:=(K+m)^(-1-alpha),
b_m:=1/a_m.                                           (MLPATH.10)
```

Finite Abel summation and its inverse give

```text
T_(K,m)
 =a_m S_(K,m)
  +sum_(ell=1)^(m-1)
    S_(K,ell)(a_ell-a_(ell+1)),                       (MLPATH.11)

S_(K,m)
 =b_m T_(K,m)
  +sum_(ell=1)^(m-1)
    T_(K,ell)(b_ell-b_(ell+1)).                       (MLPATH.12)
```

For the vectors indexed by `1<=m<K`, the strict-prefix matrix has
operator norm at most `K`. Since

```text
a_m<=K^(-1-alpha),
a_ell-a_(ell+1)<=(1+alpha)K^(-2-alpha),

b_m<=(2K)^(1+alpha),
b_(ell+1)-b_ell<=(1+alpha)(2K)^alpha,
```

(MLPATH.11)-(MLPATH.12) imply the explicit full-path comparison

```text
[2^alpha(3+alpha)]^(-1)
 K^(-1-alpha)||S_K||_2
 <=||T_K||_2
 <=(2+alpha)K^(-1-alpha)||S_K||_2.                   (MLPATH.13)
```

This does not compare their first `sqrt(K)` DCT projections. The Abel
matrix is triangular and does not commute with the low-mode projector.
Using (MLPATH.13) after a generic Cauchy-Schwarz step therefore recycles
the full Mertens-energy target rather than proving the missing
low-mode gain.

## Zero-Mean First-Mode Countermodel

The following exact scalar construction proves that coefficient size,
absolute convergence, endpoint-tail cancellation, and mean cancellation
cannot close the nonconstant modes generically.

For this fixed `alpha` and every dyadic `K>=4`, put

```text
A_K:=2^(-2-alpha)pi^(-1)K^(-alpha),

T*_(K,m):=A_K sin(2pi*m/K),                 0<=m<=K,

q*_(K+m):=T*_(K,m)-T*_(K,m-1),             1<=m<=K.
                                                            (MLPATH.14)
```

Set `q*_1=0`. The annuli `(K,2K]` partition all remaining indices. For
`n=K+m`,

```text
|q*_n|
 <=2pi*A_K/K
 =2^(-1-alpha)K^(-1-alpha)
 <=n^(-1-alpha).                                    (MLPATH.15)
```

The block `l1` norm is at most `2pi*A_K`, so the global sequence is
absolutely summable. Moreover,

```text
sum_(K<n<=2K)q*_n=T*_(K,K)-T*_(K,0)=0.
```

All later block totals also vanish. Hence, at every dyadic endpoint,

```text
r*_K=0,
tau*_(K,0)=0,
c*_(K,0)=0.                                          (MLPATH.16)
```

Nevertheless the first nonconstant coefficient is

```text
|c*_(K,1)|
 =A_K sqrt(2/K) cos(pi/K)/2
  [csc(pi/(2K))+csc(3pi/(2K))].                      (MLPATH.17)
```

Indeed, product-to-sum gives

```text
sum_(m=0)^(K-1)
 sin(2pi*m/K)cos(pi*(m+1/2)/K)

 =cos(pi/K)/2
  [csc(pi/(2K))+csc(3pi/(2K))].
```

For `K>=4`, `cos(pi/K)>=sqrt(2)/2` and `sin x<=x`, so

```text
|c*_(K,1)|>=A_K sqrt(K)/pi,

K^alpha|c*_(K,1)|^2
 >=2^(-4-2alpha)pi^(-4)K^(1-alpha).                  (MLPATH.18)
```

The dyadic series diverges although every dyadic tail endpoint and
every mean mode vanish. This is deliberately not a Mobius model. It
proves a method barrier: the actual proof must use arithmetic signed
cancellation in the nonconstant local path, not only its coefficient
envelope or the scalar tail level.

## Source And Method Audit

Davenport's theorem gives, uniformly in `theta`,

```text
sum_(n<=X)mu(n)e(n theta)<<_A X(log X)^(-A).
```

That arbitrary-log estimate is power-short for (MLPATH.9). The focused
primary-source audit also checked:

- R. C. Baker and G. Harman, *Exponential Sums Formed with the Mobius
  Function*:
  https://doi.org/10.1112/jlms/s2-43.2.193
- H. Maier and A. Sankaranarayanan, *On an Exponential Sum Involving
  the Mobius Function*:
  https://doi.org/10.46298/hrj.2005.151
- W. Zhang, *On an Exponential Sum Related to the Mobius Function*:
  https://arxiv.org/abs/2204.04613

Zhang records the unconditional Davenport estimate and quotes
Baker-Harman power estimates under the hypothesis that every Dirichlet
`L(s,chi)` is zero-free in `Re(s)>a`. At the endpoint `a=1/2`, the
quoted uniform exponent is `3/4+epsilon`. Lemma 11.22V requires

```text
sigma<(1+alpha)/2,
```

so `sigma=3/4+epsilon` closes only
`alpha>1/2+2epsilon`, never a cofinal `alpha->0` sequence. The
unconditional rational estimate of Maier-Sankaranarayanan, as quoted
in Zhang's Lemma 2.1, retains Dirichlet-L zero-density terms. The audit
found no source-backed simplification yielding the required uniform
cofinal-closing power. This is a focused search report, not a novelty,
priority, or exhaustive-literature claim.

A viable Vaughan or bilinear attack must therefore enter before
rowwise absolute values, preserve signed coupling over the local modes
and dyadic scales, and treat the affine mean defect at the same time.
Translation-averaged and almost-all interval theorems need an explicit
de-averaging step. Estimates only for nonzero phases do not treat the
mean defect. Finally, (MLPATH.9) concerns `R_alpha`, not the complete
Burnol energy `Q`.

The open gate is:

```text
for every alpha in one fixed cofinal sequence alpha_j->0,

sum_(K dyadic)K^alpha
 (
   |sqrt(K)r_K-tau_(K,0)|^2
   +sum_(1<=r<ceil(sqrt(K)))|tau_(K,r)|^2
 )
<infinity.                                            (MLPATH.19)
```

This is RH-equivalent. The exact transference and countermodel sharpen
where new arithmetic input must act, but supply no such gain. RH,
PF-infinity, and `Lambda <= 0` remain open.
