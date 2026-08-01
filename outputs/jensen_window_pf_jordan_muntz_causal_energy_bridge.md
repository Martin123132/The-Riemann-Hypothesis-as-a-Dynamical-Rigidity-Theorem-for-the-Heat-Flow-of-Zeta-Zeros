# Jensen-Window PF Jordan-Muntz Causal-Energy Bridge

Date: 2026-07-23

Status: exact Jordan-Muntz and finite-energy reductions with internally
audited fixed-shift/cofinal theorem candidates and one open natural-mollifier
gate. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_jordan_muntz_causal_energy_bridge.json
python work/rh_compute/scripts/jensen_window_pf_jordan_muntz_causal_energy_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py
```

## Generalized Muntz Identity

For `0<omega<1/2`, retain

```text
c_omega(n)
 =sum_(d*m=n) mu(d)d^(-omega)m^omega,

C_omega(x)=sum_(n<=x)c_omega(n),

A_omega=1/((1+omega)*zeta(1+2omega)),

E_omega(x)=C_omega(x)-A_omega*x^(1+omega).
```

Define the fractional-power remainder

```text
R_omega(y)
 =sum_(m<=y)m^omega-y^(1+omega)/(1+omega).           (JMCE.1)
```

For `0<y<1`, this is exactly
`-y^(1+omega)/(1+omega)`. Expanding `C_omega` by its Dirichlet
convolution and using

```text
sum_(d>=1)mu(d)d^(-(1+2omega))=1/zeta(1+2omega)
```

gives the full identity

```text
E_omega(x)
 =sum_(d>=1)mu(d)d^(-omega)R_omega(x/d).             (JMCE.2)
```

This version includes the earlier `d>x` tail rather than writing it
separately. That tail is absolutely convergent because its terms are

```text
-mu(d)*x^(1+omega)/((1+omega)*d^(1+2omega)).
```

The elementary power-sum estimate gives
`R_omega(y)=O_omega(y^omega)` at infinity. Together with its exact
behavior below one, this proves

```text
R_omega in H,
H=L2((0,infinity),dx/x^2),                           (JMCE.3)
```

precisely for `0<omega<1/2`. Its Mellin transform in the honest
convergence strip is

```text
integral_0^infinity R_omega(x)x^(-s)dx/x
 =zeta(s-omega)/s,
omega<Re(s)<1+omega.                                 (JMCE.4)
```

## Natural Mobius Dilations

The operators

```text
(U_d f)(x)=sqrt(d)f(x/d)
```

are unitary on `H`. Hence the finite natural approximants are

```text
E_(omega,N)(x)
 =sum_(d<=N)mu(d)d^(-omega)R_omega(x/d)
 =sum_(d<=N)mu(d)d^(-omega-1/2)U_dR_omega(x).        (JMCE.5)
```

Put

```text
M_N(u)=sum_(d<=N)mu(d)d^(-u).
```

Mellin-Plancherel gives the exact finite identity

```text
||E_(omega,N)||_H^2
 =1/(2*pi)*integral_R
  |zeta(1/2-omega+i*t)/(1/2+i*t)|^2
  |M_N(1/2+omega+i*t)|^2 dt.                         (JMCE.6)
```

No infinite series has been moved through this integral. The zeta functional
equation rewrites the same integrand as

```text
W_omega(t)
*|zeta(sigma+i*t)M_N(sigma+i*t)|^2,

sigma=1/2+omega,

W_omega(t)
 =|chi(1/2-omega+i*t)|^2/|1/2+i*t|^2,               (JMCE.7)
```

where `zeta(v)=chi(v)zeta(1-v)`. Thus the live arithmetic object is a
weighted all-height norm of the ordinary natural mollifier on a line
strictly to the right of one half.

For every fixed `x`, (JMCE.2) gives

```text
E_(omega,N)(x)->E_omega(x).
```

Consequently,

```text
sup_N ||E_(omega,N)||_H < infinity
=> E_omega in H                                      (JMCE.8)
```

by Fatou. Square summability of the scalar coefficients
`d^(-omega-1/2)` is not enough: the dilation orbit has not been proved to
be a Bessel sequence.

## Causal Error Criterion

The part below `x=1` is explicit:

```text
integral_0^1 |E_omega(x)|^2 dx/x^2
 =A_omega^2/(1+2omega).
```

Therefore define the only open part by

```text
J_omega
 =integral_1^infinity |E_omega(x)|^2 dx/x^2.         (JMCE.9)
```

With `s=z+1/2`, its positive-log-time Laplace transform is

```text
G_omega(s)
 =zeta(s-omega)/(s*zeta(s+omega))
  -A_omega/(s-1-omega).                              (JMCE.10)
```

The subtraction removes the residue at `s=1+omega`. Initially this is an
ordinary transform for `Re(s)>1-omega`, using the elementary error bound.

If `J_omega<infinity`, Paley-Wiener puts (JMCE.10) in `H2` of
`Re(s)>1/2`. Every possible pole from a zero of `zeta(s+omega)` must then
be removable, so Suzuki's reduced fixed-shift condition `D(omega)` holds.

Conversely, if `D(omega)` holds, the completed quotient

```text
Q_omega(s-1/2)
 =xi(s-omega)/xi(s+omega)
```

is inner. Factoring (JMCE.10) through `Q_omega`, removing its single
archimedean pole at `1+omega`, and using the standard gamma quotient
estimate `O(|t|^(omega-1))` gives `G_omega in H2`. Transform uniqueness
then recovers the actual arithmetic error. The resulting fixed-shift
theorem candidate is

```text
J_omega<infinity
iff D(omega).                                        (JMCE.11)
```

Combining this with the audited cofinal phase theorem gives

```text
RH
iff there are omega_j->0 such that
    integral_1^infinity |E_(omega_j)(x)|^2 dx/x^2
    is finite for every j.                           (JMCE.12)
```

Burnol's shifted-ratio estimate and the Balazard-Saias natural-Mobius
approximation show under RH that `E_(omega,N)->E_omega` in `H` for every
fixed `omega`. In the reverse direction, (JMCE.8) and (JMCE.12) require
only a uniform norm bound. This yields the finite natural criterion

```text
RH
iff there are omega_j->0 such that
    sup_N ||E_(omega_j,N)||_H < infinity
    for every j.                                     (JMCE.13)
```

Equations (JMCE.11)-(JMCE.13) are new corpus theorem candidates requiring
independent review. They do not assert the required bound.

## Finite Brownian Energy

For every finite `X>=1`, direct expansion inside the finite integral gives

```text
J_omega(X)
 =integral_1^X E_omega(x)^2 dx/x^2

 =sum_(n,m<=X)c_omega(n)c_omega(m)
   *(1/max(n,m)-1/X)

  -(2*A_omega/omega)
   *sum_(n<=X)c_omega(n)*(X^omega-n^omega)

  +A_omega^2*(X^(1+2omega)-1)/(1+2omega).            (JMCE.14)
```

There is also a clean signed-measure form. Put

```text
d eta_(omega,X)(u)
 =sum_(n<=X)c_omega(n)delta_n
  -A_omega*(1+omega)u^omega*1_(0,X](u)du.
```

Then

```text
J_omega(X)
 =double_integral_[0,X]^2
  K_X(u,v)deta(u)deta(v),

K_X(u,v)
 =integral_1^X 1_(u<=x)1_(v<=x)dx/x^2
 =1/max(1,u,v)-1/X.                                  (JMCE.15)
```

The kernel is positive semidefinite. Under reciprocal coordinates it is
a Brownian `min` kernel, and every ordered finite matrix is totally
nonnegative by a cumulative-matrix factorization. This is exact contact
with total positivity, but it supplies only `J_omega(X)>=0`. The
positive-positive, cross, and continuum-continuum terms in (JMCE.14) each
have the same potentially divergent scale; a uniform upper bound requires
their arithmetic cancellation.

## Countermodel Gates

Boundary energy is automatic. For `a>0`,

```text
1/(a+z),  1/(a-z)
```

have the same boundary squared norm `1/(2a)`. The first is the Laplace
transform of the causal `exp(-a*t)1_(t>0)`. The second has anti-causal
boundary inverse `exp(a*t)1_(t<0)`, while its actual positive-time
Bromwich inverse is `-exp(a*t)` and is not square integrable. Thus
Plancherel without Hardy support cannot prove (JMCE.11).

Likewise, scalar coefficients in `l2` do not make an arbitrary unit-vector
series convergent: `v_d=v` and coefficients `1/d` are the elementary
countermodel. A Bessel or cancellation estimate for the actual dilation
orbit is a theorem obligation, not a Hilbert-space formality.

## Literature Fit

Burnol's `f_epsilon` theorem has the same shifted-zeta-ratio,
square-integrability, natural-Mobius, and causality architecture. The
present `R_omega` is the fractional-power Jordan remainder naturally
selected by Suzuki's coefficients.

Baez-Duarte's general strong Nyman-Beurling theorem does not close
(JMCE.13): its audited sufficiency direction requires compact support,
whereas `R_omega` is not compactly supported, and arbitrary closure does
not prove convergence of this fixed natural Mobius sequence.

Primary sources:

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- Jean-Francois Burnol, `On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte`: https://arxiv.org/abs/math/0202166
- Luis Baez-Duarte, `A general strong Nyman-Beurling Criterion for the Riemann Hypothesis`: https://arxiv.org/abs/math/0505453

## Open Gate

For one explicit cofinal sequence, prove the all-height finite inequality

```text
sup_N integral_R W_omega(t)
 |zeta(1/2+omega+i*t)M_N(1/2+omega+i*t)|^2 dt
 <infinity.                                          (JMCE.16)
```

This is a sharper target than another finite `x` or finite `t` grid. It is
also exactly RH-strength through (JMCE.13), so no claim of an easy
mean-square estimate is made.

## Proof Boundary

Equations (JMCE.1)-(JMCE.10), (JMCE.14), and (JMCE.15) are exact
identities, convergence statements, or finite Hilbert-space formulas.
Equations (JMCE.11)-(JMCE.13) are internally audited theorem candidates
using standard Hardy-space and conditional natural-approximation
machinery. The uniform estimate (JMCE.16) is open. This artifact does not
prove that estimate, RH, or `Lambda <= 0`.
