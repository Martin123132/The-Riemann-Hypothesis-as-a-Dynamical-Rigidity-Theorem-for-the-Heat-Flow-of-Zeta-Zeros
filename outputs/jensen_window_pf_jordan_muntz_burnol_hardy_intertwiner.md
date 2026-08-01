# Jensen-Window PF Jordan-Muntz/Burnol Hardy Intertwiner

Date: 2026-07-23

Status: exact finite operator bridge with one source-backed cofinal theorem
candidate and one open natural-mollifier gate. This is not a proof of RH or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json
python work/rh_compute/scripts/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py
```

## Fractional-Part Intertwiner

Put

```text
R_0(y)=floor(y)-y=-{y}.
```

Stieltjes integration by parts gives, for `0<omega<1`,

```text
R_omega(y)
 =-y^omega*{y}
  +omega*integral_0^y {u}u^(omega-1)du.             (JMBH.1)
```

Define

```text
(T_omega f)(x)
 =x^omega*f(x)
  -omega*integral_0^x f(u)u^(omega-1)du.            (JMBH.2)
```

Then (JMBH.1) is exactly

```text
R_omega=T_omega R_0.                                (JMBH.3)
```

For `D_d f(x)=f(x/d)`, a change of variables gives

```text
T_omega D_d=d^omega D_d T_omega.                    (JMBH.4)
```

Consequently the finite Jordan-Muntz partial sum factors as

```text
E_(omega,N)(x)
 =T_omega[
   sum_(d<=N)mu(d)d^(-2omega)R_0(x/d)].             (JMBH.5)
```

No infinite sum or limiting argument occurs in (JMBH.5).

## Burnol Coordinates

Use Burnol's finite natural fractional-part approximant

```text
f_(epsilon,N)(t)
 =sum_(d<=N)mu(d)d^(-epsilon){1/(d*t)}.             (JMBH.6)
```

The reciprocal substitution `x=1/t` is unitary:

```text
integral_0^infinity |F(x)|^2 dx/x^2
 =integral_0^infinity |F(1/t)|^2 dt.                (JMBH.7)
```

Let

```text
(H_star g)(t)=integral_t^infinity g(v)dv/v.         (JMBH.8)
```

Applying (JMBH.5) at `x=1/t`, with `epsilon=2omega`, gives the exact finite
intertwiner

```text
E_(omega,N)(1/t)
 =-(I-omega*H_star)
   [t^(-omega)f_(2omega,N)(t)].                     (JMBH.9)
```

This identifies the Jordan partial sums with the weighted functions used in
Burnol's discussion of Baez-Duarte's natural approximants.

## Two-Sided Norm

The classical tail Hardy operator satisfies

```text
||H_star||_(L2(dt)->L2(dt))=2.
```

One proof sets `G=H_star g`, uses `G'=-g/t`, and integrates
`integral G^2` by parts before Cauchy-Schwarz. Therefore, for
`0<omega<1/2`,

```text
(1-2omega)||t^(-omega)f_(2omega,N)||_2
 <=||E_(omega,N)||_H
 <=(1+2omega)||t^(-omega)f_(2omega,N)||_2.          (JMBH.10)
```

Moreover `I-omega H_star` is boundedly invertible by its Neumann series.
Thus uniform boundedness, Cauchy convergence, and convergence of either
finite family are equivalent to the corresponding property of the other.
The constants improve, rather than deteriorate, along `omega->0`.

## Mellin Check

In a common convergence strip,

```text
M[H_star g](s)=M[g](s)/s,

M[t^(-omega)f_(2omega,N)](s)
 =-zeta(s-omega)M_N(s+omega)/(s-omega).             (JMBH.11)
```

The multiplier of `I-omega H_star` is `(s-omega)/s`. Combining this with
the minus sign in (JMBH.9) recovers

```text
M[E_(omega,N)(1/t)](s)
 =zeta(s-omega)M_N(s+omega)/s,                      (JMBH.12)
```

which is the independently obtained Jordan-Muntz transform.

## Cofinal Consequence

Suppose `omega_j->0` and

```text
sup_N ||E_(omega_j,N)||_H<infinity
```

for every `j`. Equation (JMBH.10) uniformly bounds
`t^(-omega_j)f_(2omega_j,N)`. The infinite fractional-part series converges
pointwise for each positive epsilon, so Fatou gives the weighted Burnol
function in `L2`. On `(0,1)`, this implies the unweighted function is in
`L2`; on `(1,infinity)` it is an explicit constant times `1/t`. Burnol's
published sequence theorem then implies RH.

Conversely, under RH the Balazard-Saias estimate used by Burnol gives `L2`
convergence of these weighted natural partial sums. The bounded operator in
(JMBH.9) transfers convergence to `E_(omega,N)`. This yields the
source-backed cofinal theorem candidate

```text
RH
iff there are omega_j->0 such that
    sup_N ||E_(omega_j,N)||_H<infinity
    for every j.                                    (JMBH.13)
```

The normalization transfer and published-theorem composition in (JMBH.13)
still warrant independent expert review, but the route no longer depends
only on the new archimedean factorization.

## What Remains Open

The operator identity does not provide the arithmetic estimate. By
(JMBH.10), the live target is equivalently

```text
sup_N ||t^(-omega)f_(2omega,N)||_2<infinity          (JMBH.14)
```

on one explicit cofinal sequence, without assuming RH. Burnol explicitly
locates these natural Mobius approximants inside the Nyman programme; the
Balazard-Saias convergence estimate used there is conditional on RH.

Baez-Duarte's general strong Nyman theorem does not bypass (JMBH.14): the
audited sufficiency direction requires compact support, and `R_omega` is not
compactly supported. Hardy boundedness also cannot be mistaken for a bound
on the arithmetic input.

Primary sources:

- Jean-Francois Burnol, `On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte`: https://arxiv.org/abs/math/0202166
- Luis Baez-Duarte, `A general strong Nyman-Beurling Criterion for the Riemann Hypothesis`: https://arxiv.org/abs/math/0505453
- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823

## Proof Boundary

Equations (JMBH.1)-(JMBH.12) are finite identities, unitary substitutions,
or consequences of the classical Hardy inequality. Equation (JMBH.13) is a
source-backed theorem candidate whose normalization chain requires external
review. Equation (JMBH.14) is open. RH, PF-infinity, and `Lambda <= 0`
remain unproved.
