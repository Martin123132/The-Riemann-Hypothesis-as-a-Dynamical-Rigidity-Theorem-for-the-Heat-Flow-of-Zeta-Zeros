# Formal Core

Date: 2026-07-03

Status: exact-lemma ledger. This is not a proof of RH or `Lambda <= 0`; it records noncircular identities and separates them from conjectural bridges.

Project:

```text
RH Dynamical Rigidity / de Bruijn-Newman Zero-Flow Programme
```

Purpose:

Collect the exact, non-circular mathematical statements that can safely be used in later proof development. Numerical observations, conjectural bridges, and RH-equivalent assumptions are deliberately separated from exact lemmas.

## 0. Conventions

We use the standard de Bruijn-Newman family

```text
H_lambda(z) = integral_0^infty exp(lambda u^2) Phi(u) cos(z u) du.
```

It satisfies the heat-type equation

```text
partial_lambda H_lambda(z) = - partial_zz H_lambda(z).
```

The Newman parameter is `lambda`.

The de Bruijn-Newman constant `Lambda` is defined so that `H_lambda` has only real zeros in the real-zero regime. With the usual convention:

```text
RH is equivalent to Lambda <= 0.
```

Rodgers-Tao proves:

```text
Lambda >= 0.
```

Therefore this programme must prove:

```text
Lambda <= 0
```

without assuming RH at `lambda = 0`.

## 1. Zero Tracking

### Lemma 1.1: Simple Zero Tracking

Let `F(lambda,z)` be analytic in both variables and real-valued on the real axis for real `lambda`. Suppose

```text
F(lambda_0, x_0) = 0,
partial_z F(lambda_0, x_0) != 0.
```

Then there is a unique differentiable real branch `x(lambda)` near `lambda_0` such that

```text
x(lambda_0) = x_0,
F(lambda, x(lambda)) = 0.
```

The branch satisfies

```text
x'(lambda) = - F_lambda(lambda,x(lambda)) / F_z(lambda,x(lambda)).
```

Proof:

Immediate from the implicit function theorem and differentiation of

```text
F(lambda, x(lambda)) = 0.
```

### Corollary 1.2: Newman Zero Velocity

For `F = H_lambda`, since

```text
F_lambda = -F_zz,
```

any simple real zero branch satisfies

```text
x_i'(lambda) = F_zz(lambda,x_i) / F_z(lambda,x_i).
```

This is exact wherever the zero remains simple.

## 2. Local Product Expansion And Singular Velocity

### Lemma 2.1: Local Hadamard/Weierstrass Velocity Form

Suppose that near a finite cluster of simple real zeros

```text
x_1(lambda), ..., x_r(lambda)
```

the function can be written as

```text
F(lambda,z) = U(lambda,z) product_{k=1}^r (z - x_k(lambda)),
```

where `U` is analytic and nonzero in the neighbourhood.

Then at a zero `x_i`,

```text
F_zz/F_z (x_i)
= 2 sum_{1 <= j <= r, j != i} 1/(x_i - x_j)
  + 2 partial_z log U(x_i).
```

For the Newman flow,

```text
x_i'
= 2 sum_{j != i} 1/(x_i - x_j)
  + E_i,
```

where

```text
E_i = 2 partial_z log U(x_i)
```

is the local external field contributed by all structure outside the finite cluster.

Proof sketch:

Write `P(z) = product_k (z - x_k)`. At a simple zero `x_i`,

```text
P_zz/P_z (x_i) = 2 sum_{j != i} 1/(x_i - x_j).
```

The logarithmic derivative of the nonvanishing analytic factor contributes `2 U_z/U`.

## 3. Adjacent Gap Equation In The Finite Pairwise Model

Consider the reduced finite model

```text
x_i' = 2 sum_{j != i} 1/(x_i - x_j),
```

with ordered real points

```text
x_1 < x_2 < ... < x_n.
```

Let

```text
g_i = x_{i+1} - x_i.
```

Then

```text
g_i'
= 4/g_i
  - 2 g_i sum_{j != i,i+1}
      1 / ((x_{i+1}-x_j)(x_i-x_j)).
```

For every outside point `x_j`, the product

```text
(x_{i+1}-x_j)(x_i-x_j)
```

is positive. Therefore the pairwise tail in this reduced model is compressive:

```text
tail <= 0.
```

Consequence:

A bridge lemma asserting nonnegative pairwise tail is false for the reduced finite model. Any useful gap comparison theorem must use additional Xi-specific structure, a controlled-negativity estimate, or a different global principle.

## 4. Finite Cluster Non-Collapse In The Repulsive Direction

### Lemma 4.1: Internal Variance Growth

For the internal reduced flow on `r` points

```text
x_i' = 2 sum_{j != i} 1/(x_i - x_j),
```

define the cluster mean

```text
m = (1/r) sum_i x_i
```

and internal variance

```text
V = sum_i (x_i - m)^2.
```

Then the internal pairwise contribution satisfies

```text
V'_internal = 2r(r-1).
```

Proof:

Since `sum_i (x_i-m)=0`,

```text
V' = 2 sum_i (x_i-m) x_i'.
```

Pairing terms `i,j` gives one unit contribution per ordered interaction pair, hence

```text
V'_internal
= 4 sum_{i<j} 1
= 2r(r-1).
```

### Lemma 4.2: Non-Collapse With Locally Lipschitz External Field

Suppose a finite cluster evolves by

```text
x_i'
= 2 sum_{j != i} 1/(x_i - x_j)
  + E(x_i,lambda),
```

and suppose `E` is locally Lipschitz in `x` across the cluster, with Lipschitz constant `L`.

Then

```text
V' >= 2r(r-1) - 2LV.
```

Hence a positive-variance cluster cannot collapse to `V=0` in finite forward `lambda` time in the repulsive direction.

Interpretation:

This is a genuine local theorem about already-real finite clusters. It does not exclude birth of a real pair from a complex conjugate pair at a positive Newman boundary.

## 5. Square-Root Boundary Normal Form

### Lemma 5.1: Local Double-Zero Normal Form

Let `F(lambda,z)` be analytic and real on the real axis. Suppose

```text
F(lambda_*, t_*) = 0,
F_z(lambda_*, t_*) = 0,
F_zz(lambda_*, t_*) != 0,
```

and the parameter is transverse:

```text
F_lambda(lambda_*, t_*) != 0.
```

Then, locally,

```text
F(lambda,z)
= U(lambda,z) [ (z - t_* - r(lambda))^2 - s(lambda) ],
```

where `U` is analytic and nonzero, and `s(lambda_*)=0` with `s'(lambda_*) != 0`.

The local zeros are

```text
z_+(lambda), z_-(lambda)
= t_* + r(lambda) +/- sqrt(s(lambda)).
```

Thus the transition between a conjugate pair and two real zeros is square-root type.

## 6. Minimal Obstruction To The Local-Repulsion Proof Strategy

Consider

```text
F(lambda,t) = (t-a)^2 - 2(lambda-lambda_*).
```

Then

```text
F_lambda = -F_tt.
```

The zeros are:

```text
lambda < lambda_*:  t = a +/- i sqrt(2(lambda_*-lambda))
lambda = lambda_*:  t = a, double real zero
lambda > lambda_*:  t = a +/- sqrt(2(lambda-lambda_*))
```

On the real-zero side,

```text
g(lambda) = 2 sqrt(2(lambda-lambda_*))
```

and

```text
g' = 4/g.
```

Therefore the local repulsive law is compatible with positive square-root birth. It cannot, by itself, prove `Lambda <= 0`.

## 7. Exact Implication Needed For RH

A non-circular proof through the Newman route must establish:

```text
There is no lambda_* > 0 and real t_* such that
H_{lambda_*}(t_*) = 0
partial_t H_{lambda_*}(t_*) = 0.
```

Equivalently, it must rule out a positive Newman boundary without assuming all zeros are real at `lambda = 0`.

## 8. Conditional Stability Statement

The local cluster and gap machinery may support conditional statements of the form:

```text
If all zeros are already real at some lambda_0,
and if the relevant external-field/comparison hypotheses hold,
then real-rootedness persists for larger lambda.
```

This is useful, but it is not a proof of RH if `lambda_0 = 0`, because that assumption is RH itself.

## 9. Numerical Evidence Status

The existing numerical material supports:

- reduced finite-window rank-2 structure;
- Hadamard alignment of the leading mode;
- density-balancing secondary mode;
- signed Hankel/Jensen patterns in tested finite ranges.

These are evidence and guideposts. They are not currently exact theorems about the full de Bruijn-Newman family.

## 10. Next Formal Target

The next non-circular theorem search should focus on one of:

```text
A. no positive Newman boundary;
B. Laguerre-Polya membership for H_lambda in the required lambda range;
C. signed total positivity / sign-regularity sufficient for real-rootedness;
D. Xi-specific comparison principle controlling the global tail.
```

Recommended first attack:

```text
C. signed total positivity / sign-regularity.
```

Reason:

The local zero-flow route has a known obstruction: positive square-root birth is locally compatible with the same repulsive law. The signed Hankel/Jensen evidence points toward a different global mechanism.
