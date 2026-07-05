# Clay Prize Readiness Audit

Date: 2026-07-03

Status: readiness audit. This is not a proof of RH or `Lambda <= 0`; it separates the current programme from a Clay-ready argument.

Source bundle reviewed:

- `The Riemann Hypothesis as a Dynamical Rigidity Theorem for the Heat Flow of Zeta Zeros`
- `The Riemann Hypothesis as a Dynamical Rigidity Problem`
- `STEP 9 - corrected reproducibility package`
- Parts 14-19 and the README

External status checks:

- Clay Mathematics Institute still lists the Riemann Hypothesis as a Millennium Prize Problem: <https://www.claymath.org/millennium/riemann-hypothesis/>
- Rodgers and Tao prove the complementary Newman bound `Lambda >= 0`; RH is equivalent to `Lambda <= 0`: <https://arxiv.org/abs/1801.05914>

## Verdict

This bundle is not a Clay-ready proof of the Riemann Hypothesis.

It contains useful exact local zero-motion identities and substantial numerical evidence for a reduced repulsive zero-flow geometry, but the argument does not prove the global statement needed for RH. The strongest later notes in the bundle already identify the central missing theorem: control of the nonlocal/global tail in the actual de Bruijn-Newman zero flow, or an equivalent theorem excluding a positive Newman boundary.

The current best framing is:

> A serious research programme around zero dynamics and Newman-flow rigidity, with a sharp global bridge problem still open.

not:

> A completed proof of RH.

## The Prize-Critical Gap

The de Bruijn-Newman setup gives a constant `Lambda` such that the deformed xi function has all real zeros precisely in the real-zero regime. RH is equivalent to

```text
Lambda <= 0.
```

Rodgers-Tao proves

```text
Lambda >= 0.
```

So a Clay-winning argument through this route must prove

```text
Lambda <= 0,
```

equivalently, it must rule out every positive Newman boundary

```text
lambda_* > 0
```

where a complex conjugate pair becomes a double real zero and then splits into two real zeros as lambda increases.

The local repulsive zero-flow law does not rule this out. It describes exactly what happens on the real-zero side after such a birth.

## Minimal Countermodel To The Local-Barrier Strategy

Consider the one-parameter entire function

```text
F_lambda(t) = (t - a)^2 - 2(lambda - lambda_*).
```

It satisfies the same Newman heat sign:

```text
partial_lambda F = - F_tt.
```

Its zeros are:

```text
lambda < lambda_*:  t = a +/- i sqrt(2(lambda_* - lambda))
lambda = lambda_*:  t = a, a double real zero
lambda > lambda_*:  t = a +/- sqrt(2(lambda - lambda_*))
```

On the real-zero side, the zero velocity is

```text
t'_+ =  1 / sqrt(2(lambda - lambda_*))
t'_- = -1 / sqrt(2(lambda - lambda_*))
```

and the real gap

```text
g(lambda) = t_+ - t_-
```

satisfies

```text
g' = 4 / g.
```

This is exactly the repulsive adjacent-pair law used in the manuscript.

Therefore, no argument based only on the local statement

```text
g' ~ 4/g
```

can exclude a positive Newman boundary. The local law is compatible with square-root birth.

## Tail Positivity Cannot Be The Reduced-Model Theorem

For the finite repulsive model

```text
x_i' = 2 sum_{j != i} 1/(x_i - x_j),
```

with ordered zeros and neighboring gap

```text
g_i = x_{i+1} - x_i,
```

one computes

```text
g_i'
= 4/g_i
  - 2 g_i sum_{j != i,i+1}
      1 / ((x_{i+1}-x_j)(x_i-x_j)).
```

For every outside zero `x_j`, the product

```text
(x_{i+1}-x_j)(x_i-x_j)
```

is positive. Hence the finite-window tail term is compressive, not positive.

So a bridge lemma of the form

```text
R_i(lambda) >= 0
```

cannot be true as a statement about the reduced pairwise tail. A viable theorem would need either:

- a full Xi-specific external term that reverses or controls this compression,
- a sharp lower bound showing the negative tail is too small near every possible global birth,
- or a different global principle entirely.

## The Circular Conditional

The "Candidate Lemma X" in the corrected reproducibility package assumes:

```text
RH at lambda = 0, so all zeros are real initially.
```

But in the de Bruijn-Newman framework, RH is already equivalent to `Lambda <= 0`. Thus this assumption cannot be part of a proof of RH. It can at most prove a conditional stability statement:

```text
If RH holds at lambda = 0, then real-rootedness persists for lambda > 0.
```

That conditional direction is not the missing Clay theorem.

## What Is Solid

The following parts look mathematically meaningful and worth preserving:

- zero tracking by the implicit function theorem for simple real zeros;
- the exact velocity identity `t_i' = -F_lambda/F_t`;
- the Hadamard/local product derivation of the singular adjacent-pair term;
- Weierstrass square-root normal form near a simple double-zero birth;
- the local statement that already-real finite isolated clusters do not collapse in the increasing-lambda repulsive direction;
- the numerical evidence that finite windows of the reduced model exhibit a robust low-dimensional geometry.

These are not enough for RH, but they are a coherent research programme.

## The Actual Theorem Needed

A route through this programme would need to prove one of the following, with no RH assumption at lambda = 0.

### Target A: No Positive Newman Boundary

Prove that there is no `lambda_* > 0` and real `t_*` such that

```text
Xi_{lambda_*}(t_*) = 0
partial_t Xi_{lambda_*}(t_*) = 0.
```

This would exclude the positive square-root birth scenario and imply `Lambda <= 0`.

### Target B: Laguerre-Polya Membership For All Positive Lambda

Prove that `Xi_lambda` belongs to the Laguerre-Polya class for every `lambda > 0`.

This is essentially the desired real-rootedness statement in entire-function language.

### Target C: Global Sign-Regularity Or Total Positivity

Turn the observed signed Hankel/Jensen patterns into a uniform theorem strong enough to imply real-rootedness of the limiting entire function, not just finite Jensen polynomial hyperbolicity.

This is the most natural continuation of the later notebook evidence.

2026-07-04 update:

The finite signed-Hankel grid has been upgraded from reproducible numerical evidence to enclosure-backed interval certificate:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
m = 0..19
s = 0..24
2,500/2,500 signed-Hankel determinants sign-separated positive
```

This improves the evidence quality but does not change the Clay-readiness conclusion: the missing ingredient is still a theorem turning finite or all-order signed-Hankel/sign-regularity into Jensen hyperbolicity, Laguerre-Polya membership, or `Lambda <= 0`.

### Target D: Xi-Specific Tail/Comparison Principle

For actual de Bruijn-Newman zeros, prove a comparison theorem for neighboring gaps that rules out square-root birth at positive lambda. It cannot be a merely local `4/g` argument; it must use global Xi structure.

## Recommended Manuscript Revision

Do not submit the current manuscript as a proof of RH.

Revise it as:

```text
The Riemann Hypothesis as a Dynamical Rigidity Programme:
Local Zero Motion, Reduced Repulsive Flow, and the Newman Boundary Gap
```

The main theorem should be local/conditional. The central open problem should be stated explicitly:

```text
Open bridge problem:
Exclude a positive de Bruijn-Newman square-root birth by proving a global
Xi-specific sign-regularity, total-positivity, or tail-comparison theorem.
```

That version would be honest, technically sharper, and far more useful to a serious reader.
