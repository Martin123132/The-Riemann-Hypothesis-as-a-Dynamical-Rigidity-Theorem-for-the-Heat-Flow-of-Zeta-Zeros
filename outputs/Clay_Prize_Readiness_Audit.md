# Clay Prize Readiness Audit

Date: 2026-07-16

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

2026-07-16 update:

The fixed-order signed-Hankel programme is now rigorous through order nine on
the complete heat segment:

```text
epsilon_m H_(m,n)(lambda)>0
for every 1<=m<=9, n>=0, and -100<=lambda<=0,
```

including the consecutive-row arbitrary-column transfer at each of those
orders. The order-nine endpoint proof combines a rigorous prefix, a global
continuous sixth-nested first-summand curvature theorem, an exact
complete-to-first-summand transfer, and a two-index finite splice.

More importantly, the repeated heat-flow step has been proved once at
arbitrary order. The all-fixed-order eventual tails, affine determinant heat
identity, and flag-Plucker recursion give the exact reduction

```text
[Q_(m,n)(-100)>0 for every m>=10,n>=0]
iff
[Q_(m,n)(lambda)>0 for every m>=10,n>=0
 and -100<=lambda<=0].
```

The thresholds in the eventual-tail theorem may depend on `m`; ordinary
induction fixes one order, uses its own threshold, completes that layer, and
then moves upward. Thus a common tail threshold uniform in order is not a
necessary part of the dynamical argument.

This is a real narrowing of the conditional dynamical argument, but the
endpoint antecedent has now been tested at its first previously open order
and found false. Rigorous `4096`-bit Arb determinant evaluations give

```text
Q_(10,n)(-100)<0 for n=0,1,2,3,
Q_(10,4)(-100)>0.
```

Direct Hankel, Jacobi-Trudi, and Toda computations agree, and all relevant
Arb balls are sign-separated. Therefore the proposed all-order
signed-Hankel endpoint hierarchy is not an unproved obligation: it is a
rejected antecedent. The endpoint-to-heat equivalence remains a valid
conditional theorem, and all fixed-order results through order nine remain
valid.

The static endpoint obligation now has an exact normalized Schur coordinate.
For

```text
h_k=A_k(-100)/A_0(-100),
```

column reversal and Jacobi-Trudi prove

```text
Q_(m,n)(-100)=A_0(-100)^m s_((n+m-1)^m)(h).
```

Arbitrary columns correspond bijectively to partitions
`lambda_1>=...>=lambda_m>=m-1`. Consequently the rejected endpoint
antecedent is equivalent, using the completed lower base and fixed-order
initial-minor transfer, to

```text
s_((N^m))(h)>0 for every m>=10 and N>=m-1,
```

or to positivity on the entire associated deep partition cone. At order ten,
the required rectangles `(9^10),(10^10),(11^10),(12^10)` are all strictly
negative; `(13^10)` is positive. A certified finite scan finds positivity for
every shift `4<=n<=1240`, but four failures are enough to refute the
all-shift statement.

This is not an ordinary PF-infinity problem. Rigorous endpoint coefficient
enclosures prove

```text
s_(1,1,1)(h)<0,
```

so the actual endpoint sequence is not even `PF_3`. That earlier failed shape
lies outside the deep cone at order three and rules out the stronger Edrei/PF
shortcut; the new order-ten failures lie inside the deep cone and reject the
rectangle hierarchy itself. The Clay-readiness verdict remains negative,
but the reason is sharper: this route cannot be completed by proving its
former endpoint obligation because that obligation is false.

The newest exact gate narrows both obligations further. For
`tau_(m,N)=s_((N^m))(h)`, Desnanot-Jacobi gives

```text
tau_(m+1,N) tau_(m-1,N)
 =tau_(m,N)^2-tau_(m,N-1) tau_(m,N+1).
```

Once the lower row is positive, the next rectangle is therefore equivalent
to strict width-log-concavity of the current row. This is a useful Toda
coordinate, but not a propagation theorem: the decisive right-hand side is a
difference of positive terms. At the first previously open step it is
strictly negative for four admissible widths.

Two plausible promotions are now excluded exactly. First, resetting negative
indices in the ordinary shifted tail changes shallow Jacobi-Trudi
determinants; at `r=2,s=1,mu=(0,0)` the defect is
`h_0 h_2/h_1^2>0`. Deep positivity therefore does not imply moving-tail PF by
that translation. Second,

```text
H(z)=exp(z/100)/((1-z)(1-2z))
```

is strictly Schur-positive for every partition but has a cubic Jensen
polynomial with exact negative discriminant. Thus even strict full
unweighted Schur/PF positivity cannot serve as a generic Jensen bridge. Any
successful conversion must exploit additional structure specific to the
actual Xi/Phi coefficient sequence. Together with the direct order-ten
failure, these are route eliminations, not negative evidence against RH or
Jensen hyperbolicity itself.

The remaining live research burden is therefore to find a weaker
Xi/Phi-specific Jensen mechanism, a different determinant/kernel condition,
or a direct Newman-flow comparison theorem that does not assume all-shift
signed-Hankel positivity. No such closing theorem is currently proved.

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
Xi-specific Jensen-window, determinant/kernel, or tail-comparison theorem
whose hypotheses are weaker than the rejected all-shift signed-Hankel cone.
```

That version would be honest, technically sharper, and far more useful to a serious reader.
