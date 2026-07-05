# Bridge Theorem Options

Date: 2026-07-03

Status: theorem-route map. This is not a proof of RH or `Lambda <= 0`; it lists noncircular bridge targets and failure modes.

Purpose:

Identify non-circular theorem routes that could, in principle, connect the current dynamical-rigidity programme to a proof of RH through the de Bruijn-Newman family.

This file starts from `formal_core.md`.

## Required Endgame

The route must prove

```text
Lambda <= 0
```

without assuming RH at `lambda = 0`.

Equivalently, it must rule out a positive Newman boundary:

```text
lambda_* > 0,
H_{lambda_*}(t_*) = 0,
partial_t H_{lambda_*}(t_*) = 0.
```

Any theorem that assumes all zeros are real at `lambda = 0` is conditional, not a proof of RH.

## Option A: No Positive Newman Boundary

### Statement To Prove

There is no positive `lambda_*` at which `H_lambda` has a real multiple zero:

```text
not exists lambda_* > 0, t_* real:
H_{lambda_*}(t_*) = 0
and partial_t H_{lambda_*}(t_*) = 0.
```

### Why It Would Work

A positive Newman boundary is exactly the obstruction to `Lambda <= 0`. Excluding such a boundary directly proves the missing Newman direction.

### What We Have

- Local square-root normal form.
- Local cluster non-collapse on the already-real side.
- Numerical local no-escape/winding experiments in finite windows.

### Main Obstacle

Local repulsion does not exclude positive square-root birth. The toy model

```text
F(lambda,t) = (t-a)^2 - 2(lambda-lambda_*)
```

satisfies the same heat sign and has the same `g' = 4/g` law on the real side after birth.

### Verdict

This is the cleanest target but probably too direct for the next sprint. We need a global invariant or sign structure first.

## Option B: Laguerre-Polya Membership

### Statement To Prove

Show that the relevant entire function belongs to the Laguerre-Polya class in the needed lambda range, for example:

```text
H_lambda in LP for all lambda > 0,
with enough closure at lambda = 0 to imply Lambda <= 0.
```

### Why It Would Work

Real entire functions in the Laguerre-Polya class are locally uniform limits of real-rooted polynomials and have only real zeros.

### What We Have

- The de Bruijn-Newman family is already framed as a real entire heat deformation.
- The signed Hankel/Jensen evidence may be a coefficient-side shadow of an LP mechanism.

### Main Obstacle

LP membership is essentially the desired real-rootedness statement. We need a more checkable sufficient criterion, not just a restatement.

### Verdict

Good as a language and endpoint. Not yet good as an attack route unless paired with a concrete positivity/sign-regularity criterion.

## Option C: Signed Total Positivity / Sign-Regularity

### Statement To Prove

Find a transformed coefficient sequence or kernel associated with `H_lambda` whose signed minors obey a uniform sign-regularity theorem strong enough to imply LP membership or real-rooted Jensen polynomial hyperbolicity for all degrees and shifts.

Schematic target:

```text
for all admissible minors M:
sign(det M) = prescribed_sign(M)
and this sign-regularity implies H_lambda in LP / no positive boundary.
```

### Why It Might Work

The later notebook evidence points away from ordinary positivity and toward a signed Hankel pattern. Classical total positivity and Polya-frequency theory connect coefficient matrices, variation-diminishing properties, and real-rootedness.

### What We Have

- Finite signed Hankel/Jensen tests.
- Lambda-grid Sturm audits in the bundle.
- Evidence that ordinary Stieltjes positivity is not the right formulation.

### Main Obstacles

- The exact coefficient sequence must be defined cleanly.
- The tested signed Hankel matrices may not be the right matrices for known PF/TP theorems.
- Finite tests do not imply uniform sign-regularity.
- Known Jensen asymptotics are not enough for RH; all degrees and all shifts matter.

### Verdict

Recommended first serious route.

It is the most aligned with the new evidence while avoiding the known local-repulsion obstruction. But it must become a theorem about all minors/degrees/shifts, not just larger numerical sweeps.

The current formal version of this target is maintained in:

```text
outputs/signed_hankel_jensen_bridge_target.md
python work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py
```

## Option D: Xi-Specific Tail / Gap Comparison Principle

### Statement To Prove

For actual `H_lambda` zeros, prove a gap comparison theorem that rules out positive boundary birth.

Naive target:

```text
g_i' >= 4/g_i
```

is not viable for the finite pairwise reduced model because the pairwise tail is compressive. A viable version must be Xi-specific:

```text
g_i' >= 4/g_i - controlled_tail_i,
```

with a bound strong enough to prevent positive square-root birth globally, or a different comparison principle bypassing this form.

### Why It Might Work

The dynamical project naturally wants a gap theorem. If such a comparison principle exists, it would connect directly to zero motion.

### Main Obstacle

The local law and finite pairwise model do not give it. The proof would have to use global structure of the Xi kernel, not just local zero repulsion.

### Verdict

Important but probably second priority. Use countermodels aggressively before trusting any proposed comparison statement.

## Literature Anchors

The following references shape the next search.

- Clay Mathematics Institute: RH is a Millennium Prize Problem and asserts all nontrivial zeta zeros lie on the critical line.
  <https://www.claymath.org/millennium/riemann-hypothesis/>

- Rodgers-Tao prove `Lambda >= 0`, while RH is equivalent to `Lambda <= 0`.
  <https://arxiv.org/abs/1801.05914>

- Griffin-Ono-Rolen-Zagier prove major asymptotic Jensen polynomial hyperbolicity results and recall Polya's Jensen criterion.
  <https://arxiv.org/abs/1902.07321>

- Griffin-Ono-Rolen-Thorner-Tripp-Wagner make effective Jensen polynomial results for the xi function.
  <https://arxiv.org/abs/1910.01227>

- Farmer argues that Jensen-polynomial asymptotic routes are not a plausible path to RH by themselves.
  <https://arxiv.org/abs/2008.07206>

- Classical Polya-frequency / total-positivity theory is the likely theorem ecosystem for the signed-minor route. Search anchors include Aissen-Schoenberg-Whitney, Edrei, Schoenberg, Karlin, and modern total positivity preservers.

## Recommendation

Proceed with Option C first:

```text
Signed Total Positivity / Sign-Regularity
```

but with a strict discipline:

1. Define the coefficient sequence exactly.
2. Define the tested matrices exactly.
3. Reproduce the finite evidence cleanly.
4. Search for a theorem that turns the exact signed pattern into real-rootedness.
5. If no known theorem fits, state the missing theorem precisely.
6. Only then expand computation.

## Next File To Build

```text
signed_hankel_jensen_audit.md
```

It should extract the exact objects from the notebooks:

- coefficient sequence;
- Jensen polynomial normalization;
- lambda-grid definition;
- Hankel matrix definition;
- shifted Hankel matrix definition;
- sign convention;
- finite pass/fail ranges;
- numerical precision and certification method.
