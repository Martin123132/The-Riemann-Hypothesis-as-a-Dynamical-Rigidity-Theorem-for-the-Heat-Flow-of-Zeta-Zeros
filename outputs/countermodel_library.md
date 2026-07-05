# Countermodel Library

Date: 2026-07-05

Status: proof-safety artifact. This is not evidence against RH. It records small models, finite-prefix traps, and finite-grid traps that block invalid bridge lemmas in the RH dynamical-rigidity programme.

Executable gate:

```text
python work/rh_compute/scripts/countermodel_gate_examples.py
```

Current result:

```text
validated 11 countermodel gate examples
```

## Purpose

The programme is trying to prove the missing Newman direction:

```text
Lambda <= 0.
```

The dangerous failure mode is to turn one of these into a proof:

```text
local zero repulsion
finite signed-Hankel evidence
finite Jensen hyperbolicity
finite Jensen-window PF/Sturm rectangles
finite Toeplitz/PF certificates
finite Schur/Toeplitz shape-prefix evidence
finite Edrei moment-recurrence evidence
Stieltjes/Hankel moment positivity without Toeplitz total positivity
```

Each item is useful, but none is logically enough by itself. This file gives concrete gates that proposed proof steps must pass before they can enter a manuscript as more than evidence or a conditional claim.

## Gate 1: Local Heat Birth

The exact model:

```text
F_tau(t) = t^2 - 2 tau
```

satisfies the same Newman heat sign:

```text
partial_tau F = -2
partial_tt F = 2
partial_tau F = - partial_tt F
```

Its zeros are:

```text
tau < 0:  t = +/- i sqrt(-2 tau)
tau = 0:  double real zero at t = 0
tau > 0:  t = +/- sqrt(2 tau)
```

On the real-zero side, the gap is:

```text
g(tau) = 2 sqrt(2 tau)
```

and:

```text
g'(tau) = 4/g(tau).
```

Blocked proof step:

```text
The local law g' = 4/g excludes positive Newman birth.
```

Correct use:

```text
The local law describes the already-real side after a square-root birth.
```

Therefore any proof of `Lambda <= 0` must use global Xi structure, a sign-regularity theorem, or another nonlocal invariant. Local repulsion alone cannot do the job.

## Gate 2: Finite Prefix Is Not All-Order PF

Let:

```text
c_k(lambda) = A_k(lambda) / k!
```

be the ordinary coefficient sequence used in the coefficient-PF route. Suppose a finite prefix through `c_K` has passed every certified Toeplitz/PF test we have run.

That still cannot prove PF-infinity. Preserve the whole prefix through `K`, and define one positive next coefficient by:

```text
c_{K+1} = 2 c_K^2 / c_{K-1}.
```

Then the order-2 Toeplitz minor:

```text
det [[c_K,   c_{K+1}],
     [c_{K-1}, c_K]]
= c_K^2 - c_{K-1} c_{K+1}
= -c_K^2
< 0.
```

This does not say the actual zeta coefficient sequence fails. It says a proof step using only a finite prefix is invalid unless it adds a structural all-order theorem.

Blocked proof step:

```text
The certified finite Toeplitz/PF ledger makes c_k PF-infinity.
```

Correct use:

```text
The finite Toeplitz/PF ledger is falsification pressure and theorem-search evidence.
```

## Gate 2A: Finite Schur Prefix Is Not A Positive Specialization

The positive Schur-specialization target studies:

```text
h_k -> d_k(0)
```

and tries to prove all skew-Schur evaluations nonnegative. A finite set of
Schur or Toeplitz checks still cannot prove such a specialization exists.

The exact gate starts with:

```text
h_k = 1/k!
```

whose generating function is `exp(z)`, a clean restricted PF-infinity model.
It preserves:

```text
h_0, h_1, ..., h_6
```

and exactly checks a finite Toeplitz/Schur grid:

```text
N = 7
orders <= 4
2,940 finite tests
1,204 structurally nonzero positive minors
1,736 structural zeros
```

Then it chooses the next positive complete-homogeneous coordinate:

```text
h_7 = 1/2160
```

At the first untested Jacobi-Trudi shape:

```text
lambda = (6,6)
mu = (0,0)
```

the determinant becomes:

```text
s_(6,6) = det [[h_6, h_7],
               [h_5, h_6]]
        = -1/518400 < 0.
```

Blocked proof step:

```text
A finite Schur/Toeplitz shape ledger proves h_k -> d_k is a positive
specialization.
```

Correct use:

```text
Finite Schur/Toeplitz checks test proposed formulas. A proof needs an
all-order positive specialization, planar network, production matrix,
continued fraction, positive determinant integral, or equivalent theorem.
```

## Gate 3: Finite Signed-Hankel Is Not All-Order Signed Regularity

For:

```text
D_{m,s} = det(A_{i+j+s})_{i,j=0}^m
sigma(m) = (-1)^(m(m+1)/2)
```

the observed signed-Hankel condition is:

```text
sigma(m) D_{m,s} > 0.
```

For `m = 1`, this says:

```text
-(A_s A_{s+2} - A_{s+1}^2) > 0.
```

Given any positive prefix through `A_K`, preserve it and define:

```text
A_{K+1} = 2 A_K^2 / A_{K-1}.
```

At shift `s = K-1`:

```text
D_{1,K-1} = A_{K-1} A_{K+1} - A_K^2
          = A_K^2
          > 0
```

so:

```text
sigma(1) D_{1,K-1} = -A_K^2 < 0.
```

Blocked proof step:

```text
The finite signed-Hankel certificate grid proves all-order signed regularity.
```

Correct use:

```text
The finite signed-Hankel grid supports a conjectural all-order signed-regularity target.
```

## Gate 4: Finite Jensen Hyperbolicity Is Not Jensen Criterion

For degree 2 and shift `K-1`, the Jensen polynomial has the form:

```text
P_{2,K-1}(x)
= A_{K-1} + 2 A_K x + A_{K+1} x^2.
```

With the same positive one-term extension:

```text
A_{K+1} = 2 A_K^2 / A_{K-1}
```

the discriminant becomes:

```text
4(A_K^2 - A_{K-1} A_{K+1})
= -4 A_K^2
< 0.
```

So the next degree-2 Jensen polynomial is not hyperbolic, despite every earlier finite check being preserved.

Blocked proof step:

```text
Many finite Jensen/Sturm passes imply all Jensen polynomials are hyperbolic.
```

Correct use:

```text
Finite Jensen/Sturm passes guide theorem search and catch numerical failures.
```

## Gate 4A: Finite Consecutive Signed-Hankel Grid Is Not All Shifts

The executable gate also includes an exact finite-grid trap independent of
the zeta coefficients. Start with:

```text
a_k = 1/k!
```

For the shifted-principal signed-Hankel determinants:

```text
sigma(m) det(a_{i+j+s})_{i,j=0}^m
sigma(m) = (-1)^(m(m+1)/2)
```

the exact rational check validates the whole grid:

```text
m = 0..4
s = 0..8
45/45 signed determinants > 0
coefficients used: a_0..a_16
minimum signed grid value about 4.048066611965355946e-53
```

Now preserve all those coefficients and choose the next positive coefficient:

```text
a_17 = 1/167382319104000.
```

At the next untested shift:

```text
s = 15
```

the `m = 1` signed-Hankel value becomes negative:

```text
-2.284340357080483672e-27
```

and the degree-2 Jensen discriminant at the same shift is negative:

```text
-9.137361428321934688e-27.
```

Blocked proof step:

```text
The finite shifted-principal signed-Hankel grid proves all shifts, all-order
sign-regularity, or Jensen hyperbolicity.
```

Correct use:

```text
The finite grid is a certified finite diagnostic. A promotion needs a known
Hankel sign-consistency reduction plus an all-order proof, or a new bridge
theorem matching the zeta coefficient sequence.
```

## Gate 4B: Current Jensen-Window Rectangle Is Not All Shifts

The current Arb Jensen-window diagnostics are deliberately finite:

```text
PF obligation checks:
  degrees d = 3,4
  shifts n = 0..20
  coefficients used through A_24

Sturm/root-count checks:
  degrees d = 3,4,5
  shifts n = 0..20
  coefficients used through A_25
```

The executable gate preserves every coefficient those finite Jensen-window
checks can see:

```text
A_0, A_1, ..., A_25
```

for all five lambda rows, then chooses a positive next coefficient:

```text
A_26 = 2 A_25^2 / A_24.
```

At the next untested degree-2 Jensen window, shift `24`, the discriminant is:

```text
4(A_25^2 - A_24 A_26) = -4 A_25^2 < 0.
```

So all existing finite Jensen-window PF/Sturm inputs are preserved, but the
next degree-2 Jensen polynomial is not hyperbolic.

Blocked proof step:

```text
The finite Jensen-window PF/Sturm rectangle proves all-shift Jensen
hyperbolicity.
```

Correct use:

```text
The finite Jensen-window rectangle is a strong stress test and normalization
check. A proof still needs an all-degree/all-shift theorem for the actual
coefficient sequence.
```

## Gate 5: Finite Moment Recurrence Is Not An Edrei Representation

The Edrei reconstruction scout checks finite recurrence data for the shifted
moment sequence:

```text
a_n = p_{n+1}.
```

Positive recurrence data through any fixed order still does not prove an
all-order positive measure or Edrei zero-parameter representation.

The executable gate uses exact rational arithmetic. Start with the factorial
moments:

```text
m_n = n!
```

These are genuine Stieltjes moments of `exp(-x) dx` on `[0, infinity)`, so the
finite recurrence/Hankel prefix is honestly positive. For the current default
order `12`, preserve:

```text
m_0, m_1, ..., m_23
```

Then choose the next even moment `m_24` as a positive number below the exact
Schur-complement threshold for the next Hankel matrix. The gate reports:

```text
all preserved leading Hankel determinants > 0
adversarial m_24 > 0
next Hankel determinant < 0
```

Blocked proof step:

```text
The finite Arb recurrence scout through order 12 proves the all-order Edrei
moment representation.
```

Correct use:

```text
The recurrence scout is a constructive finite diagnostic and precision
frontier. It becomes a proof only after an all-order moment theorem,
positive parameter construction, or analytic recurrence formula is supplied.
```

## Gate 6: Stieltjes Moment Positivity Is Not Coefficient PF

The coefficient-PF route uses:

```text
c_k = mu_k / (2k)!.
```

A tempting but invalid bridge is:

```text
mu_k is a Stieltjes moment sequence
and 1/(2k)! has a restricted Laguerre-Polya generating function
therefore c_k is PF-infinity.
```

The executable gate blocks this with an exact positive measure:

```text
10 delta_0 + delta_1 + delta_2 on [0, infinity).
```

Its moments begin:

```text
mu_0..mu_6 = 12, 3, 5, 9, 17, 33, 65
```

and the leading Hankel determinants are:

```text
size 1..4 = 12, 51, 40, 0
```

so the moment sequence is Stieltjes/Hankel-nonnegative. But after the
coefficient-route normalization:

```text
c_0 = 12
c_1 = 3/2
c_2 = 5/24
```

the first order-2 Toeplitz/PF minor is:

```text
c_1^2 - c_0 c_2 = -1/4 < 0.
```

Blocked proof step:

```text
Stieltjes/Hankel positivity of the moments mu_k, together with the
factor 1/(2k)!, proves coefficient PF-infinity.
```

Correct use:

```text
Moment positivity can motivate the route, but a valid proof needs a
Toeplitz-total-positivity theorem, a positive determinant integral formula,
or an explicit restricted Laguerre-Polya factorization.
```

## Current Cache-Based Gate Run

The executable gate reads:

```text
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
```

and applies the finite-prefix traps to the five lambda rows:

```text
0
1e-6
1e-4
1e-2
1e-1
```

With the current cache, the prefix is preserved through:

```text
K = 32
```

For each lambda, a positive one-term extension at `K+1` breaks:

```text
order-2 Toeplitz/PF
m = 1 signed-Hankel
degree-2 Jensen hyperbolicity
```

The same run also validates the exact rational moment-recurrence trap:

```text
recurrence order <= 12
moments 0..23 preserved
positive edited moment 24
next Hankel/moment gate breaks
```

It validates the exact finite shifted-principal signed-Hankel grid trap:

```text
base a_k = 1/k!
m <= 4, shifts <= 8 preserved
positive a_17 breaks the next shifted m = 1 signed-Hankel/Jensen gate
```

It validates the exact finite Schur-prefix trap:

```text
base h_k = 1/k!
h_0..h_6 preserved
2,940 finite Toeplitz/Schur tests preserved
positive h_7 breaks s_(6,6)
```

It validates the finite Jensen-window rectangle trap:

```text
current Jensen-window PF/Sturm coefficient inputs A_0..A_25 preserved
positive A_26 breaks degree-2, shift-24 Jensen hyperbolicity
```

It also validates the Stieltjes multiplier trap:

```text
positive measure = 10 delta_0 + delta_1 + delta_2
leading Hankel determinants = 12, 51, 40, 0
c_1^2 - c_0 c_2 = -1/4
```

Again, these are proof-safety models, not claims about the actual zeta coefficients. Their role is logical: they show that finite-prefix evidence, finite shifted-principal grids, finite recurrence evidence, and generic Stieltjes/Hankel positivity cannot be promoted into an all-order theorem by wording alone.

## Manuscript Rule

Before accepting any proposed bridge lemma, ask:

```text
Does the lemma fail on the local heat-birth model?
Does the lemma rely only on a finite coefficient prefix?
Does the lemma rely only on a finite Schur/Toeplitz shape prefix?
Does the lemma rely only on a finite shifted-principal signed-Hankel grid?
Does the lemma rely only on a finite moment or recurrence prefix?
Does the lemma prove only Stieltjes/Hankel positivity when Toeplitz PF is required?
Does the lemma silently assume PF-infinity, Jensen hyperbolicity, or Laguerre-Polya membership?
Does the lemma assume RH at lambda = 0?
```

If yes, it cannot be used to prove `Lambda <= 0`.

The only acceptable promotions are:

```text
finite certificate
conditional theorem with explicit hypotheses
all-order theorem with noncircular structural proof
```

Executable result-language scan:

```text
python work/rh_compute/scripts/check_output_reference_integrity.py
python work/rh_compute/scripts/check_output_status_manifest.py
python work/rh_compute/scripts/check_proof_claim_ledger.py
python work/rh_compute/scripts/check_result_language_boundaries.py
```

Current result:

```text
validated output references: scanned 40 markdown files, 396 path references, 0 missing required paths, 3 planned missing deliverables
validated output artifact statuses: scanned 40 markdown files, 0 status issues
validated proof-claim ledger: 24 claims, 0 issues, 6 open theorem targets
validated result-language boundaries: scanned 40 markdown files, 0 overclaims
```
