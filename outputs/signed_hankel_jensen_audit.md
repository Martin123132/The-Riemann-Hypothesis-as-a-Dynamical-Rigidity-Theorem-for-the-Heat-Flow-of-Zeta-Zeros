# Signed Hankel / Jensen Audit

Date: 2026-07-03

Status: finite-evidence audit. This is not a proof of RH or `Lambda <= 0`; it records signed-Hankel/Jensen diagnostics and their missing bridge theorem.

Purpose:

Extract the exact kernel, coefficient, Jensen-polynomial, Sturm, and signed-Hankel objects from the existing bundle. This is a definition audit, not an independent rerun of the numerical claims.

Source notes used:

- `Kernel and Jensen Diagnostics Checkpoint`
- `Checkpoint 5th may`
- `## 12. Lambda-Grid Sturm Audit for Degrees 21-30`
- `## 13. Lambda-Grid Sturm Audit for Degrees 21-30, n=21,...,40`
- relevant code/output cells in `Untitled247 (5).ipynb`

## 1. Kernel Convention

The working kernel is

```text
Phi(u)
= sum_{r >= 1}
  (2*pi^2*r^4*exp(9u) - 3*pi*r^2*exp(5u))
  * exp(-pi*r^2*exp(4u)).
```

The even extension is used:

```text
Phi(-u) = Phi(u).
```

For positive Newman deformation parameter `lambda`, the tested deformed kernel is

```text
Phi_lambda(u) = exp(lambda*u^2) Phi(u).
```

## 2. Moment And Coefficient Sequence

To remove an indexing ambiguity in the notes, use `mu_{2k}` for the raw even moment:

```text
mu_{2k}(lambda)
= integral_R u^(2k) Phi_lambda(u) du.
```

The Jensen coefficient sequence is

```text
A_k(lambda)
= mu_{2k}(lambda) * k! / (2k)!.
```

This matches the notebook code:

```python
A.append(M * mp.factorial(k) / mp.factorial(2 * k))
```

where `M` is the kth stored even moment, i.e. `mu_{2k}`.

## 3. Jensen Polynomials

For degree `d` and shift `n`, the finite Jensen polynomial is

```text
P_{d,n}(x)
= sum_{k=0}^d binom(d,k) A_{n+k} x^k.
```

The tested polynomial is

```text
Q_{d,n}(y) = P_{d,n}(-y).
```

The pass criterion is:

```text
Q_{d,n}(y) has exactly d positive real roots.
```

Equivalently:

```text
P_{d,n}(x) has d real nonpositive roots.
```

The high-degree audits use Sturm positive-root counting on `Q_{d,n}` rather than floating root classification.

## 4. Recorded Jensen/Sturm Evidence

These are recorded results from the bundle, not independently rerun in this audit.

### Lambda = 0

Recorded clean ranges:

```text
degrees 2..20,  n = 0..80: pass
degrees 21..30, n = 0..40: pass
```

The degree `16`, `n=54`, `lambda=0` floating-root failure was reclassified in the notes as a numerical conditioning artifact after Sturm/root verification.

### Positive Lambda Spot Checks

Recorded clean range:

```text
lambda in {1e-4, 1e-2, 1e-1}
degrees 16..20
n = 0..80
pass: 1215/1215
```

### Full Lambda Grid For Degrees 21..30

Recorded lambda grid:

```text
0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1
```

Recorded tested range:

```text
degrees 21..30
n = 0..40
```

Recorded result:

```text
3280/3280 pass
failures: 0
```

## 5. Hankel Objects

For the coefficient sequence `A_k(lambda)`, define the shifted Hankel determinant:

```text
D_{m,s}(lambda)
= det( A_{i+j+s}(lambda) )_{i,j=0}^m.
```

The observed signed pattern is

```text
sigma(m) D_{m,s}(lambda) > 0,
```

where

```text
sigma(m) = (-1)^(m(m+1)/2).
```

Equivalently, the determinant sign pattern by `m mod 4` is:

```text
m = 0 mod 4: positive
m = 1 mod 4: negative
m = 2 mod 4: negative
m = 3 mod 4: positive
```

## 6. Recorded Signed-Hankel Evidence

### Checkpoint 15B

Recorded test:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts s = 0,1
Hankel sizes m = 0..14
```

Recorded result:

```text
150/150 signed determinants positive
```

### Checkpoint 15C

Recorded test:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts s = 0..8
Hankel sizes m = 0..12
```

Recorded result:

```text
585/585 signed determinants positive
failures: 0
```

This is the main finite evidence for multi-shift signed Hankel sign-regularity.

Certificate update added 2026-07-04:

The reproduced grid was first interval-certified through the original Checkpoint 15C range, then extended after adding rigorous coefficient enclosures through `k = 64`. The currently promoted certificate grid is:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
m = 0..19
s = 0..24
```

It is interval-certified using the rigorous `A_ball` coefficient enclosures through `k = 64`.

```text
work/rh_compute/scripts/arb_hankel_enclosure_sign_probe.py
work/rh_compute/scripts/check_hankel_certificate_manifest.py
```

The manifest checker validates:

```text
2,500 signed-Hankel finite certificates
```

This upgrades the finite signed-Hankel grid from high-precision reproduction to finite interval certificate. It still does not prove the all-order signed-Hankel hypothesis or any implication to Jensen hyperbolicity.

Countermodel update added 2026-07-04:

```text
work/rh_compute/scripts/countermodel_gate_examples.py
outputs/countermodel_library.md
```

Using the current five-lambda coefficient cache, the gate preserves each finite prefix through `k = 32` and constructs positive one-term extensions that break `m = 1` signed-Hankel and degree-2 Jensen hyperbolicity. These are not models of the zeta coefficients; they are logical guards showing that finite-prefix signed-Hankel/Jensen success cannot be promoted to an all-order bridge theorem without additional structure.

The same executable gate now includes an exact finite-grid trap independent of
the zeta coefficients: `a_k = 1/k!` satisfies the shifted-principal
signed-Hankel grid `m = 0..4`, `s = 0..8`, but a positive `a_17` extension
preserving all coefficients used by that grid breaks the next `m = 1`
signed-Hankel determinant and the matching degree-2 Jensen discriminant. This
blocks the narrower promotion from a finite shifted-principal grid to all
shifts or all-order sign consistency.

## 7. Failed Or Damaged Routes

### Ordinary Log-Concavity

The sampled `Q_2` diagnostic in the notes did not support naive ordinary log-concavity.

Status:

```text
weak / not current route
```

### Simple Coefficient PF-Infinity Toeplitz Positivity

The notes report robust negative order-3 Toeplitz minors from the coefficient sequence `A_k`.

Representative described pattern:

```text
rows = (0,1,2)
cols = (1,2,3)
det approx -1.7e-11
```

Status:

```text
damaged as a naive route
```

Correction added 2026-07-03:

This damage applies to the exponential/Jensen coefficients `A_k`. The classical ASW/Edrei PF-sequence theorem is formulated for ordinary Taylor coefficients. Since

```text
F_lambda(s) = sum A_k(lambda) s^k / k!,
c_k(lambda) = A_k(lambda) / k!,
```

the PF-infinity Toeplitz route must be re-audited for `c_k = A_k/k!`. A new exact-rational finite audit found:

```text
lambda = 0, N = 10, orders <= 4
A_k:      4,589 negative minors among 60,625
A_k/k!:  0 negative minors among 60,625
```

So the corrected status is: `A_k` naive PF is damaged; `A_k/k!` PF is an active theorem-search route, still unproved.

### Ordinary Stieltjes / Orthogonal-Polynomial Bridge

Checkpoint 16 reports that the ordinary Hankel-to-orthogonal-polynomial bridge failed:

```text
polynomial root failures: 225/250
interlacing failures: 225/225
```

Status:

```text
ordinary Stieltjes route failed
```

Interpretation:

The observed structure is not ordinary positive moment theory. It appears, if real, to be a rotated, signed, or indefinite sign-regular structure.

## 8. Active Hypothesis

The active coefficient-side hypothesis is:

```text
For all m >= 0, s >= 0, and lambda in the required range,
(-1)^(m(m+1)/2)
det( A_{i+j+s}(lambda) )_{i,j=0}^m
> 0.
```

Call this:

```text
Signed Hankel Sign-Regularity Hypothesis.
```

## 9. Why This Is Not Yet A Proof

The missing theorem is not merely another finite test.

The needed theorem would be something like:

```text
Signed Hankel Sign-Regularity Hypothesis
  + Xi/Phi-specific analytic hypotheses
  => Jensen hyperbolicity for all d,n
  => H_lambda in Laguerre-Polya / no positive Newman boundary.
```

At present, the bundle has:

```text
finite evidence for signed Hankel signs
finite Sturm evidence for Jensen hyperbolicity
ordinary positivity routes eliminated
```

It does not yet have:

```text
a theorem connecting the signed Hankel pattern to all-degree/all-shift Jensen hyperbolicity
or to de Bruijn-Newman real-rootedness.
```

## 10. Next Compute Sprint

Build a clean command-line reproduction script with stable outputs:

```text
work/rh_compute/scripts/jensen_hankel_runner.py
```

Minimum requirements:

- parameters: lambda grid, max degree, n range, max Hankel size, shift range, precision;
- output JSONL rows for every test;
- output CSV summaries;
- write coefficient cache for `A_k(lambda)`;
- implement Sturm root counting for `Q_{d,n}`;
- implement signed Hankel determinant tests;
- record failures and suspicious near-zero determinants.

Initial reproduction targets:

```text
Jensen:
  lambda grid = 0,1e-6,1e-5,1e-4,1e-3,1e-2,5e-2,1e-1
  d = 21..30
  n = 0..40

Hankel:
  lambda grid = 0,1e-6,1e-4,1e-2,1e-1
  m = 0..12
  s = 0..8
```

Only after clean reproduction should the upgraded laptop be used for larger sweeps.

## 11. Next Theory Sprint

Search specifically for the correct theorem ecosystem:

```text
signed Hankel matrices
sign-regular Hankel matrices
oscillatory kernels
indefinite moment problems
J-fractions / continued fractions with alternating signs
total positivity after diagonal sign conjugation
Pade / Stieltjes transforms with signed measures
```

Question to answer:

```text
Is there a known sign-regular Hankel criterion that implies real-rootedness,
interlacing, or Laguerre-Polya membership after a transformation matching A_k?
```

If yes, try to match its hypotheses.

If no, state the missing theorem explicitly.
