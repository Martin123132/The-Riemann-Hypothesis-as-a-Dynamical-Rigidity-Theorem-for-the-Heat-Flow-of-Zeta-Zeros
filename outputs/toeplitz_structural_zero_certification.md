# Toeplitz Structural-Zero Certification Gate

Date: 2026-07-04

Status: finite certification-engine improvement. This is not a proof of RH. The promoted ranges listed in the certificate status note use coefficient enclosures plus Arb sign separation; active frontier attempts remain unpromoted until they pass the manifest gate.

## Purpose

The corrected PF route tests the ordinary Taylor coefficients:

```text
c_k(lambda) = A_k(lambda) / k!.
```

Previous Arb propagation found no negative intervals for `c_k`, but produced many zero-containing inconclusive intervals. Those were suspected to be structural zero minors of the upper-triangular Toeplitz matrix, not numerical sign failures.

This pass adds an exact structural-zero classifier.

## Structural-Zero Test

For a Toeplitz minor with row set `R` and column set `C`, each determinant term is represented symbolically as a monomial in variables:

```text
a_d, where d = column - row.
```

Entries with `d < 0` are forced zero. The classifier expands the determinant by permutations:

```text
sum_{permutations pi} sign(pi) product_i a_{C_{pi(i)} - R_i}
```

and combines identical monomials. If all monomials cancel, or all permutations hit forced zeros, the minor is structurally zero for every coefficient sequence.

This is exact combinatorics; it does not use zeta coefficients.

Fast criterion update:

For sorted row and column sets,

```text
rows = r_1 < ... < r_m
cols = c_1 < ... < c_m
```

the upper-triangular Toeplitz minor is structurally nonzero iff:

```text
c_i >= r_i for every i.
```

The scripts now use this criterion by default and keep the original permutation expansion as `--structural-mode exact`; `--structural-mode validate` checks both classifiers against each other. Validation on the `N = 10`, orders `<= 4` range matched the original structural counts exactly.

## Scripts Updated

```text
work/rh_compute/scripts/toeplitz_pf_audit.py
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

The exact audit now reports:

```text
structural_zero
point_zero_nonstructural
nonzero
```

The Arb probe now skips interval determinant evaluation for structural-zero minors and records them as:

```text
classification = structural_zero
```

## Exact Rational Audit Check

Run:

```text
toeplitz_pf_beta_taylor_N10_o4_struct
lambda = 0
N = 10
orders <= 4
rational_digits = 70
```

For both `beta = A_k` and `taylor = A_k/k!`, every exact point-zero minor in this range was structurally zero:

```text
order 1: structural_zero = 45, point_zero_nonstructural = 0
order 2: structural_zero = 1,200, point_zero_nonstructural = 0
order 3: structural_zero = 9,450, point_zero_nonstructural = 0
order 4: structural_zero = 30,240, point_zero_nonstructural = 0
```

For the corrected `taylor` sequence:

```text
positive = 19,690
structural_zero = 40,935
negative = 0
point_zero_nonstructural = 0
```

For the uncorrected `beta` sequence:

```text
negative = 4,589
```

So the zero classifier does not hide the known `A_k` Toeplitz failure.

## Arb Propagation After Structural-Zero Classification

All runs below use coefficient balls centered on the cached `repro_hankel_15c` coefficients with assumed absolute radius `1e-80` on each `A_k`.

### Taylor, Lambda = 0

```text
N = 10, orders <= 4:
  tests = 60,625
  positive = 19,690
  structural_zero = 40,935
  negative = 0
  inconclusive_contains_zero = 0

N = 10, orders <= 5:
  tests = 124,129
  positive = 39,094
  structural_zero = 85,035
  negative = 0
  inconclusive_contains_zero = 0

N = 12, orders <= 4:
  tests = 297,925
  positive = 88,309
  structural_zero = 209,616
  negative = 0
  inconclusive_contains_zero = 0
```

### Taylor, Lambda Grid

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 10
orders <= 4
```

the structural Arb result was:

```text
tests = 60,625
positive = 19,690
structural_zero = 40,935
negative = 0
inconclusive_contains_zero = 0
```

### Beta Red-Flag Control

For the uncorrected `beta = A_k` sequence:

```text
lambda = 0
N = 10
orders <= 3
positive = 5,127
structural_zero = 10,695
negative = 703
inconclusive_contains_zero = 0
```

This is the desired control: the improved classifier removes zero-minor interval noise without suppressing genuine negative minors.

## Interpretation

The finite PF certification workflow is now cleaner:

```text
structural-zero classifier
  -> ignore identically zero minors
  -> require Arb sign separation for all nonzero minors
```

Within the tested finite ranges, the corrected `c_k = A_k/k!` Toeplitz route has no detected negative or unresolved nonzero minors.

## Optimized Enumeration Update

The Arb propagation script now supports:

```text
--enumeration-mode nonzero
```

When paired with:

```text
--structural-mode fast
```

it enumerates only structurally nonzero minors, i.e. sorted row/column sets satisfying:

```text
c_i >= r_i for every i.
```

The skipped structural-zero minors are still included in the reported `tests` count and in `counts.structural_zero`; the summary also reports:

```text
evaluated_nonstructural_tests
```

as the number of determinant evaluations actually performed.

Planning utility:

```text
work/rh_compute/scripts/toeplitz_frontier_counts.py
```

This counts `tests`, `evaluated_nonstructural_tests`, and `structural_zero` by prefix DP before launching a large Arb run. It is not a sign test and not certificate evidence; it only estimates compute size.

Validation smoke check:

```text
lambda = 0
N = 10
orders <= 4
all mode and nonzero mode matched exactly:
tests = 60,625
positive = 19,690
structural_zero = 40,935
negative = 0
inconclusive = 0
```

Optimized lambda-grid frontier run:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 24
orders <= 4
tests = 117,085,204
evaluated_nonstructural_tests = 28,075,480
positive = 28,075,480
structural_zero = 89,009,724
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

Optimized matrix-size lambda-grid frontier run:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 26
orders <= 4
tests = 230,368,801
evaluated_nonstructural_tests = 54,414,126
positive = 54,414,126
structural_zero = 175,954,675
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

Optimized matrix-size lambda-grid frontier run:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 28
orders <= 4
tests = 430,101,469
evaluated_nonstructural_tests = 100,304,533
positive = 100,304,533
structural_zero = 329,796,936
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

Optimized matrix-size lambda-grid frontier extension:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 30
orders <= 4
tests = 767,707,750
evaluated_nonstructural_tests = 177,089,980
positive = 177,089,980
structural_zero = 590,617,770
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

Optimized higher-order lambda-grid frontier extension:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 22
orders <= 5
tests = 749,414,226
evaluated_nonstructural_tests = 161,341,895
positive = 161,341,895
structural_zero = 588,072,331
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

Optimized higher-order lambda-grid frontier extension:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 24
orders <= 5
tests = 1,923,675,220
evaluated_nonstructural_tests = 404,448,400
positive = 404,448,400
structural_zero = 1,519,226,820
negative = 0
inconclusive_contains_zero = 0
unknown = 0
zero = 0
recorded_problem_rows = 0
problem-row file size = 0
stderr size = 0
```

All five `N = 24`, orders `<= 5` summaries are promoted by the Toeplitz certificate manifest. This extends the structural-zero certificate frontier, but it remains a finite certificate range rather than a PF-infinity proof.

## GPU/Torch Prefilter

An auxiliary non-rigorous prefilter now exists:

```text
work/rh_compute/scripts/torch_toeplitz_prefilter.py
```

It uses the same structurally nonzero enumeration condition, but floating-point Torch determinants are not certificates. Its role is to find suspicious minors and prioritize Arb follow-up. Boundary and validation notes are recorded in:

```text
outputs/gpu_prefilter_and_certificate_boundary.md
```

## Manifest Checker

Promoted certificate claims are validated by:

```text
python work/rh_compute/scripts/check_toeplitz_certificate_manifest.py
```

The checker covers the finite ranges advertised in the status notes plus the
beta negative-control run, so future frontier expansions have a single
referee-facing audit gate.

## Remaining Hard Gate

For the promoted ranges, the coefficient balls are backed by the retained-integral and tail enclosure workflow. For any new active frontier, promotion still requires:

```text
interval-enclosed moment computation for c_k = mu_2k/(2k)!
Arb determinant sign separation
manifest validation
```

Only after those conditions are met can a structural Arb probe become a finite certificate. No finite certificate by itself proves PF-infinity or RH.
