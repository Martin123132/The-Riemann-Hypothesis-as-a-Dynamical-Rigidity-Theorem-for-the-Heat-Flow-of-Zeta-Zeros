# Finite Toeplitz/PF Certificate Status

Date: 2026-07-04

Status: selected finite Toeplitz/PF sign certificates for the corrected ordinary Taylor sequence `c_k = A_k/k!`. This is not a proof of RH, not an all-order PF-infinity proof, and not a Clay-ready argument.

## What Changed

Previously, the Arb Toeplitz/PF probe used coefficient balls:

```text
A_k(cache) +/- 1e-80
```

as propagation input, but those balls were assumed.

Now we have a coefficient-enclosure runner:

```text
work/rh_compute/scripts/acb_coefficient_enclosures.py
```

It combines:

```text
python-flint acb.integral for the retained finite integral;
analytic n-series tail bounds;
analytic u-cutoff tail bounds;
decimal cache-radius containment checks.
```

The retained integral is:

```text
2 integral_0^cutoff u^(2k) exp(lambda u^2) Phi_{n_sum}(u) du,
```

where `Phi_{n_sum}` is the finite `n = 1..n_sum` kernel sum.

## Coefficient Enclosure Runs

### Lambda Grid, k = 0..9

Run:

```text
acb_enclosures_repro_hankel_15c_lamgrid_k0_k9
lambdas = {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 0..9
n_sum = 100
cutoff = 8
dps = 120
abs_tol = 1e-90
cache_radius = 1e-80
```

Result:

```text
rows = 50
all cache balls A_k(cache) +/- 1e-80 contain the rigorous A_k balls
max A_k radius about 3.51e-93
```

This justifies the `1e-80` coefficient balls used by the structural Arb Toeplitz/PF probe for all `N = 10` Toeplitz matrices in the lambda grid, since those matrices only use coefficients `A_0..A_9`.

### Lambda Grid, k = 10..20

Run:

```text
acb_enclosures_repro_hankel_15c_lamgrid_k10_k20
lambdas = {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 10..20
n_sum = 100
cutoff = 8
dps = 120
abs_tol = 1e-90
cache_radius = 1e-80
```

Result:

```text
rows = 55
all cache balls A_k(cache) +/- 1e-80 contain the rigorous A_k balls
max A_k radius about 3.03e-102
```

Together with the previous `k = 0..9` run, this extends coefficient coverage to:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 0..20
```

### Lambda Grid, k = 21..32

Run:

```text
acb_enclosures_repro_hankel_15c_lamgrid_k21_k32
lambdas = {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 21..32
n_sum = 100
cutoff = 8
dps = 120
abs_tol = 1e-90
cache_radius = 1e-80
```

Result:

```text
rows = 60
all cache balls A_k(cache) +/- 1e-80 contain the rigorous A_k balls
max A_k radius about 1.16e-121
```

Together with the previous runs, this extends coefficient coverage to:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
k = 0..32
```

Thus all Toeplitz matrices with `N <= 33` in this lambda grid are covered at the coefficient-enclosure level for the current `repro_hankel_15c` cache.

## Certified Toeplitz/PF Finite Ranges

The structural-zero-aware Arb Toeplitz/PF probe already showed that, with coefficient balls `A_k(cache) +/- 1e-80`, the corrected sequence:

```text
c_k = A_k/k!
```

has no negative nonstructural Toeplitz minors in the following ranges.

Because the coefficient balls are now justified for the coefficients used, these ranges can be treated as finite certificates.

Performance/correctness update:

The structural-zero classifier now has a fast criterion:

```text
rows = r_1 < ... < r_m
cols = c_1 < ... < c_m
minor is structurally nonzero iff c_i >= r_i for all i.
```

This was validated against the original exact permutation expansion on the known `N = 10`, orders `<= 4` range. The fast classifier keeps the certified counts unchanged and makes larger Toeplitz probes practical.

Implementation update:

```text
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

now supports:

```text
--enumeration-mode all
--enumeration-mode nonzero
```

The default `all` mode preserves the original audit path by looping over every row/column pair. The `nonzero` mode requires `--structural-mode fast`, enumerates only column sets satisfying `c_i >= r_i`, and counts skipped forced structural zeros combinatorially. A smoke check on `lambda = 0`, `N = 10`, orders `<= 4` matched the all-pairs counts exactly:

```text
tests = 60,625
evaluated_nonstructural_tests = 19,690
positive = 19,690
structural_zero = 40,935
negative = 0
inconclusive = 0
```

Future larger probes should use `nonzero` mode after preserving occasional `all` or `validate` spot checks as audit controls.

### Lambda Grid, N = 10, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 10
orders <= 4
```

the certified structural Arb result is:

```text
tests = 60,625
positive = 19,690
structural_zero = 40,935
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 10, Orders <= 5

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 10
orders <= 5
```

the certified structural Arb result is:

```text
tests = 124,129
positive = 39,094
structural_zero = 85,035
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 12, Orders <= 5

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 12
orders <= 5
```

the certified structural Arb result is:

```text
tests = 925,189
positive = 258,193
structural_zero = 666,996
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 14, Orders <= 5

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 14
orders <= 5
```

the certified structural Arb result is:

```text
tests = 5,150,978
positive = 1,319,969
structural_zero = 3,831,009
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 16, Orders <= 5

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 16
orders <= 5
```

the certified structural Arb result is:

```text
tests = 22,720,080
positive = 5,471,960
structural_zero = 17,248,120
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 18, Orders <= 5

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 18
orders <= 5
```

the certified structural Arb result is:

```text
tests = 83,463,813
positive = 19,183,464
structural_zero = 64,280,349
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 20, Orders <= 5

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 20
orders <= 5
```

the certified structural Arb result is:

```text
tests = 265,184,141
positive = 58,773,841
structural_zero = 206,410,300
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 22, Orders <= 5

Using the optimized structural-zero enumeration:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 22
orders <= 5
enumeration_mode = nonzero
```

for each grid point, the certified structural Arb result is:

```text
tests = 749,414,226
evaluated_nonstructural_tests = 161,341,895
positive = 161,341,895
structural_zero = 588,072,331
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

### Lambda Grid, N = 24, Orders <= 5

Using the optimized structural-zero enumeration:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 24
orders <= 5
enumeration_mode = nonzero
```

for each grid point, the certified structural Arb result is:

```text
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

All five summaries are now promoted by `check_toeplitz_certificate_manifest.py`. This is still a finite certificate range, not an all-order PF-infinity theorem.

### Lambda Grid, N = 12, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 12
orders <= 4
```

the certified structural Arb result is:

```text
tests = 297,925
positive = 88,309
structural_zero = 209,616
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 14, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 14
orders <= 4
```

the certified structural Arb result is:

```text
tests = 1,142,974
positive = 317,968
structural_zero = 825,006
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 16, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 16
orders <= 4
```

the certified structural Arb result is:

```text
tests = 3,640,656
positive = 967,096
structural_zero = 2,673,560
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 18, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 18
orders <= 4
```

the certified structural Arb result is:

```text
tests = 10,053,189
positive = 2,578,680
structural_zero = 7,474,509
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 20, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 20
orders <= 4
```

the certified structural Arb result is:

```text
tests = 24,810,125
positive = 6,192,025
structural_zero = 18,618,100
negative = 0
inconclusive = 0
```

### Lambda Grid, N = 22, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 22
orders <= 4
```

the certified structural Arb result is:

```text
tests = 55,934,670
positive = 13,656,434
structural_zero = 42,278,236
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

### Lambda Grid, N = 24, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 24
orders <= 4
enumeration_mode = nonzero
```

the certified structural Arb result is:

```text
tests = 117,085,204
evaluated_nonstructural_tests = 28,075,480
positive = 28,075,480
structural_zero = 89,009,724
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

### Lambda Grid, N = 26, Orders <= 4

For each:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 26
orders <= 4
enumeration_mode = nonzero
```

the certified structural Arb result is:

```text
tests = 230,368,801
evaluated_nonstructural_tests = 54,414,126
positive = 54,414,126
structural_zero = 175,954,675
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

### Lambda Grid, N = 28, Orders <= 4

Using the optimized structural-zero enumeration:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 28
orders <= 4
enumeration_mode = nonzero
```

for each grid point, the certified structural Arb result is:

```text
tests = 430,101,469
evaluated_nonstructural_tests = 100,304,533
positive = 100,304,533
structural_zero = 329,796,936
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

### Lambda Grid, N = 30, Orders <= 4

Using the optimized structural-zero enumeration:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
N = 30
orders <= 4
enumeration_mode = nonzero
```

for each grid point, the certified structural Arb result is:

```text
tests = 767,707,750
evaluated_nonstructural_tests = 177,089,980
positive = 177,089,980
structural_zero = 590,617,770
negative = 0
inconclusive = 0
recorded_problem_rows = 0
```

## Negative-Control Check

The uncorrected exponential/Jensen coefficient sequence:

```text
beta_k = A_k
```

still fails under the same structural-zero machinery:

```text
lambda = 0
N = 10
orders <= 3
negative = 703
```

So the certification machinery is not hiding sign failures.

## Manifest Validation

The promoted finite certificate ledger is checked by:

```text
python work/rh_compute/scripts/check_toeplitz_certificate_manifest.py
```

This validates the advertised positive Toeplitz/PF summaries, their zero-byte
problem-row files, any available stderr logs, and the beta negative control.
As of the latest run it validates:

```text
95 promoted positive certificate summaries
1 negative control
```

## What This Proves

It proves selected finite Toeplitz/PF minor sign statements for the corrected ordinary Taylor sequence `c_k = A_k/k!`, within the ranges listed above.

It does not prove:

```text
c_k is PF-infinity;
all Toeplitz minors are nonnegative;
F_0 is Laguerre-Polya;
RH;
Lambda <= 0.
```

## Next Certification Expansion

1. Promote the remaining nonzero-lambda `N = 24`, orders `<= 5` exact Arb runs if they finish with the same clean counts, then decide whether to push `N = 26`, orders `<= 5`, or matrix size beyond `N = 30` at orders `<= 4`.
2. Preserve occasional `all` or `validate` spot checks while using `nonzero` mode for large matrix-size probes.
3. Use the first failure or first computational bottleneck to decide whether the PF route is plausibly structural or merely finite-pattern evidence.
4. Start theorem search for why `c_k = mu_2k/(2k)!` might be PF-infinity, now that the normalized finite evidence is certified across larger finite ranges.

Active frontier attempts and non-rigorous prefilters are tracked separately in:

```text
outputs/toeplitz_frontier_run_ledger.md
```
