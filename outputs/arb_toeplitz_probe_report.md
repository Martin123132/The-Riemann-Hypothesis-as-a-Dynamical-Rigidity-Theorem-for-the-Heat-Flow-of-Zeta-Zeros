# Arb Toeplitz/PF Propagation Probe

Date: 2026-07-03

Status: propagation diagnostic only. This is not a rigorous coefficient certificate, because the input coefficient balls are centered on cached numerical coefficients with an assumed absolute radius.

Update:

The earlier zero-containing inconclusive rows have now been resolved for the tested finite ranges by adding an exact structural-zero classifier to the Arb probe. See `outputs/toeplitz_structural_zero_certification.md`.

## Script

```text
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

Purpose:

```text
Given cached A_k(lambda) values and an assumed absolute radius on each A_k,
build Arb balls for either:

  beta:   a_k = A_k
  taylor: a_k = A_k/k!

Then propagate the balls through finite upper-triangular Toeplitz minors.
```

The Toeplitz matrix is:

```text
T[i,j] = a[j-i] if j >= i, otherwise 0.
```

## Source Cache

```text
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
dps = 90
n_sum = 100
cutoff = 8
```

## Probe Results

### Taylor Sequence, Small Clean Probe

Run:

```text
sequence = taylor = A_k/k!
lambda = 0
matrix size N = 6
orders <= 3
absolute radius on A_k = 1e-80
Arb dps = 120
```

Result:

```text
tests: 661
positive: 301
zero: 360
negative: 0
inconclusive_contains_zero: 0
```

This finite propagation gate is clean.

### Taylor Sequence, Wider Probe Before Structural-Zero Classification

Run:

```text
sequence = taylor = A_k/k!
lambda = 0
matrix size N = 10
orders <= 4
absolute radius on A_k = 1e-80
Arb dps = 140
```

Result:

```text
tests: 60,625
positive: 19,690
zero: 37,635
negative: 0
inconclusive_contains_zero: 3,300
```

The same count occurs with assumed radii `1e-90` and `1e-100`.

Interpretation:

The wider `taylor` probe finds no negative intervals. The 3,300 inconclusive rows line up with minors that are zero in the exact-rational point audit but become interval-wide under independent coefficient balls. This looks like interval dependency/overestimation around zero minors, not evidence of negative Toeplitz minors.

For comparison, the exact-rational point audit at the same range had:

```text
taylor, N = 10, orders <= 4:
  positive: 19,690
  zero: 40,935
  negative: 0
```

So:

```text
40,935 exact point zeros
-37,635 Arb exact zeros
= 3,300 dependency-sensitive zero minors
```

### Taylor Sequence, Wider Probe After Structural-Zero Classification

The updated probe classifies determinant-polynomial structural zeros before calling Arb determinant arithmetic.

For the same run:

```text
sequence = taylor = A_k/k!
lambda = 0
matrix size N = 10
orders <= 4
absolute radius on A_k = 1e-80
Arb dps = 140
```

the result is now:

```text
tests: 60,625
positive: 19,690
structural_zero: 40,935
negative: 0
inconclusive_contains_zero: 0
```

Additional structural Arb passes:

```text
taylor, lambda = 0, N = 10, orders <= 5:
  tests = 124,129
  positive = 39,094
  structural_zero = 85,035
  negative = 0
  inconclusive_contains_zero = 0

taylor, lambda = 0, N = 12, orders <= 4:
  tests = 297,925
  positive = 88,309
  structural_zero = 209,616
  negative = 0
  inconclusive_contains_zero = 0
```

For `lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}`, `N = 10`, orders `<= 4`, every structural Arb run has:

```text
positive = 19,690
structural_zero = 40,935
negative = 0
inconclusive_contains_zero = 0
```

### Beta Sequence Red Flag

Run:

```text
sequence = beta = A_k
lambda = 0
matrix size N = 10
orders <= 3
absolute radius on A_k = 1e-80
Arb dps = 140
```

Result:

```text
tests: 16,525
positive: 5,127
zero: 10,695
negative: 703
inconclusive_contains_zero: 0
```

This confirms that the Arb probe detects the known `A_k` Toeplitz failure cleanly.

## Meaning For The Proof Programme

The corrected ordinary Taylor sequence remains alive:

```text
c_k = A_k/k!
```

No finite Arb propagation run found a negative Toeplitz minor for `c_k`.

The structural-zero part of the certification hygiene is now handled for the tested ranges. The next blocker is coefficient enclosure:

```text
cached numerical coefficients + assumed coefficient radius
  -> not yet a proof that the balls enclose exact analytic coefficients.
```

## Next Compute Tasks

1. Replace assumed coefficient balls with interval-enclosed moment computation.
2. For nonzero minors, require Arb sign separation from proved coefficient balls.
3. Preserve the structural-zero classifier as a front gate for all Toeplitz/PF audits.
4. Consider determinant-polynomial evaluation methods only if nonstructural interval dependency appears in higher ranges.
