# Jensen/Hankel Runner Smoke Report

Date: 2026-07-03

Status: smoke-test report. This is not a proof artifact; it records reproducibility checks for the Jensen/Hankel runner.

Runner:

```text
work/rh_compute/scripts/jensen_hankel_runner.py
```

Purpose:

Verify that the first clean command-line runner can compute moments, build `A_k(lambda)`, test signed Hankel determinants, count positive roots for small Jensen cases, and write JSONL/CSV outputs.

## Smoke 1: Signed Hankel

Command shape:

```text
mode = hankel
lambda = 0
dps = 50
n_sum = 30
cutoff = 6
m = 0..2
shift = 0..1
```

Result:

```text
6/6 passed
failures: 0
```

Output files:

```text
work/rh_compute/results/smoke_hankel.jsonl
work/rh_compute/results/smoke_hankel.csv
work/rh_compute/results/smoke_hankel_summary.json
```

## Smoke 2: Jensen

Command shape:

```text
mode = jensen
lambda = 0
dps = 50
rational_digits = 35
n_sum = 30
cutoff = 6
d = 2..3
n = 0..2
```

Result:

```text
6/6 passed
failures: 0
```

Output files:

```text
work/rh_compute/results/smoke_jensen.jsonl
work/rh_compute/results/smoke_jensen.csv
work/rh_compute/results/smoke_jensen_summary.json
```

## Important Limitation

The Jensen root count currently uses exact Sturm counting after rationalizing high-precision numerical coefficients. This is useful for reproducible testing, but it is not yet a rigorous interval certificate for the original transcendental moments.

Rigorous certification still requires coefficient error bounds, interval arithmetic, or Arb/FLINT-style ball arithmetic.
