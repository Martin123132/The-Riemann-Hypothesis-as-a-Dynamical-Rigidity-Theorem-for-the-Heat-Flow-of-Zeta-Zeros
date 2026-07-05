# Coefficient Tail-Bound Status

Date: 2026-07-03

Status: analytic truncation/cutoff tail bounds. This is not yet a full coefficient enclosure, because quadrature error for the retained finite integral is not certified.

## Purpose

The structural Toeplitz/PF Arb probe currently assumes coefficient balls around cached `A_k(lambda)` values, often with absolute radius:

```text
1e-80
```

To turn that into a real finite certificate, we need to justify the coefficient balls. This pass handles two components:

```text
1. omitted n-series tail, n > n_sum, on 0 <= u <= cutoff;
2. omitted u-integral tail, u > cutoff, with the infinite n-series.
```

It does not yet handle:

```text
3. quadrature error for the retained finite n_sum terms on 0 <= u <= cutoff.
```

## Script

```text
work/rh_compute/scripts/coefficient_tail_bounds.py
```

The script bounds:

```text
mu_2k(lambda) = 2 integral_0^infty u^(2k) exp(lambda u^2) Phi(u) du
```

and converts the resulting bound into:

```text
A_k(lambda) = mu_2k(lambda) k! / (2k)!
c_k(lambda) = A_k(lambda) / k! = mu_2k(lambda) / (2k)!.
```

## Bound Shape

For the omitted `n > n_sum` series on `[0, C]`, the script uses:

```text
exp(-pi n^2 exp(4u)) <= exp(-pi n^2)
```

plus a Gaussian integral bound for:

```text
sum_{n > n_sum} n^p exp(-pi n^2).
```

For the omitted `u > C` tail, it uses:

```text
sum_{n >= 1} n^p exp(-pi n^2 exp(4u))
  <= K_p exp(-pi exp(4u)),
```

where:

```text
K_p = sum_{n >= 1} n^p exp(-pi(n^2 - 1)).
```

Then it applies an endpoint bound to:

```text
integral_C^infty u^(2k) exp(lambda u^2 + a u - pi exp(4u)) du.
```

The monotonicity gate passed in all reported runs.

## Existing Hankel/Toeplitz Cache Settings

Run:

```text
coefficient_tail_bounds_repro_hankel_15c
lambdas = {0, 1e-6, 1e-4, 1e-2, 1e-1}
max_k = 32
n_sum = 100
cutoff = 8
```

Summary:

```text
worst mu_2k tail:
  lambda = 0.1
  k = 32
  bound about 9.36e-13547

worst A_k tail:
  lambda = 0.1
  k = 15
  bound about 1.91e-13597

worst c_k tail:
  lambda = 0.1
  k = 3
  bound about 5.04e-13601
```

Comparison:

```text
target Arb radius: 1e-80
tail bounds:       < 1e-13596 for A_k
```

So the omitted series and cutoff tails are negligible compared with the Arb radius used in the Toeplitz propagation probes.

## Jensen Fragile-Region Cache Settings

Run:

```text
coefficient_tail_bounds_jensen_d16_d20_n50_n60
lambda = 0
max_k = 80
n_sum = 120
cutoff = 9
```

Summary:

```text
worst mu_2k tail:
  k = 80
  bound about 4.17e-19454

worst A_k tail:
  k = 19
  bound about 1.53e-19597

worst c_k tail:
  k = 4
  bound about 1.67e-19602
```

Again, the omitted tails are far below `1e-80`.

## What This Certifies

Conditional on the retained finite integral values being accurate to the cached precision, the already-used `1e-80` coefficient radius comfortably covers:

```text
decimal output rounding at roughly 80 printed digits;
omitted n-series tail;
omitted u-cutoff tail.
```

## What It Does Not Certify Yet

By itself, this tail-bound script does not prove that the cached coefficients enclose the exact coefficients, because the original retained integral was computed by `mpmath.quad` without a rigorous quadrature error certificate.

The remaining coefficient-enclosure gate is:

```text
finite sum over n = 1..n_sum
finite interval u in [0, cutoff]
validated quadrature / interval integration
```

Update:

That gate now has a selected-range implementation:

```text
work/rh_compute/scripts/acb_coefficient_enclosures.py
```

It uses `python-flint`'s rigorous `acb.integral` for the retained finite integral, then adds the analytic tail bounds from this file.

Certified coefficient-enclosure runs now cover:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}, k = 0..32,
n_sum = 100, cutoff = 8
```

In these runs, the existing cache balls `A_k(cache) +/- 1e-80` contain the rigorous `A_k` balls. See `outputs/finite_toeplitz_certificate_status.md`.

Once that is done, the pipeline becomes:

```text
validated finite integral enclosure
  + analytic n-tail and cutoff-tail bounds
  + decimal rounding bound
  -> proved coefficient balls
  -> structural-zero-aware Arb Toeplitz/PF probe
  -> certified finite Toeplitz minor signs.
```

## Immediate Next Options

1. Extend validated `acb.integral` coefficient enclosures beyond the currently certified ranges.
2. Feed those proved coefficient balls into larger structural Toeplitz/PF probes.
3. Keep treating `1e-80` as a propagation-only radius for any coefficient range not yet covered by `acb_coefficient_enclosures.py`.
4. For larger `k`, investigate a transformed integral in `x = exp(4u)` if direct `acb.integral` becomes slow or loses accuracy.
