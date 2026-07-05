# Progress Report: Toeplitz Normalization Correction

Date: 2026-07-03

Status: normalization progress report. This is not a proof artifact; it records a correction to the coefficient-PF diagnostic route.

## Summary

The prior conclusion that the "simple PF-infinity Toeplitz route is damaged" was too broad.

The negative Toeplitz minors are real for the tested exponential/Jensen coefficients:

```text
beta_k = A_k.
```

But the classical ASW/Edrei coefficient theorem for the restricted Laguerre-Polya class uses ordinary Taylor coefficients:

```text
c_k = A_k / k!.
```

After correcting the normalization, the finite Toeplitz diagnostics become sign-clean in the ranges tested.

## New Script

```text
work/rh_compute/scripts/toeplitz_pf_audit.py
```

Purpose:

```text
Build upper-triangular Toeplitz matrices T[i,j] = a[j-i] for j >= i,
rationalize cached coefficients, and compute finite minors exactly over QQ.
```

This is still not a rigorous analytic certificate, because the input coefficients come from high-precision numerical moment computation.

Additional propagation script:

```text
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

This uses Arb balls centered on cached coefficients with an assumed absolute radius. It is a propagation gate, not a proof of coefficient enclosure.

Structural-zero update:

```text
work/rh_compute/scripts/toeplitz_pf_audit.py
work/rh_compute/scripts/arb_toeplitz_pf_probe.py
```

now classify Toeplitz minors whose determinant polynomial is identically zero before evaluating the zeta coefficients. This removes the earlier zero-minor interval overestimation from the Arb probe.

## Main Result

Source cache:

```text
work/rh_compute/results/repro_hankel_15c_coefficients.jsonl
dps = 90
n_sum = 100
cutoff = 8
```

At `lambda = 0`, `N = 10`, order `<= 4`:

```text
beta = A_k:
  tests = 60,625
  negatives = 4,589
  first known negative:
    rows = 0,1,2
    cols = 1,2,3
    det = -1.71119821042992898147325e-11

taylor = A_k/k!:
  tests = 60,625
  negatives = 0
  same minor:
    det = +7.59800565413986900638776e-11
```

Additional `taylor` checks:

```text
lambda = 0, N = 10, order <= 5:
  124,129 minors, 0 negative

lambda = 0, N = 12, order <= 4:
  297,925 minors, 0 negative

lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}, N = 10, order <= 4:
  60,625 minors per lambda, 0 negative
```

## Interpretation

This revives a classical coefficient route:

```text
prove all Toeplitz minors of c_k = A_k/k! are nonnegative
  -> c_k is a PF-infinity sequence
  -> F_0(s) = sum c_k s^k is restricted Laguerre-Polya
  -> Xi has only real zeros
  -> RH.
```

This is not yet a proof route with a missing calculation; it is a correctly normalized theorem-search target.

## Hard Caveats

- Finite minors do not prove PF-infinity.
- Rationalized numerical coefficients do not certify exact analytic signs.
- Some positive minors are extremely small; in the `lambda=0`, `N=10`, order `<=4` `taylor` run, the smallest nonzero minor was about `9.69e-75`.
- Arb propagation now finds no negative or inconclusive intervals for `A_k/k!` in the tested structural-zero-classified ranges, but the coefficient balls are still assumed rather than analytically certified.
- Schoenberg's direct zeta condition concerns an inverse transform of `1/Xi` or `1/Xi_1`, not the original heat kernel moments.

## Next Work

1. Add interval-enclosed moment computation for `c_k = mu_2k/(2k)!`.
2. Feed proved coefficient balls into the structural-zero-classified Arb Toeplitz probe.
3. Search ASW/Edrei/Karlin literature for structural sufficient conditions beyond brute-force minors.
4. Compare the corrected PF route with the existing signed-Hankel/Jensen route and keep them separate until a theorem connects them.
