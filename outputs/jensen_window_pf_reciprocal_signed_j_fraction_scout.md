# Jensen-Window PF Reciprocal Signed J-Fraction Scout

Date: 2026-07-06

Status: signed J-fraction Hankel-signature diagnostic. This is not a proof of a signed continued-fraction theorem,
the all-order column recurrence, Schur
positivity, Jensen-window PF-infinity, Jensen hyperbolicity, RH, or
`Lambda <= 0`.

## Purpose

The standard positive Stieltjes/Jacobi fraction route is blocked for:

```text
E(t) = 1/H(-t)
```

because the first ordinary quadratic parameter is negative. The remaining
fraction route is therefore a signed or modified J-fraction. This scout tests
the determinant signature that such a route would naturally require.

For:

```text
mu_m = [t^m] E(t)
Delta_r = det(mu_{i+j})_{0<=i,j<r}
```

the finite signed-Hankel target is:

```text
(-1)^(r(r-1)/2) Delta_r > 0
```

The ordinary Jacobi parameter is:

```text
lambda_n = Delta_{n+1} Delta_{n-1} / Delta_n^2
```

So this signature implies:

```text
lambda_n < 0
kappa_n = -lambda_n > 0
```

because:

```text
n(n+1)/2 + (n-1)(n-2)/2 = n^2 - n + 1
```

is always odd.

## Outputs

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_reciprocal_signed_j_fraction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py
```

Current result:

```text
validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues
```

## Finite Grid

The Arb grid checks:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
degrees d=2..8
dps=520
```

It validates:

```text
3,675 / 3,675 signed reciprocal-Hankel determinant rows positive
2,940 / 2,940 ordinary Jacobi lambda rows negative
2,940 / 2,940 signed kappa_n=-lambda_n rows positive
```

This is evidence for:

```text
rp_09_signed_or_modified_continued_fraction
```

and a sharper target for a signed J-fraction theorem. It is finite evidence only.

## Interpretation

The finite result says the signed-fraction route is not dead in the same way
the standard positive fraction route is dead. Instead, the remaining burden is
now very specific:

```text
find or prove a signed/oscillatory/production-matrix theorem where
(-1)^(r(r-1)/2) Delta_r > 0
or kappa_n=-lambda_n>0
implies coefficient positivity of E(t)=1/H(-t)
for the actual zeta heat-flow Jensen windows
```

This still would prove only the column recurrence unless it is lifted to a
positive specialization, planar network, or all-Schur/Toeplitz theorem.

## Boundary

Passing this checker means the finite signed reciprocal-Hankel and signed
J-fraction parameter pattern is reproducible on the stated grid. It does not
prove the all-order determinant signature, a signed continued-fraction theorem,
the column recurrence for all windows, Schur positivity, Jensen-window
PF-infinity, or `Lambda <= 0`.
