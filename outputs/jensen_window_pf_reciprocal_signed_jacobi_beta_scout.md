# Jensen-Window PF Reciprocal Signed Jacobi Beta Scout

Date: 2026-07-06

Status: signed Jacobi beta-parameter diagnostic. This is not a proof of a signed continued-fraction theorem,
a production-matrix model, the all-order column recurrence, Schur
positivity, Jensen-window PF-infinity, Jensen hyperbolicity, RH, or
`Lambda <= 0`.

## Purpose

The signed J-fraction scout established a finite pattern for the reciprocal
Hankel determinants and the off-diagonal Jacobi parameters. This note adds the
diagonal Jacobi parameters for:

```text
E(t) = 1/H(-t)
mu_m = [t^m] E(t)
```

For:

```text
Delta_r = det(mu_{i+j})_{0<=i,j<r}
```

the diagonal parameter is not the fully shifted determinant
`det(mu_{i+j+1})`. Instead define `Delta_r^*` by shifting only the final
Hankel column:

```text
Delta_r^* = det(mu_{i+j}) with the final column replaced by mu_{i+r}
Q_r = Delta_r^*/Delta_r
Q_0 = 0
beta_n = Q_{n+1} - Q_n
```

This matters already in degree 2: for `H(t)=1+g_1*t+g_2*t^2`, the formula gives
`beta_0=g_1`, `lambda_1=-g_2`, and terminal `beta_1=0`.

## Outputs

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py
```

Current result:

```text
validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues
```

## Finite Grid

The Arb grid checks:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
degrees d=2..8
dps=520
```

It validates the finite signature:

```text
beta_0 > 0
beta_1 < 0 for d >= 3
beta_n > 0 for n >= 2
terminal degree-2 beta_1 contains zero
```

Counts:

```text
2,940 / 3,675 beta rows positive
630 / 3,675 beta rows negative
105 / 3,675 terminal degree-2 rows contain zero
```

This sharpens:

```text
rp_09_signed_or_modified_continued_fraction
```

and makes the candidate signed-Jacobi model more concrete. It is finite evidence only.

## Interpretation

The ordinary positive J-fraction route is still blocked, but the signed route
now has a more precise finite target. Any viable signed Jacobi, oscillatory
matrix, or production-matrix theorem must explain not only:

```text
lambda_n < 0
kappa_n = -lambda_n > 0
```

but also the diagonal pattern:

```text
beta_0 > 0, beta_1 < 0, beta_n > 0 for n >= 2
```

with the degree-2 terminal zero handled explicitly.

This does not prove coefficientwise positivity of `E(t)=1/H(-t)`. It narrows
the missing theorem and adds a formula guard against using the wrong shifted
Hankel determinant.

## Boundary

Passing this checker means the finite signed Jacobi beta-parameter pattern is
reproducible on the stated grid. It does not prove an all-order beta
signature, a signed continued-fraction theorem, a production-matrix model, the
column recurrence for all windows, Schur positivity, Jensen-window
PF-infinity, or `Lambda <= 0`.
