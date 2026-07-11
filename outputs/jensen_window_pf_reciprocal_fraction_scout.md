# Jensen-Window PF Reciprocal Fraction Scout

Date: 2026-07-06

Status: reciprocal continued-fraction sign diagnostic. This is not a proof of the all-order column recurrence,
Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The reciprocal route matrix left a possible continued-fraction route:

```text
rp_03_positive_stieltjes_or_j_fraction
rp_04_companion_or_production_matrix_total_positivity
```

This scout checks the first formal parameters for:

```text
E(t) = 1 / H(-t)
H(t) = 1 + g_1*t + g_2*t^2 + ...
```

Equivalently:

```text
E(t) = 1 / (1 - g_1*t + g_2*t^2 - ...)
```

## Result

For the standard S-fraction convention:

```text
F(t)=1/(1-a_1*t/(1-a_2*t/(1-...)))
```

the first two parameters are:

```text
a_1 = g_1
a_2 = -g_2/g_1
```

For the standard J-fraction convention:

```text
F(t)=1/(1-beta_0*t-lambda_1*t^2*F_1(t))
F_1(0)=1
```

the first parameters are:

```text
beta_0 = g_1
lambda_1 = -g_2
```

So when `g_1>0` and `g_2>0`, the ordinary positive Stieltjes/Jacobi continued
fraction route has the wrong first nontrivial sign. This does not rule out a
signed or modified continued fraction, but it does reject the standard
positive S-fraction/J-fraction route for `E(t)=1/H(-t)`.

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json
work/rh_compute/results/jensen_window_pf_reciprocal_fraction_sign_lamgrid_n0_n20_d2_d8_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_reciprocal_fraction_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py
```

Current result:

```text
validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues
```

## Finite Sign Check

The Arb grid checks:

```text
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
degrees d=2..8
dps=520
```

All `735` checked finite rows have:

```text
h_0 > 0
h_1 > 0
h_2 > 0
standard S-fraction a_1 > 0
standard S-fraction a_2 < 0
standard J-fraction beta_0 > 0
standard J-fraction lambda_1 < 0
```

This finite sign check supports the symbolic obstruction on the current
zeta-window grid. It is not an all-degree or all-shift theorem.

## Consequence

The continued-fraction route is narrower now:

```text
rejected:
  standard positive Stieltjes/Jacobi S-fraction or J-fraction for E(t)

still possible:
  a signed or modified continued fraction with its own total-positivity,
  oscillatory-matrix, or production-matrix theorem
```

Any future fraction proof must explain why the negative ordinary
`a_2`/`lambda_1` parameter is harmless under a different theorem, rather than
silently citing the standard positive Stieltjes/Jacobi machinery.

## Boundary

Passing this checker means the first continued-fraction sign obstruction is
reproducible and the finite zeta-window grid lies in that sign regime. It does
not prove or disprove a signed modified continued-fraction theorem, a positive
production matrix, the all-order column recurrence, Schur positivity, or
`Lambda <= 0`.
