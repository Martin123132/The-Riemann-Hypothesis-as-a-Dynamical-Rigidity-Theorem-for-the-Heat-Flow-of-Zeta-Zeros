# Jensen-Window PF Reciprocal Coefficient Extended Stress

Date: 2026-07-06

Status: extended finite Arb reciprocal-coefficient diagnostic. This is not a proof of the all-order
column recurrence theorem, Schur positivity, Jensen-window PF-infinity,
Jensen hyperbolicity, RH, or `Lambda <= 0`.

## Purpose

The reciprocal route matrix reduces the column recurrence to positivity of:

```text
E(t) = 1 / H(-t)
H(t) = 1 + g_1*t + ... + g_d*t^d
g_j = binom(d,j) A_{n+j}(lambda) / A_n(lambda)
```

The target is:

```text
[t^m] E(t) >= 0
```

for every `m`, degree `d`, and shift `n`. This diagnostic tests a larger
finite grid directly in normalized reciprocal coefficients.

## Outputs

Machine-readable outputs:

```text
work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json
work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_lamgrid_n0_n32_d2_d12_m1_m40_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_reciprocal_coefficient_extended_stress.py
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_coefficient_extended_stress.py
```

Current result:

```text
validated Jensen-window PF reciprocal coefficient extended stress: 72600 rows, 0 issues
```

## Checked Grid

```text
degrees d=2..12
sizes m=1..40
shifts n=0..32
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
dps = 520
needed max coefficient index = 44
```

This gives:

```text
summary rows = 1,815
stress rows = 72,600
positive rows = 72,600
failed or inconclusive rows = 0
72,600 / 72,600 checked reciprocal coefficients are positive and separated from zero
```

Per degree:

```text
d=2: 6,600 positive rows
d=3: 6,600 positive rows
d=4: 6,600 positive rows
d=5: 6,600 positive rows
d=6: 6,600 positive rows
d=7: 6,600 positive rows
d=8: 6,600 positive rows
d=9: 6,600 positive rows
d=10: 6,600 positive rows
d=11: 6,600 positive rows
d=12: 6,600 positive rows
```

## Interpretation

This is finite evidence only. It gives wider support for the exact
reciprocal coefficient target than the earlier recurrence stress grid:

```text
old stress: d=3..8, n=0..20, m=1..20, 12,600 rows
extended stress: d=2..12, n=0..32, m=1..40, 72,600 rows
```

The result supports the signed/modified continued-fraction, positive renewal,
and production-matrix theorem search, but it does not identify the missing
all-order mechanism.

## Boundary

Passing this checker means the stated finite Arb grid has positive normalized
reciprocal coefficients. It is finite evidence only and not an all-order recurrence theorem.
It does not prove Schur positivity, a positive specialization, a planar
network, a production matrix, a determinant integral, or `Lambda <= 0`.
