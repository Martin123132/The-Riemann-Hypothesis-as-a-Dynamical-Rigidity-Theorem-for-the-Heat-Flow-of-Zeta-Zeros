# Jensen-Window PF Column Recurrence Finite Coverage

Date: 2026-07-06

Status: finite coverage map. This is not a proof of the all-order column
recurrence theorem, Jensen-window PF-infinity, Jensen hyperbolicity,
Laguerre-Polya membership, RH, or `Lambda <= 0`; it maps existing finite
evidence onto the column recurrence contract. It is finite evidence only and
not an all-order recurrence theorem.

## Purpose

The column recurrence contract isolates the necessary condition:

```text
C[m] >= 0
for every m >= 0, every degree d, and every shift n.
```

This note records which already-promoted finite diagnostics support that
condition for checked zeta heat-flow windows.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_column_recurrence_finite_coverage.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py
```

Current result:

```text
validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues
```

## Direct Arb Determinant Coverage

The direct Arb determinant manifest checks the column/contiguous determinants:

```text
d=3, m=1..8
d=4, m=1..6
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
```

Source files:

```text
work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520_summary.json
work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520.jsonl
python work/rh_compute/scripts/check_arb_jensen_window_pf_obligation_manifest.py
```

Coverage:

```text
direct checked rows = 1470
direct positive rows = 1470
failed or inconclusive rows = 0
```

The two hard recurrence sizes are included:

```text
d=3,m=8:
  checked rows = 105
  positive rows = 105

d=4,m=6:
  checked rows = 105
  positive rows = 105
```

So the exact rational countermodel kills generic bridge conditions, but the
actual checked zeta-window grid is positive on these hard column recurrence
rows.

## Sturm/PF Window Coverage

The finite Sturm-to-PF consequence supplies a different finite route:

```text
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json
work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json
outputs/jensen_window_sturm_pf_consequence.md
python work/rh_compute/scripts/check_jensen_window_sturm_pf_consequence.py
```

For each checked window, Sturm-certified finite PF-infinity implies all finite
Toeplitz minors of that one binomial window are nonnegative, including all
column recurrence minors.

Coverage:

```text
degree d = 3,4,5
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
shifts n=0..20
315 checked Sturm/PF windows
ok Sturm/PF windows = 315
```

This is a finite window-by-window consequence, not an all-order recurrence
theorem.

## Recurrence Stress Grid

A wider recurrence-only Arb stress diagnostic is:

```text
outputs/arb_jensen_window_column_recurrence_stress.md
work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json
work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520.jsonl
python work/rh_compute/scripts/arb_jensen_window_column_recurrence_stress.py
python work/rh_compute/scripts/check_arb_jensen_window_column_recurrence_stress.py
```

It validates:

```text
validated 12600 Arb Jensen-window column recurrence stress rows with 0 issues
```

This checks degrees `d=3..8`, sizes `m=1..20`, five lambdas, and shifts
`n=0..20`. It is finite necessary-condition evidence, not an all-order
recurrence theorem.

## Remaining Gap

Covered:

```text
direct Arb determinant positivity for d=3, m=1..8 and d=4, m=1..6
direct Arb determinant positivity for d=3,m=8 and d=4,m=6 hard recurrence sizes
finite Sturm/PF consequence for d=3,4,5 windows on the checked grid
```

Not covered:

```text
all degrees d
all shifts n
all lambda values or an interval in lambda
an analytic proof that the recurrence is nonnegative for the actual zeta windows
all non-column skew shapes beyond the finite checked windows
```

The current reciprocal theorem-mechanism matrix is:

```text
outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md
work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json
python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py
```

It keeps the 12,600-row stress grid as finite evidence only, while separating
the live positive renewal, continued-fraction, and production-matrix routes
from endpoint/circular and countermodel-rejected shortcuts.

## Boundary

Passing this checker means the finite evidence is correctly mapped onto the
column recurrence contract. It does not prove the all-order recurrence, Schur
positivity, a positive specialization, a planar network, a production matrix,
a determinant integral, or an all-order sign-regularity-to-Jensen-window PF
theorem.
