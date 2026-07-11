# Arb Jensen-Window Column Recurrence Stress

Date: 2026-07-06

Status: finite Arb stress diagnostic. This is not a proof of the all-order
column recurrence theorem, Jensen-window PF-infinity, Jensen hyperbolicity,
Laguerre-Polya membership, RH, or `Lambda <= 0`; it stress-tests the necessary
column recurrence condition on a wider finite grid.

## Purpose

The column recurrence contract requires:

```text
C[m] >= 0
for every m >= 0, every degree d, and every shift n.
```

The earlier finite coverage map records direct determinant evidence for
`d=3,4` and Sturm/PF consequence evidence for `d=3,4,5`. This stress diagnostic
uses the recurrence itself to test a larger finite grid cheaply.

Machine-readable results:

```text
work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json
work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520.jsonl
```

Generator and checker:

```text
python work/rh_compute/scripts/arb_jensen_window_column_recurrence_stress.py
python work/rh_compute/scripts/check_arb_jensen_window_column_recurrence_stress.py
```

Current result:

```text
validated 12600 Arb Jensen-window column recurrence stress rows with 0 issues
```

## Checked Grid

```text
degrees d=3..8
sizes m=1..20
shifts n=0..20
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}
dps = 520
needed max coefficient index = 28
```

This gives:

```text
summary rows = 630
stress rows = 12600
positive rows = 12600
failed or inconclusive rows = 0
```

Per degree:

```text
d=3: 2100 positive rows
d=4: 2100 positive rows
d=5: 2100 positive rows
d=6: 2100 positive rows
d=7: 2100 positive rows
d=8: 2100 positive rows
```

## Interpretation

This is finite necessary-condition evidence for the actual checked zeta
heat-flow coefficient enclosures. It gives broader support for the column
recurrence target than the direct determinant manifest alone, especially for
degrees `d=6,7,8`.

It does not test all skew shapes, does not prove Schur positivity, and is not
an all-order recurrence theorem. In short: not an all-order recurrence theorem.

## Integration Points

This diagnostic supports:

```text
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_column_recurrence_finite_coverage.md
work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json
work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py
python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py
```

## Boundary

Passing this checker means the finite Arb recurrence stress grid is
reproducible and positive. It does not prove the all-order recurrence, Schur
positivity, a positive specialization, a planar network, a production matrix,
a determinant integral, or an all-order sign-regularity-to-Jensen-window PF
theorem.
