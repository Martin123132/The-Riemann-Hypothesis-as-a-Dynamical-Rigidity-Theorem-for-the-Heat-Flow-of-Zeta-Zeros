# Jensen-Window PF Monotone-Contraction Column Extension Scout

Date: 2026-07-06

Status: exact bounded column-extension diagnostic. This is not a proof
of Jensen-window PF-infinity, all-shape Schur positivity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_column_extension_scout`.

Proof boundary: this artifact proves finite exact Bernstein certificates
for bounded column rows under monotone contractions. It does not prove
that the zeta coefficients satisfy monotone contractions, and it does
not prove all column rows or all Schur shapes.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_column_extension_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_column_extension_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_column_extension_scout.py
```

Current result:

```text
validated Jensen-window PF monotone-contraction column extension scout: 25 column rows, 3329 Bernstein coefficients, 3 beyond-frontier rows, 0 negative Bernstein rows, 0 issues
```

## Certified Region

The scout works on the monotone-contraction cones:

```text
degree 3: 0 <= x1 <= x2 <= 1
degree 4: 0 <= x1 <= x2 <= x3 <= 1
degree 5: 0 <= x1 <= x2 <= x3 <= x4 <= 1
```

where `x_i=(A_{n+i+1}/A_{n+i})/(A_{n+i}/A_{n+i-1})`.

## Column Rows

```text
mccx_d3_m1: degree 3, m=1, Bernstein count=1, min=3, beyond_frontier=False
mccx_d3_m2: degree 3, m=2, Bernstein count=2, min=6, beyond_frontier=False
mccx_d3_m3: degree 3, m=3, Bernstein count=8, min=10, beyond_frontier=False
mccx_d3_m4: degree 3, m=4, Bernstein count=8, min=15, beyond_frontier=False
mccx_d3_m5: degree 3, m=5, Bernstein count=10, min=21, beyond_frontier=False
mccx_d3_m6: degree 3, m=6, Bernstein count=21, min=28, beyond_frontier=False
mccx_d3_m7: degree 3, m=7, Bernstein count=21, min=36, beyond_frontier=False
mccx_d3_m8: degree 3, m=8, Bernstein count=24, min=45, beyond_frontier=False
mccx_d3_m9: degree 3, m=9, Bernstein count=40, min=55, beyond_frontier=True
mccx_d3_m10: degree 3, m=10, Bernstein count=40, min=66, beyond_frontier=True
mccx_d4_m1: degree 4, m=1, Bernstein count=1, min=4, beyond_frontier=False
mccx_d4_m2: degree 4, m=2, Bernstein count=2, min=10, beyond_frontier=False
mccx_d4_m3: degree 4, m=3, Bernstein count=8, min=20, beyond_frontier=False
mccx_d4_m4: degree 4, m=4, Bernstein count=56, min=35, beyond_frontier=False
mccx_d4_m5: degree 4, m=5, Bernstein count=56, min=56, beyond_frontier=False
mccx_d4_m6: degree 4, m=6, Bernstein count=64, min=84, beyond_frontier=False
mccx_d4_m7: degree 4, m=7, Bernstein count=100, min=120, beyond_frontier=True
mccx_d5_m1: degree 5, m=1, Bernstein count=1, min=5, beyond_frontier=False
mccx_d5_m2: degree 5, m=2, Bernstein count=2, min=15, beyond_frontier=False
mccx_d5_m3: degree 5, m=3, Bernstein count=8, min=35, beyond_frontier=False
mccx_d5_m4: degree 5, m=4, Bernstein count=56, min=70, beyond_frontier=False
mccx_d5_m5: degree 5, m=5, Bernstein count=616, min=126, beyond_frontier=False
mccx_d5_m6: degree 5, m=6, Bernstein count=616, min=210, beyond_frontier=False
mccx_d5_m7: degree 5, m=7, Bernstein count=672, min=330, beyond_frontier=False
mccx_d5_m8: degree 5, m=8, Bernstein count=896, min=495, beyond_frontier=False
```

The extension beyond the original first hard frontier is:

```text
degree 3: m=9, m=10
degree 4: m=7
```

The higher-degree extension is:

```text
degree 5: m=1..8
```

## Consequence

This strengthens the monotone-contraction route from a two-row escape
from the rational countermodel into a bounded column-family theorem
search target that now includes a first degree-5 band. It still leaves two hard problems open: proving
monotone contractions for the actual zeta heat-flow coefficients and
lifting from bounded column rows to all column rows/all Schur shapes.

Integration:

```text
outputs/jensen_window_pf_monotone_contraction_frontier_scout.md
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
```

Summary:

The monotone-contraction cone certifies more than the first two hard column-frontier rows: degree 3 is certified through m=10 and degree 4 through m=7, with three rows beyond the original hard frontier; the same method also certifies degree 5 through m=8. This strengthens the monotone-contraction route as a column-family target, while leaving the all-m/all-shape and zeta cone-entry theorems open.
