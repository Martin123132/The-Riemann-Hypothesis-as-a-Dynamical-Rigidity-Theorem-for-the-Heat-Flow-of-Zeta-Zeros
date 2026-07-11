# Jensen-Window PF Monotone-Contraction All-M Counterexample

Date: 2026-07-10

Status: exact countermodel gate. This is not evidence against RH,
not evidence against `Lambda <= 0`, and not a zeta-window computation.

Artifact kind: `jensen_window_pf_monotone_contraction_all_m_counterexample`.

Proof boundary: this artifact blocks one generic proof step only:
the full static ratio cone does not imply all-`m` column recurrence
positivity. A zeta-specific strengthened theorem may still exist.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_all_m_counterexample.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_all_m_counterexample.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_all_m_counterexample.py
```

Current result:

```text
validated Jensen-window PF monotone-contraction all-m counterexample: degree 7, m=11, exact full-cone witness, 6 lower walls, negative normalized value, 0 issues
```

## Exact Witness

For degree `d=7` and column row `m=11`, after factoring the positive
monomial `A^m*rho^m`, evaluate the normalized recurrence at:

```text
s = (19/20, 0, 1, 0, 0, 0)
x = (19/20, 19/20, 1, 1, 1, 1)
```

The full static-cone checks are:

```text
all_in_unit_interval = True
weakly_increasing = True
upper_bound_ok = True
lower_walls = (1/3, 3/5, 5/7, 7/9, 9/11, 11/13)
lower_wall_margins = (37/60, 7/20, 2/7, 2/9, 2/11, 2/13)
all_lower_walls_ok = True
tail = x_k=1 for every k>=7
tail_lower_wall_ok = True
```

The exact normalized column value is:

```text
Q_11 = -2611049410212561054670871/163840000000000000000
Q_11 approximately -15936.580872879402
```

This value is strictly negative.

## Consequence

The exact shift-0 witness x=(19/20,19/20,1,1,1,1), extended by x_k=1 for k>=7, satisfies the full propagated static ratio cone, including every pointwise lower wall, but makes the normalized degree-7 m=11 column recurrence negative. Thus the static cone is not an all-order PF theorem by itself.

Integration:

```text
outputs/jensen_window_pf_column_recurrence_contract.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The exact shift-0 witness x=(19/20,19/20,1,1,1,1), extended by x_k=1 for k>=7, satisfies the full propagated static ratio cone, including every pointwise lower wall, but makes the normalized degree-7 m=11 column recurrence negative. Thus the static cone is not an all-order PF theorem by itself.
