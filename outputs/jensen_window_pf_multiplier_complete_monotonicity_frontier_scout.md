# Jensen-Window PF Multiplier Complete-Monotonicity Frontier Scout

Date: 2026-07-11

Status: finite high-precision interval frontier diagnostic. This is not a
proof of all-k complete monotonicity, a counting-measure factorization,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_multiplier_complete_monotonicity_frontier_scout`.

```text
work/rh_compute/results/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.json
python work/rh_compute/scripts/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py
```

Current result:

```text
validated Jensen-window PF multiplier complete-monotonicity frontier scout: 7980 positive intervals, 0 inconclusive, orders 0..55, 5 lambdas, 0 issues
```

## Certified Frontier

Using the dps220 coefficient enclosures `A_0..A_57` and evaluating all
derived Arb operations at 250 decimal digits proves

```text
(-1)^m*Delta^m y_k>0,
y_k=-log(x_k),
0<=m<=55, 1<=k<=56-m,
lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}.
```

All `7980` interval signs are strict. The weakest
certified row is

```text
lambda=0.0, order=29, k=27
ball=[1.061559125732763522638512075777491877176364411342642089348417289925890346466165e-19 +/- 4.94e-98]
```

The old order-9 frontier came from evaluating derived Arb arithmetic at
insufficient working precision; the tighter source balls themselves already
support the complete finite order-55 triangle.

## Consequence

This removes finite complete monotonicity as an immediate falsification of the
multiplier counting-measure ansatz. It does not construct the required unit-atomic
measure. The next discriminating obligation is the integer-to-continuous
interpolation/uniqueness step or another unit-multiplicity constraint stronger
than ordinary Hausdorff complete monotonicity.

```text
outputs/jensen_window_pf_defect_complete_monotonicity_scout.md
outputs/jensen_window_pf_multiplier_counting_measure_target.md
outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md
outputs/signed_hankel_jensen_dependency_graph.md
```
