# Jensen-Window PF Negative-Lambda Scaled-Defect Frontier Scout

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
defect-tail theorem, cone entry, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_scaled_defect_frontier_scout`.

Proof boundary: this artifact separates finite scaled-defect barriers.
It does not prove an all-`k` tail theorem, does not prove `jwpf_06`,
and does not establish `Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k60.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.json --note outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.md
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 177 scaled rows, 177 cone rows, 177 half-width rows, 139 one-third rows, 38 one-third failures, 174 scaled-increase rows, 0 issues
```

## Scaled Defect

Write:

```text
d_k = 1 - x_k
s_k = ((2*k+1)/2) * d_k
```

The exact lower wall is `0 <= s_k <= 1`. The previously useful
one-third-width buffer is `s_k <= 1/3`; the weaker half-width buffer
is `s_k <= 1/2`.

## Frontier

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_60
checked contractions: x_1..x_59
exact cone rows: 177 / 177
half-width rows: 177 / 177
one-third rows: 139 / 177
one-third failure rows: 38
scaled-defect increase rows: 174 / 174
```

Per-lambda one-third frontier:

```text
lambda=-100.0: one-third 59/59, first failure none, max s=2.539780558950884235E-1 at k=59
lambda=-50.0: one-third 50/59, first failure k=51, margin=-2.719207370662713847E-4, max s=3.533907112847882969E-1 at k=59
lambda=-25.0: one-third 30/59, first failure k=31, margin=-1.768768251686068314E-3, max s=4.124052882015568604E-1 at k=59
```

Global extrema:

```text
min exact-cone scaled margin: 1.346067912947116963E-2 at lambda=-100.0, k=1
min half-width margin: 8.759471179844313957E-2 at lambda=-25.0, k=59
min one-third margin: -7.907195486822352710E-2 at lambda=-25.0, k=59
max scaled defect: 4.124052882015568604E-1 at lambda=-25.0, k=59
min scaled-defect increase: 1.988804647851078841E-3 at lambda=-25.0, k=58
```

## Consequence

On the k60 negative-lambda prefix, the exact scaled-defect cone 0<=s_k<=1 and the half-width buffer s_k<=1/2 hold on every checked row, while the one-third buffer holds on only 139/177 rows. Thus the previous one-third sufficient route is too strong for the observed prefix; any analytic tail theorem should target the exact cone or a weaker buffer such as one-half.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_k60_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```
