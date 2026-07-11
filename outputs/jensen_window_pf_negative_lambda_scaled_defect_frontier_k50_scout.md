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
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k50.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.json --note outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.md
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 147 scaled rows, 147 cone rows, 147 half-width rows, 128 one-third rows, 19 one-third failures, 144 scaled-increase rows, 0 issues
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
coefficient range: A_0..A_50
checked contractions: x_1..x_49
exact cone rows: 147 / 147
half-width rows: 147 / 147
one-third rows: 128 / 147
one-third failure rows: 19
scaled-defect increase rows: 144 / 144
```

Per-lambda one-third frontier:

```text
lambda=-100.0: one-third 49/49, first failure none, max s=2.276063817976405457E-1 at k=49
lambda=-50.0: one-third 49/49, first failure none, max s=3.281709510724553816E-1 at k=49
lambda=-25.0: one-third 30/49, first failure k=31, margin=-1.768768251686068314E-3, max s=3.905562462337887786E-1 at k=49
```

Global extrema:

```text
min exact-cone scaled margin: 1.346067912947116963E-2 at lambda=-100.0, k=1
min half-width margin: 1.094437537662112214E-1 at lambda=-25.0, k=49
min one-third margin: -5.722291290045544523E-2 at lambda=-25.0, k=49
max scaled defect: 3.905562462337887786E-1 at lambda=-25.0, k=49
min scaled-defect increase: 2.455985080948463422E-3 at lambda=-25.0, k=48
```

## Consequence

On the k50 negative-lambda prefix, the exact scaled-defect cone 0<=s_k<=1 and the half-width buffer s_k<=1/2 hold on every checked row, while the one-third buffer holds on only 128/147 rows. Thus the previous one-third sufficient route is too strong for the observed prefix; any analytic tail theorem should target the exact cone or a weaker buffer such as one-half.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_k50_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```
