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
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k150.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.json --note outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues
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
coefficient range: A_0..A_150
checked contractions: x_1..x_149
exact cone rows: 447 / 447
half-width rows: 430 / 447
one-third rows: 179 / 447
one-third failure rows: 268
scaled-defect increase rows: 444 / 444
```

Per-lambda one-third frontier:

```text
lambda=-100.0: one-third 99/149, first failure k=100, margin=-9.894755840781583542E-4, max s=3.960755868223805985E-1 at k=149
lambda=-50.0: one-third 50/149, first failure k=51, margin=-2.719207370662713847E-4, max s=4.722996954985768865E-1 at k=149
lambda=-25.0: one-third 30/149, first failure k=31, margin=-1.768768251686068314E-3, max s=5.110735937007886999E-1 at k=149
```

Global extrema:

```text
min exact-cone scaled margin: 1.346067912947116963E-2 at lambda=-100.0, k=1
min half-width margin: -1.107359370078869990E-2 at lambda=-25.0, k=149
min one-third margin: -1.777402603674553666E-1 at lambda=-25.0, k=149
max scaled defect: 5.110735937007886999E-1 at lambda=-25.0, k=149
min scaled-defect increase: 6.439210603944736609E-4 at lambda=-25.0, k=148
```

## Consequence

On the k150 negative-lambda prefix, the exact scaled-defect cone 0<=s_k<=1 holds on every checked row and the half-width buffer s_k<=1/2 fails on 17/447 checked rows, while the one-third buffer holds on only 179/447 rows. Thus the previous one-third sufficient route is too strong for the observed prefix; any analytic tail theorem should target the exact cone or a lambda/k-dependent buffer above the observed frontier.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```
