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
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.json --note outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-defect frontier scout: 597 scaled rows, 597 cone rows, 521 half-width rows, 179 one-third rows, 418 one-third failures, 594 scaled-increase rows, 0 issues
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
coefficient range: A_0..A_200
checked contractions: x_1..x_199
exact cone rows: 597 / 597
half-width rows: 521 / 597
one-third rows: 179 / 597
one-third failure rows: 418
scaled-defect increase rows: 594 / 594
```

Per-lambda one-third frontier:

```text
lambda=-100.0: one-third 99/199, first failure k=100, margin=-9.894755840781583542E-4, max s=4.390063006663628216E-1 at k=199
lambda=-50.0: one-third 50/199, first failure k=51, margin=-2.719207370662713847E-4, max s=5.049575716217948385E-1 at k=199
lambda=-25.0: one-third 30/199, first failure k=31, margin=-1.768768251686068314E-3, max s=5.376643171065356005E-1 at k=199
```

Global extrema:

```text
min exact-cone scaled margin: 1.346067912947116963E-2 at lambda=-100.0, k=1
min half-width margin: -3.766431710653560048E-2 at lambda=-25.0, k=199
min one-third margin: -2.043309837732022671E-1 at lambda=-25.0, k=199
max scaled defect: 5.376643171065356005E-1 at lambda=-25.0, k=199
min scaled-defect increase: 4.449655594664470276E-4 at lambda=-25.0, k=198
```

## Consequence

On the k200 negative-lambda prefix, the exact scaled-defect cone 0<=s_k<=1 holds on every checked row and the half-width buffer s_k<=1/2 fails on 76/597 checked rows, while the one-third buffer holds on only 179/597 rows. Thus the previous one-third sufficient route is too strong for the observed prefix; any analytic tail theorem should target the exact cone or a lambda/k-dependent buffer above the observed frontier.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_k200_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```
