# Jensen-Window PF Negative-Lambda Finite-Collar Contract

Date: 2026-07-06

Status: finite-collar accounting contract. This is not a proof of zeta cone
entry, monotone contractions, Jensen-window PF-infinity, RH, or the
Newman-direction goal.

Artifact kind: `jensen_window_pf_negative_lambda_finite_collar_contract`.

Proof boundary: this artifact extracts the finite active depth currently
supported by the negative-lambda prefix and the exact ratio-cone collar
requirements. It does not prove an infinite cone-entry theorem, does not prove
an all-`k` tail theorem, does not prove a finite-collar flow theorem beyond the
stated depth, does not prove `jwpf_06`, and does not settle `Lambda <= 0`.

Machine-readable contract:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k80_contract.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_finite_collar_contract.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k80_contract.json --note outputs/jensen_window_pf_negative_lambda_finite_collar_k80_contract.md
```

Current result:

```text
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=77, 231 active lower rows, 231 active upper rows, 231 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
```

## Collar Rule

The exact ratio-cone invariance lemma states:

```text
checking x_1..x_K requires collar variables x_{K+1} for lower/upper walls
and x_{K+2} for monotone wall tests
```

For the current negative-lambda prefix:

```text
lambdas: -100.0, -50.0, -25.0
available coefficient range: A_0..A_80
available ratio range: x_1..x_79
lower/upper range: x_1..x_79
monotone-gap range: 1 <= k <= 78

certified active depth: K=77
first collar: x_78
second collar: x_79
```

## Certified Rows

```text
active lower rows: 231 / 231
active upper rows: 231 / 231
active monotone rows: 231 / 231
collar lower/upper rows: 6 / 6
collar monotone rows: 3 / 3

minimum active lower margin:
  7.190992973645317614E-3 at lambda=-25.0, k=77

minimum active upper margin:
  3.791406490745955089E-3 at lambda=-100.0, k=77

minimum active monotone gap:
  2.311788760235680901E-5 at lambda=-100.0, k=77

minimum collar lower/upper margin:
  3.768288603143598280E-3 at lambda=-100.0, k=78

minimum collar monotone gap:
  2.281394370544295027E-5 at lambda=-100.0, k=78
```

## Meaning

The negative-lambda prefix is now precisely collar-accounted: it can seed the
active finite depth `K=77` with collars
`x_78` and `x_79`.
It does not seed `K=78` because the
monotone-wall collar would require `x_80`,
hence additional coefficient enclosures through
`A_81`, or an analytic tail theorem.

Related artifacts:

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```
