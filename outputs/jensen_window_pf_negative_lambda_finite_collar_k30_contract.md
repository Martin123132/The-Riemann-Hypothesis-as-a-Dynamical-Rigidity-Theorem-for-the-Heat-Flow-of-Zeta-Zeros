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
work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k30_contract.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_finite_collar_contract.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py --target work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k30_contract.json --note outputs/jensen_window_pf_negative_lambda_finite_collar_k30_contract.md
```

Current result:

```text
validated Jensen-window PF negative-lambda finite-collar contract: active depth K=27, 81 active lower rows, 81 active upper rows, 81 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues
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
available coefficient range: A_0..A_30
available ratio range: x_1..x_29
lower/upper range: x_1..x_29
monotone-gap range: 1 <= k <= 28

certified active depth: K=27
first collar: x_28
second collar: x_29
```

## Certified Rows

```text
active lower rows: 81 / 81
active upper rows: 81 / 81
active monotone rows: 81 / 81
collar lower/upper rows: 6 / 6
collar monotone rows: 3 / 3

minimum active lower margin:
  2.479278584529633053E-2 at lambda=-25.0, k=27

minimum active upper margin:
  5.609413704642829646E-3 at lambda=-100.0, k=27

minimum active monotone gap:
  5.917631559591166666E-5 at lambda=-100.0, k=27

minimum collar lower/upper margin:
  5.550237389046917980E-3 at lambda=-100.0, k=28

minimum collar monotone gap:
  5.743772179246203517E-5 at lambda=-100.0, k=28
```

## Meaning

The negative-lambda prefix is now precisely collar-accounted: it can seed the
active finite depth `K=27` with collars
`x_28` and `x_29`.
It does not seed `K=28` because the
monotone-wall collar would require `x_30`,
hence additional coefficient enclosures through
`A_31`, or an analytic tail theorem.

Related artifacts:

```text
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md
outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```
