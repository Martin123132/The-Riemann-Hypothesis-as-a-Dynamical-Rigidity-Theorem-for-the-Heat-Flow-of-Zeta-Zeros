# Jensen-Window PF Negative-Lambda k300 Precision-Repair Audit

Date: 2026-07-07

Status: finite precision-repair theorem-search diagnostic. This is not
a proof of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_k300_precision_repair_audit`.

Proof boundary: this artifact documents a k300 stress extension and a
local high-precision repair. It does not prove any all-`k` theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_k300_precision_repair_audit.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_k300_precision_repair_audit.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py
```

Current result:

```text
validated Jensen-window PF negative-lambda k300 precision-repair audit: 7 rows, 0 issues, 894 repaired decrement-corridor rows, 891 repaired theta-k-monotone rows, 0 ready-to-apply rows
```

## Inputs

Broad k300 run:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300_summary.json
python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-25,-50,-100 --k-min 0 --k-max 300 --n-sum 50 --cutoff 6 --dps 160 --digits 80 --abs-tol 1e-110 --constant-terms 30 --output-dir work/rh_compute/results --run-id acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300 --overwrite
```

Local repair runs:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220_summary.json
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220_summary.json
python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-100 --k-min 220 --k-max 250 --n-sum 70 --cutoff 7 --dps 220 --digits 100 --abs-tol 1e-150 --constant-terms 40 --output-dir work/rh_compute/results --run-id acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220 --overwrite
python work/rh_compute/scripts/acb_coefficient_enclosures.py --lambdas=-100 --k-min 245 --k-max 320 --n-sum 70 --cutoff 7 --dps 220 --digits 100 --abs-tol 1e-150 --constant-terms 40 --output-dir work/rh_compute/results --run-id acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220 --overwrite
```

## Broad-Run Precision Alarm

The broad dps160/cutoff6 run alone reports lambda=-100 failures:

```text
first raw decrease failure: k=252
first lower decrement failure: k=252
first upper decrement failure: k=253
first raw upper-wall failure: k=271
first theta-k monotone failure: k=233
```

These broad-run failures are treated as precision alarms, not mathematical
counterexamples.

## Repaired k300 Stress

After the local lambda=-100 repair:

```text
raw lower wall rows: 897 / 897
raw upper wall rows: 897 / 897
raw decrease rows: 894 / 894
lower decrement rows: 894 / 894
upper decrement rows: 894 / 894
decrement-corridor rows: 894 / 894
theta unit rows: 894 / 894
theta-k monotone rows: 891 / 891
theta lambda-order rows: 596 / 596
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The k300 stress extension supports the raw-ratio decrement route only after local precision repair: the broad dps160/cutoff6 run falsely reports lambda=-100 high-k failures, while the dps220/cutoff7 repair restores 897/897 raw wall rows, 894/894 decrement-corridor rows, 894/894 theta-unit rows, 891/891 theta-k monotone rows, and 596/596 theta lambda-order rows. This is finite stress evidence, not an all-k theorem.
