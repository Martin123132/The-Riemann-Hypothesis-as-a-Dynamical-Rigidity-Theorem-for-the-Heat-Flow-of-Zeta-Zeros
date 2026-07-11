# Jensen-Window PF Negative-Lambda Bounded Log-Curvature k300 Obstruction

Date: 2026-07-07

Status: finite obstruction gate. This is not a proof of the replacement
raw-corridor theorem, cone entry, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction`.

Proof boundary: this artifact retires the fixed `2/3` scaled-curvature
wall as a live target. It does not prove the newer linear curvature
barrier or any all-`k` theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py
```

Current result:

```text
validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows
```

## Exact Rewrite

Let:

```text
B_k = -log(((2*k-1)/(2*k+1))*R_k)
C_k = (2*k+1)*B_k
```

Then the former bounded log-curvature wall is exactly:

```text
B_k <= 2/(3*(2*k+1))  iff  C_k <= 2/3
```

## Repaired k300 Obstruction

Inputs:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl
```

Classified rows:

```text
scaled-curvature rows: 897
C_k <= 2/3 rows: 179 / 897
C_k > 2/3 rows: 718 / 897
inconclusive rows: 0
B_k > 0 rows: 897 / 897
```

First checked failure:

```text
lambda=-25.0, k=31
C_k = 6.737945594953250214E-1
2/3 - C_k = -7.127892828658354747E-3
B_k = 1.069515173802103209E-2
```

Worst checked slack:

```text
max C_k = 1.144219413064916367E+0 at lambda=-25.0, k=299
min 2/3-C_k = -4.775527463982497008E-1 at lambda=-25.0, k=299
```

Replacement-route contrast:

```text
C_(k+1)-C_k positive rows: 894 / 894
B_k-B_(k+1) positive rows: 894 / 894
min C increase margin: 5.226725385445290147E-4 at lambda=-25.0, k=298
min B decrease margin: 4.098948243972025449E-6 at lambda=-100.0, k=298
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The fixed bounded log-curvature wall B_k<=2/(3*(2*k+1)), equivalently C_k=(2*k+1)*B_k<=2/3, is finite-rejected by the repaired k300 data: only 179/897 checked rows satisfy it, 718/897 fail it, and zero rows are inconclusive. The live replacement direction is the zeta-specific linear/scaled curvature corridor, not a fixed 2/3 wall.
