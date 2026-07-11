# Jensen-Window PF Negative-Lambda Raw-Log Decrement Bridge

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_raw_log_decrement_bridge`.

Proof boundary: this artifact rewrites the raw-ratio decrement
corridor as a log-ratio recurrence and checks repaired k300 evidence.
It does not prove any all-`k` theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda raw-log decrement bridge: 8 rows, 0 issues, 894 log-corridor rows, 894 log-decrease rows, 2 exact counterexamples, 0 ready-to-apply rows
```

## Exact Log Form

Let:

```text
M_k = mu_{2k}
R_k = M_(k+1)*M_(k-1)/M_k^2
p_k = log(R_k)
delta_k = p_(k+1)-p_k = log(R_(k+1)/R_k)
```

The decrement corridor is exactly equivalent to:

```text
log(1-4/(2*k+1)^2) <= delta_k
delta_k <= log(1-2*(1-exp(-p_k))/(2*k+1))
```

The log-corridor width is nonnegative exactly when the raw upper wall
`R_k <= (2*k+1)/(2*k-1)` holds.

## Repaired k300 Stress

Inputs:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl
```

Repaired log diagnostics:

```text
raw-log lower wall rows: 897 / 897
raw-log upper wall rows: 897 / 897
log-decrease rows: 894 / 894
log lower-bound rows: 894 / 894
log upper-bound rows: 894 / 894
log-corridor rows: 894 / 894
log-width rows: 894 / 894
```

Extrema:

```text
min log lower margin: 4.098948243972025449E-6 at lambda=-100.0, k=298
min log upper margin: 8.786855849482974303E-7 at lambda=-25.0, k=298
min log width: 5.536125464765993858E-6 at lambda=-100.0, k=298
delta range: -5.874561412680720753E-1 at lambda=-100.0, k=1 to -5.699238767442957333E-6 at lambda=-25.0, k=298
```

## Shortcut Gate

Raw-log decrease is not enough:

```text
R_1=2, R_2=3/2: raw cone and log decrease hold, but the lower decrement wall fails
R_1=2, R_2=1: raw cone and log decrease hold, but the upper decrement wall fails
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.md
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The raw-ratio decrement corridor is exactly equivalent to a log-ratio recurrence for delta_k=log(R_(k+1)/R_k). On the repaired k300 data the log walls hold on 897/897 raw rows and the two-sided log corridor holds on 894/894 adjacent rows, but raw-log decrease alone is blocked by two exact cone counterexamples.
