# Jensen-Window PF Negative-Lambda Coefficient-Curvature Corridor Bridge

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge`.

Proof boundary: this artifact rewrites the raw-ratio decrement
corridor as a coefficient-curvature corridor and checks repaired k300
evidence. It does not prove any all-`k` theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: 9 rows, 0 issues, 894 curvature-corridor rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows
```

## Exact Curvature Form

Let:

```text
R_k = M_(k+1)*M_(k-1)/M_k^2
x_k = ((2*k-1)/(2*k+1))*R_k
B_k = -log(x_k)
```

Then:

```text
R_k = ((2*k+1)/(2*k-1))*exp(-B_k)
p_k = log(R_k) = log((2*k+1)/(2*k-1))-B_k
```

The raw decrement corridor is exactly equivalent to:

```text
log((2*k+3)/(2+(2*k+1)*exp(-B_k))) <= B_(k+1)
B_(k+1) <= B_k
```

The upper side `B_(k+1)<=B_k` is the monotone-contraction side
already isolated in log-curvature language; the new missing side is
the lower curvature barrier.

## Repaired k300 Stress

Inputs:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl
```

Repaired curvature diagnostics:

```text
B positive rows: 897 / 897
B upper-wall rows: 897 / 897
monotone-curvature rows: 894 / 894
lower-barrier rows: 894 / 894
curvature-corridor rows: 894 / 894
curvature-width rows: 894 / 894
```

Extrema:

```text
min B lower-barrier margin: 8.786855849482974303E-7 at lambda=-25.0, k=298
min B monotone margin: 4.098948243972025449E-6 at lambda=-100.0, k=298
min B corridor width: 5.536125464765993858E-6 at lambda=-100.0, k=298
B range: 1.652602130709587426E-3 at lambda=-100.0, k=299 to 3.524954269540071057E-2 at lambda=-25.0, k=1
```

## Shortcut Gates

Simple curvature shortcuts are blocked:

```text
R_1=2, R_2=3/2: raw walls and B_(k+1)<=B_k hold, but the lower curvature barrier fails
R_1=2, R_2=1: raw walls hold, but B_(k+1)<=B_k fails
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The raw-corridor theorem is exactly equivalent to the coefficient-curvature corridor log((2*k+3)/(2+(2*k+1)*exp(-B_k))) <= B_(k+1) <= B_k. On repaired k300 data the B wall holds on 897/897 rows and the curvature corridor holds on 894/894 adjacent rows; monotone curvature alone is blocked by an exact cone counterexample.
