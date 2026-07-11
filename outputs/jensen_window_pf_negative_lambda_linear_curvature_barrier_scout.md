# Jensen-Window PF Negative-Lambda Linear Curvature-Barrier Scout

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_linear_curvature_barrier_scout`.

Proof boundary: this artifact isolates a linear sufficient theorem for
the lower coefficient-curvature barrier. It does not prove any all-`k`
theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda linear curvature-barrier scout: 8 rows, 0 issues, 894 linear-barrier rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows
```

## Exact Linear Sufficient Condition

Let:

```text
B_k = -log(((2*k-1)/(2*k+1))*R_k)
alpha_k = (2*k+1)/(2*k+3)
L_k(B) = log((2*k+3)/(2+(2*k+1)*exp(-B)))
```

For `B>=0`, the exact calculus lemma is:

```text
L_k(B) <= alpha_k*B
```

Therefore the simpler linear theorem

```text
B_(k+1) >= ((2*k+1)/(2*k+3))*B_k
```

implies the nonlinear lower curvature barrier.

## Repaired k300 Stress

Inputs:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl
```

Repaired linear-barrier diagnostics:

```text
B wall rows: 897 / 897
exact lower-barrier rows: 894 / 894
linear lower-barrier rows: 894 / 894
monotone-curvature rows: 894 / 894
linear B-corridor rows: 894 / 894
```

Extrema:

```text
min linear lower-barrier margin: 8.725751895568097074E-7 at lambda=-25.0, k=298
min exact lower-barrier margin: 8.786855849482974303E-7 at lambda=-25.0, k=298
min linear-to-exact slack: 4.569274603858794776E-9 at lambda=-100.0, k=298
B_(k+1)/B_k range: 9.053367907812444210E-1 at lambda=-25.0, k=1 to 9.975258371615469160E-1 at lambda=-100.0, k=298
```

## Shortcut Distinction

The analogous defect-width recurrence is already rejected:

```text
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
```

Thus the live route is specifically the log-curvature linear barrier, not
a defect-linear shortcut.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The nonlinear lower curvature barrier is implied by the simpler linear theorem B_(k+1)>=((2*k+1)/(2*k+3))*B_k. On repaired k300 data the linear lower barrier holds on 894/894 adjacent rows and the B wall holds on 897/897 rows, while the analogous defect-width recurrence is already rejected.
