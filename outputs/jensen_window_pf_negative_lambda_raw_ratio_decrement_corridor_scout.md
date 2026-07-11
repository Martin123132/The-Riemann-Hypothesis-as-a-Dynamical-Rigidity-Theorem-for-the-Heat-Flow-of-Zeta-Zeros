# Jensen-Window PF Negative-Lambda Raw-Ratio Decrement-Corridor Scout

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of the raw-corridor theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout`.

Proof boundary: this artifact rewrites the zeta-specific raw corridor
as a decrement recurrence and validates finite shape evidence. It does
not prove an all-`k` recurrence.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: 9 rows, 0 issues, 594 decrement-corridor rows, 591 theta-k-monotone rows, 2 exact counterexamples, 0 ready-to-apply rows
```

## Exact Decrement Form

Let:

```text
M_k = mu_{2k}
R_k = M_(k+1)*M_(k-1)/M_k^2
D_k = R_k - R_(k+1)
```

The adaptive raw corridor is exactly:

```text
2*(R_k-1)/(2*k+1) <= D_k <= 4*R_k/(2*k+1)^2
```

Its width is:

```text
4*R_k/(2*k+1)^2 - 2*(R_k-1)/(2*k+1)
= 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2
```

Finite diagnostics:

```text
lambdas: -25.0, -50.0, -100.0
coefficient range: A_0..A_200
checked raw ratios: R_1..R_199
raw decrease rows: 594 / 594
lower decrement rows: 594 / 594
upper decrement rows: 594 / 594
decrement-corridor rows: 594 / 594
theta unit rows: 594 / 594
theta-k monotone rows: 591 / 591
theta lambda-order rows: 396 / 396
```

Extrema:

```text
min lower decrement margin: 2.241640098067743212E-6 at lambda=-25.0, k=198
min upper decrement margin: 7.465690906027845625E-6 at lambda=-100.0, k=198
min decrement corridor width: 1.112322974466921331E-5 at lambda=-100.0, k=198
theta range: 1.644126617921189501E-1 at lambda=-25.0, k=198 to 9.087309076587821269E-1 at lambda=-100.0, k=1
min theta k-drop: 3.722452214840226847E-4 at lambda=-25.0, k=197
min theta lambda gap: 5.048218901591782104E-2 at pair=-25.0->-50.0, k=198
```

Exact shortcut blockers:

```text
R_1=2, R_2=3/2: raw cone and decrease hold, but the lower decrement wall fails
R_1=2, R_2=1: raw cone and decrease hold, but the upper decrement wall fails
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md
outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The zeta-specific raw-corridor target can be attacked as a decrement recurrence: prove 2*(R_k-1)/(2*k+1) <= R_k-R_(k+1) <= 4*R_k/(2*k+1)^2. The k200 prefix satisfies this on 594 adjacent rows, and the normalized theta coordinate is k-monotone on 591 rows and lambda-ordered on 396 rows. Exact raw-cone monotone counterexamples show that the recurrence must use zeta-specific structure.
