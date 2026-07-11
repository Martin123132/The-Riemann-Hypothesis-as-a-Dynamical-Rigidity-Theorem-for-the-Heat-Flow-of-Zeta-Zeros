# Jensen-Window PF Negative-Lambda Scaled-Curvature Monotonicity Target

Date: 2026-07-07

Status: open theorem target. This is not a proof of the scaled-curvature
monotonicity theorem, the raw-corridor theorem, cone entry, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target`.

Proof boundary: this artifact names the replacement theorem target after
the fixed `2/3` wall was retired. It does not prove any all-`k` theorem.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: 10 rows, 0 issues, 2 live routes, 894 scaled-curvature increase rows, 0 ready-to-apply rows
```

## Target Statement

Let:

```text
R_k = M_(k+1)*M_(k-1)/M_k^2
B_k = -log(((2*k-1)/(2*k+1))*R_k)
C_k = (2*k+1)*B_k
```

The replacement target is:

```text
C_(k+1) >= C_k
```

Equivalently:

```text
B_(k+1) >= ((2*k+1)/(2*k+3))*B_k
```

This is a lower-side target for the coefficient-curvature corridor. It
must still be paired with the separate upper-side target `B_(k+1)<=B_k`.

## Repaired k300 Anchor

```text
B_k > 0 rows: 897 / 897
C_(k+1)-C_k positive rows: 894 / 894
C_(k+1)-C_k failures: 0
C_(k+1)-C_k inconclusive: 0
retired C_k<=2/3 failures: 718 / 897
min C increase margin: 5.226725385445290147E-4 at lambda=-25.0, k=298
```

## Exact Handoff

```text
C_(k+1)>=C_k
<=> B_(k+1)>=((2*k+1)/(2*k+3))*B_k
=> nonlinear lower curvature barrier, assuming B_k>=0
plus B_(k+1)<=B_k gives the two-sided curvature corridor
```

Rejected shortcuts:

```text
fixed C_k<=2/3 wall
generic raw-moment positivity or raw walls alone
monotone B curvature alone
```

Live routes:

```text
1. zeta-specific log-ratio recurrence
2. signed/tilted saddle analysis with uniform remainders
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md
outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```

Summary:

The replacement for the retired fixed 2/3 wall is the scaled-curvature monotonicity theorem target C_(k+1)>=C_k, exactly equivalent to the linear barrier B_(k+1)>=((2*k+1)/(2*k+3))*B_k. Repaired k300 data support this on 894/894 adjacent rows, while C_k<=2/3 fails on 718/897 rows. Together with B_k>=0 and the separate monotone-contraction upper side B_(k+1)<=B_k, this would feed the zeta-specific raw-corridor target.
