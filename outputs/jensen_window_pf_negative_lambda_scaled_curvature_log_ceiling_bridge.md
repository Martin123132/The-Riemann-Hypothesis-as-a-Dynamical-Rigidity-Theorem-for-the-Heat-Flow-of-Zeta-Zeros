# Jensen-Window PF Negative-Lambda Scaled-Curvature Log-Ceiling Bridge

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of scaled-curvature monotonicity, the raw-corridor theorem, cone entry,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge`.

Proof boundary: this artifact rewrites scaled-curvature monotonicity
as an affine log-ratio ceiling and checks repaired k300 evidence. It
does not prove any all-`k` theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: 8 rows, 0 issues, 894 scaled-ceiling rows, 894 scaled-log-corridor rows, 894 ceiling-dominance rows, 0 ready-to-apply rows
```

## Exact Affine Ceiling

Let:

```text
p_k = log(R_k)
delta_k = p_(k+1)-p_k
h_k = log((2*k+1)/(2*k-1))
B_k = h_k-p_k
C_k = (2*k+1)*B_k
```

Then:

```text
C_(k+1)>=C_k
<=> delta_k <= h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)
```

On the raw upper wall `B_k>=0`, this affine ceiling is stronger than
the nonlinear raw-log upper wall:

```text
raw_log_upper(p_k)-scaled_ceiling(p_k)=alpha_k*B_k-L_k(B_k)>=0
```

## Repaired k300 Stress

```text
p-wall rows: 897 lower and 897 upper / 897
scaled-ceiling rows: 894 / 894
log-lower rows: 894 / 894
scaled-log-corridor rows: 894 / 894
ceiling-dominance rows: 894 / 894
```

Sharpness:

```text
min scaled ceiling margin: 8.725751895568097074E-7 at lambda=-25.0, k=298
min raw upper margin: 8.786855849482974303E-7 at lambda=-25.0, k=298
min raw-upper minus scaled-ceiling slack: 4.569274603858794776E-9 at lambda=-100.0, k=298
min scaled log width: 5.531556190162135064E-6 at lambda=-100.0, k=298
```

Live recurrence target:

```text
log(1-4/(2*k+1)^2) <= delta_k
delta_k <= h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md
outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md
outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```

Summary:

Scaled-curvature monotonicity is exactly equivalent to the affine log-ratio ceiling delta_k<=h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3). On repaired k300 data the p-wall holds on 897/897 rows, the affine ceiling holds on 894/894 adjacent rows, and the ceiling is stronger than the nonlinear raw-log upper wall on 894/894 rows.
