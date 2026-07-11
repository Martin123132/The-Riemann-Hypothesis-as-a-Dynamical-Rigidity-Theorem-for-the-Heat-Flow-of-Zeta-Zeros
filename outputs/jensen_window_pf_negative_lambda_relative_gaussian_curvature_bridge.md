# Jensen-Window PF Negative-Lambda Relative-Gaussian Curvature Bridge

Date: 2026-07-07

Status: exact finite theorem-search diagnostic. This is not a proof
of scaled-curvature monotonicity, the raw-corridor theorem, cone entry,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge`.

Proof boundary: this artifact rewrites the scaled-curvature target in
relative-Gaussian log-moment coordinates and checks repaired k300
evidence. It does not prove any all-`k` theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: 8 rows, 0 issues, 897 B-positive rows, 894 B-decrease rows, 894 C-increase rows, 598 C-lambda-order rows, 0 ready-to-apply rows
```

## Exact Relative-Gaussian Form

Let `f_k` be the log moment sequence after subtracting any Gaussian
raw-moment baseline with:

```text
G_(k+1)*G_(k-1)/G_k^2 = (2*k+1)/(2*k-1)
```

Then:

```text
B_k = h_k-log(R_k) = 2*f_k-f_(k-1)-f_(k+1)
C_k = (2*k+1)*B_k
```

The scaled-curvature target is exactly:

```text
C_(k+1)-C_k
 = (2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)
 >= 0
```

The companion monotone upper side is:

```text
B_k-B_(k+1) = -f_(k-1)+3*f_k-3*f_(k+1)+f_(k+2) >= 0
```

## Repaired k300 Stress

```text
B-positive rows: 897 / 897
B-decrease rows: 894 / 894
C-increase rows: 894 / 894
C lambda-order rows: 598 / 598
```

Sharp finite margins:

```text
min B: 1.652602130709587426E-3 at lambda=-100.0, k=299
min B decrease: 4.098948243972025449E-6 at lambda=-100.0, k=298
min C increase: 5.226725385445290147E-4 at lambda=-25.0, k=298
min C lambda gap: 3.376275223279466485E-2 at lambda-pair=-50.0 to -100.0, k=1
C range: 2.704287906143097886E-2 to 1.144219413064916367E+0
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md
outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md
outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
```

Summary:

The scaled-curvature target can be rewritten as a weighted four-point discrete-curvature inequality for the log moment sequence relative to the Gaussian baseline. Repaired k300 data validate B_k>0 on 897/897 rows, B_(k+1)<=B_k on 894/894 rows, C_(k+1)>=C_k on 894/894 rows, and C_k lambda-ordering on 598/598 rows.
