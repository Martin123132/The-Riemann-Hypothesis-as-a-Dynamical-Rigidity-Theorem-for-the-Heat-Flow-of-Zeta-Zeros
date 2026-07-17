# Jensen-Window PF Negative-Lambda -100 Raw-Corridor Certificate

Date: 2026-07-11

Status: exact interval-analytic raw-corridor theorem at lambda=-100.
This is not a proof of PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_m100_raw_corridor_certificate`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda -100 raw-corridor certificate: 6 rows, 0 issues, 2 theorem inputs, 1 raw-corridor theorem, 0 open requirements
```

## Exact Composition

Put

```text
R_k=M_(k+1)*M_(k-1)/M_k^2,
B_k=-log(((2*k-1)/(2*k+1))*R_k).
```

Full ratio-cone entry at lambda=-100 proves

```text
B_k>=0, B_(k+1)<=B_k.
```

The scaled-curvature continuous bridge proves

```text
B_(k+1)>=((2*k+1)/(2*k+3))*B_k.
```

For `B>=0`, the exact calculus inequality

```text
log((2*k+3)/(2+(2*k+1)*exp(-B)))
<=((2*k+1)/(2*k+3))*B
```

follows because equality holds at zero and the left derivative decreases from
`(2*k+1)/(2*k+3)`. Therefore

```text
log((2*k+3)/(2+(2*k+1)*exp(-B_k)))<=B_(k+1)<=B_k.
```

Returning to raw ratios gives, for every `k>=1`,

```text
((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k<=R_(k+1)
R_(k+1)<=(2+(2*k-1)*R_k)/(2*k+1).
```

This closes the zeta-specific raw-corridor target at lambda=-100. Higher-degree
minor cones and the all-order Jensen/PF bridge remain open.

```text
outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
