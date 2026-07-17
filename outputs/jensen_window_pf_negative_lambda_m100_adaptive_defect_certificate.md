# Jensen-Window PF Negative-Lambda -100 Adaptive-Defect Certificate

Date: 2026-07-11

Status: exact interval-analytic adaptive-defect theorem at lambda=-100.
This is not a proof of PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate`.

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.json
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda -100 adaptive-defect certificate: 8 rows, 0 issues, 2 theorem inputs, 4 defect conclusions, 0 open requirements
```

## Exact Composition

Put

```text
x_k=((2*k-1)/(2*k+1))*R_k,
d_k=1-x_k,
s_k=((2*k+1)/2)*d_k.
```

Full ratio-cone entry at lambda=-100 proves

```text
(2*k-1)/(2*k+1)<=x_k<=1,
x_(k+1)>=x_k.
```

Therefore, for every `k>=1`,

```text
0<=d_k<=2/(2*k+1),
d_(k+1)<=d_k,
0<=s_k<=1.
```

The upper raw-corridor wall gives the additional exact identity

```text
s_(k+1)-s_k
  =(2+(2*k-1)*R_k-(2*k+1)*R_(k+1))/2
  >=0.
```

Thus the lambda=-100 entry route satisfies the defect-tail theorem from
`k=1` and the adaptive exact-cone plus monotone-bridge target from `k=1`.
The older targets asked for the stronger simultaneous statement at
`lambda=-25,-50,-100`; that simultaneous theorem is not proved, but it is
not needed once one full entry parameter has been established at `lambda=-100`.

Higher-degree minor cones and the all-order Jensen/PF bridge remain open.

```text
outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md
outputs/jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.md
outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```
