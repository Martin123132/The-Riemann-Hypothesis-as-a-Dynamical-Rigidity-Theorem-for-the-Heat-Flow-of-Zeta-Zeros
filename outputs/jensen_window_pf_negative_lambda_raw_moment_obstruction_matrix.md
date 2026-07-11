# Jensen-Window PF Negative-Lambda Raw-Moment Obstruction Matrix

Date: 2026-07-06

Status: exact countermodel gate. This is not a proof of the adaptive
scaled-defect target, cone entry, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix`.

Proof boundary: this artifact blocks generic Stieltjes/raw-log-convexity
proofs of the raw upper wall or adaptive corridor. It is not evidence
against the actual zeta finite prefix.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda raw-moment obstruction matrix: 7 matrix rows, 0 issues, 3 exact counterexamples, 0 ready-to-apply rows
```

## Exact Witnesses

All witnesses are positive two-atom Stieltjes moment sequences:

```text
M_k = 1 + w*a^k
R_k = M_(k+1)*M_(k-1)/M_k^2
```

Upper raw wall failure:

```text
measure: delta_1 + (1/16)*delta_16
R_1 = 289/64
upper wall at k=1 = 3
x_1 = 289/192
s_1 = -97/128
```

Scaled-upper corridor failure while the raw upper wall holds:

```text
measure: delta_1 + (1)*delta_2
R_1 = 10/9
R_2 = 27/25
corridor upper = 28/27
upper margin = -29/675
s_2-s_1 = -29/450
```

Monotone-bridge corridor failure while the raw upper wall holds:

```text
measure: delta_1 + (1/3)*delta_9
R_1 = 7/3
R_2 = 61/49
corridor lower = 35/27
lower margin = -68/1323
x_2-x_1 = -68/2205
```

Exact identity:

```text
corridor width = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2
```

Zeta contrast:

```text
outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md
validated 597/597 raw-cone rows and 594/594 corridor rows on the k200 prefix
```

Summary:

Generic positive Stieltjes/raw-moment structure does not prove the adaptive raw wall or corridor. Exact two-atom measures can violate the upper raw wall, the scaled-upper corridor side, or the monotone-bridge lower corridor side, even while raw log-convexity holds. The actual zeta k200 prefix remains finite-compatible, so the missing proof must be zeta-specific and all-k.
