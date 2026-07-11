# Jensen-Window PF Negative-Lambda Raw-Moment Bridge Matrix

Date: 2026-07-06

Status: exact finite theorem-search diagnostic. This is not a proof of
the adaptive scaled-defect target, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_raw_moment_bridge_matrix`.

Proof boundary: this artifact translates the adaptive route into raw
moment-ratio inequalities and checks the finite negative-lambda prefix.
It does not prove the all-`k` raw moment-growth or corridor theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.json
```

Input enclosures:

```text
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda raw-moment bridge matrix: 8 matrix rows, 0 issues, 597 raw-cone rows, 594 corridor rows, 76 half-width failures
```

## Raw-Moment Translation

Let:

```text
M_k = mu_{2k}
A_k = M_k*k!/(2*k)!
R_k = M_(k+1)*M_(k-1)/M_k^2
x_k = ((2*k-1)/(2*k+1))*R_k
s_k = ((2*k+1)-(2*k-1)*R_k)/2
```

Then the exact cone is:

```text
0 <= s_k <= 1 iff 1 <= R_k <= (2*k+1)/(2*k-1)
```

The adaptive corridor is:

```text
((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)
corridor width = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2
```

Finite diagnostics:

```text
lambdas: -25.0, -50.0, -100.0
coefficient range: A_0..A_200
checked raw ratios: R_1..R_199
raw lower wall rows: 597 / 597
raw upper wall rows: 597 / 597
corridor occupancy rows: 594 / 594
corridor width rows: 594 / 594
half-width failure rows: 76
one-third failure rows: 418
```

Extrema:

```text
max raw ratio: 2.973078641741057661E+0 at lambda=-100.0, k=1
min raw lower margin: 2.329147017095538537E-3 at lambda=-25.0, k=199
min raw upper slack: 2.211618643155480210E-3 at lambda=-100.0, k=199
min bridge lower margin: 7.465690906027845625E-6 at lambda=-100.0, k=198
min scaled upper margin: 2.241640098067743212E-6 at lambda=-25.0, k=198
min corridor width: 1.112322974466921331E-5 at lambda=-100.0, k=198
```

Fixed-buffer translations:

```text
s_k <= 1/2 iff R_k >= 2*k/(2*k-1)
s_k <= 1/3 iff R_k >= (6*k+1)/(3*(2*k-1))
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md
outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md
outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md
```

Summary:

The adaptive exact-cone route becomes a raw-moment ratio corridor. On the k200 negative-lambda prefix, the raw cone holds on 597/597 lower rows and 597/597 upper rows, while the adaptive corridor contains R_(k+1) on 594/594 adjacent rows. The exact corridor is nonempty precisely under the upper raw wall, so the remaining all-k burden is a zeta-specific upper raw moment-growth theorem plus a corridor occupancy theorem.
