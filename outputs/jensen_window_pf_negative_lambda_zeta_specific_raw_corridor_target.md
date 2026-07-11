# Jensen-Window PF Negative-Lambda Zeta-Specific Raw-Corridor Target

Date: 2026-07-06

Status: open theorem target. This is not a proof of cone entry,
Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target`.

Proof boundary: this artifact names the missing zeta-specific raw-wall
and adaptive-corridor theorem. It does not prove that theorem.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: 9 rows, 0 issues, 2 live routes, 2 rejected shortcuts, 0 ready-to-apply rows
```

## 2026-07-10 Upper-Wall Handoff

The kernel Mellin theorem certificate
`outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md` proves the
raw upper wall `R_k<=(2*k+1)/(2*k-1)` for every real heat parameter and every
`k>=1`. The original target statement below is retained for audit history;
the current missing theorem is only occupancy of the adjacent raw-ratio
corridor, equivalently the monotone-contraction and scaled-envelope sides.

## 2026-07-10 Parameter-Scope Correction

The actual-kernel certificate
`outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md`
proves that corridor occupancy fails at `lambda=-1156`, `k=120`. Therefore
the raw corridor cannot be demanded uniformly over every negative heat
parameter. The viable target is occupancy at one selected moderate negative
entry parameter, with a certified finite collar and a separate eventual-`k`
zeta saddle/tail theorem.

For the candidate `lambda=-100` handoff, the exact summand-shift lemma
`outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md`
reduces the eventual-`k` theorem to a dominant `n=1` saddle inequality and a
shifted `v>=3/2` far-tail budget. Its compact `n=2..20` contribution is
already certified below `2.122e-29` for every `k>=300`.

The subsequent first-summand dominance certificate closes the complete
kernel-tail source, proving relative moment error at most `2/k^6` and
adjacent log-wall error at most `16/(k-1)^6`. Together with the finite collar
through `k=318`, the raw-corridor tail now reduces to the single open target

```text
L_k^(1)>=1/(4*k^2), k>=319.
```

The exact cumulant bridge sharpens that target to the sufficient estimate

```text
kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318.
```

Its Gamma and rational-transfer algebra is proved; the uniform special-kernel
cumulant estimate remains open. The leading-saddle certificate proves the
leading cap `13/20`, cubic-correction cap `1/100`, and fifth-order cap
`1/1000`. The paired theorem proves the seventh-order normalized remainder
floor `-79/1000` on `0.9264<=u<=5`, reducing the tail to the ray `u>=5`.
The analytic ray certificate proves `H_t>=-3/250` there, closing the
lambda=-100 adjacent raw wall; the broader zeta-specific corridor remains open.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
```

## Target Statement

```text
M_k = mu_{2k}
R_k = M_(k+1)*M_(k-1)/M_k^2
1 <= R_k <= (2*k+1)/(2*k-1)
((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)
```

Finite support:

```text
k200 raw-cone rows: 597 / 597
k200 corridor rows: 594 / 594
```

Rejected shortcuts:

```text
generic Stieltjes/raw-log-convexity proof
positive Gaussian scale-mixture proof of the upper wall
```

Live routes:

```text
1. signed Gaussian perturbation plus two-scale uniform remainders
2. direct zeta-specific ratio recurrence compatible with increasing scaled defect
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md
outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md
outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
```

Summary:

The zeta-specific raw-corridor theorem target is now explicit: prove the upper raw wall and adaptive corridor for actual negative-lambda zeta moments. The k200 prefix supports the target with 597 raw-cone rows and 594 corridor rows, but generic Stieltjes positivity and positive Gaussian mixtures are blocked; the live routes are signed Gaussian perturbation with two-scale remainders or a direct zeta-specific ratio recurrence.
