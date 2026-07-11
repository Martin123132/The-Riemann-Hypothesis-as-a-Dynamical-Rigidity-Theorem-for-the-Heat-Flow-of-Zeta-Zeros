# Jensen-Window PF T=1156 Monotone-Wall Counterexample Certificate

Date: 2026-07-10

Status: interval zeta-kernel counterexample certificate. This is not a
proof or disproof of RH, and it does not prove `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate`.

Machine-readable result and source enclosures:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.json
work/rh_compute/results/acb_enclosures_lambda_m1156_k119_k122_dps250.jsonl
work/rh_compute/results/acb_enclosures_lambda_m1156_k119_k122_dps250_summary.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py
```

Current result:

```text
validated Jensen-window PF T=1156 monotone-wall counterexample certificate: 7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation
```

## Certified Violation

For the actual Newman kernel at `lambda=-1156`, define

```text
x_120=A_121*A_119/A_120^2
x_121=A_122*A_120/A_121^2.
```

ACB integration of `A_119..A_122`, composed with analytic n-series and
u-tail bounds, certifies

```text
x_120 = [0.99988134408769009708577622685369921860954640 +/- 4.16E-45]
x_121 = [0.9998813272439776665767885382673292264809601 +/- 4.01E-44]
x_121-x_120 = [-1.68437124305089876885863699921285863E-8 +/- 3.82E-44]
log(x_121/x_120) = [-1.68457114156376793080407845822737199E-8 +/- 8.14E-44]
```

In particular, `x_121-x_120 < -1.684371243050898768858636999212858627089898559583466809242237561908850E-8`. Both contractions
are strictly between zero and one, so this is specifically a failure of
the adjacent-k wall `x_(k+1)>=x_k`, not of the Mellin upper wall.

## Route Consequence

The fixed-k=22 certificate for every real `T>=1156` remains valid. This
counterexample shows that it cannot be extended into an all-k cone theorem
at `T=1156`: the actual sequence leaves the monotone wall at `k=120`.

The surviving cone-entry programme is now:

```text
1. choose a moderate fixed negative lambda with a deep certified finite collar;
2. prove eventual-k adjacent monotonicity there by a zeta saddle/tail theorem;
3. splice the finite collar to that tail and then invoke cone invariance.
```

The existing `lambda=-25,-50,-100` k300 enclosures support this revised
search, but they remain finite evidence and do not prove the tail theorem.

```text
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

Rigorous ACB coefficient enclosures for the actual Newman kernel at lambda=-1156 give x_121-x_120<0, with an upper enclosure below -1.68e-8. Thus the adjacent-k monotone wall fails at k=120 even though both contractions remain below one. The fixed-k=22 T>=1156 certificate is valid, but T=1156 cannot be an all-k cone-entry time.
