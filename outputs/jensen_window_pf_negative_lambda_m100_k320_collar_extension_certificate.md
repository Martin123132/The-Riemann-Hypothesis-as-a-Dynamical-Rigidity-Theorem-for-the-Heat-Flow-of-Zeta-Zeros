# Jensen-Window PF Negative-Lambda -100 k320 Collar Extension Certificate

Date: 2026-07-10

Status: finite Arb cone-collar extension certificate. This is not a proof
of an eventual-k theorem, cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate`.

Machine-readable result and source:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.json
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl
work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220_summary.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: 6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows
```

## Certified Extension

The repaired dps220 source already covered `A_245..A_320`. Direct Arb
ratio arithmetic certifies

```text
(2*k-1)/(2*k+1) < x_k < 1,       k=246..319
x_(k+1)>x_k,                      k=246..318.
```

The second line promotes 19 previously unused adjacent rows
`k=300..318`. The weakest new logarithmic gap is

```text
L_318 > 3.70927257372186804057703586992264447737838968980644238576986E-6
```

where `L_k=log(x_(k+1)/x_k)`.

Selected extension rows:

| k | Arb log gap |
|---:|---:|
| `300` | `[4.05724319967749486814482971269169621E-6 +/- 2.03E-42]` |
| `301` | `[4.03663227944111442068783914814591570E-6 +/- 2.35E-42]` |
| `310` | `[3.858089346138218497504010429984818E-6 +/- 9.37E-40]` |
| `317` | `[3.72739268371590992292160965714084E-6 +/- 5.12E-39]` |
| `318` | `[3.70927257372186804057703586992265E-6 +/- 5.53E-39]` |

## Handoff

The finite lambda=-100 collar now reaches `k=318`. The analytic tail
must start at `k=319`; the all-k first-summand dominance certificate
already controls the full-kernel perturbation there, but the dominant
`n=1` adjacent-wall lower bound remains open.

```text
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The existing repaired dps220 lambda=-100 source contains more rigorous information than the k300 audit exposed. Arb certifies 76 positive coefficients, both cone walls for x_246..x_319, and the adjacent wall for k=246..318. The 19 newly promoted rows k=300..318 have minimum adjacent log gap above 3.709e-6 at k=318, extending the finite handoff while leaving the analytic k>=319 tail open.
