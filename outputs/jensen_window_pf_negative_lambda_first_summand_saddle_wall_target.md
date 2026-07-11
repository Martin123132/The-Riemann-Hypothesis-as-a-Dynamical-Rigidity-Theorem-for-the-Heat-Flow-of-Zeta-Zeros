# Jensen-Window PF Negative-Lambda First-Summand Saddle-Wall Target

Date: 2026-07-10

Status: exact saddle geometry with analytic wall closure; not a proof of the remaining all-order bridge.
This proves the lambda=-100 adjacent wall, not the remaining all-order bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_first_summand_saddle_wall_target`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda first-summand saddle-wall closure: 9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, 0 open requirements, 2 ready-to-apply rows
```

## Exact Geometry

For the first Newman summand at `lambda=-100`, let

```text
S_k(u)=2*k*log(u)-100*u^2+log(phi_1(u)).
```

With `q=pi*exp(4u)`,

```text
S_k'(u)=2*k/u-200*u+5+8*q/(2*q-3)-4*q,
S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0.
```

Thus there is one saddle `s_k`. Exact endpoint propagation proves

```text
(log(k)+22/25)/8 < s_k < log(k)/4, k>=300.
```

## High-Precision Scout

The scout evaluates the first-summand moments with 60-digit mpmath
saddle-centered quadrature. The finite upper truncation is not certified,
so these rows are theorem-search evidence only.

| k | saddle | L_k^(1) | k^2 L_k^(1) | L_k^(1)-1/(4k^2) |
|---:|---:|---:|---:|---:|
| `300` | `0.91270928796298985056807425462483823` | `4.0572431996774948681448297126916962e-6` | `0.36515188797097453813303467414225266` | `1.2794654218997170903670519349139184e-6` |
| `400` | `0.97888094701689181843228254728997258` | `2.5810854757523237610846318406240925e-6` | `0.4129736761203718017735410944998548` | `1.0185854757523237610846318406240925e-6` |
| `500` | `1.030077594236378203800516720387092` | `1.7900980615579707918656351323123603e-6` | `0.44752451538949269796640878307809008` | `7.9009806155797079186563513231236033e-7` |
| `700` | `1.1069452994078673013038750783415861` | `1.0098324178528385723969950897527977e-6` | `0.49481788474789090047452759397887087` | `4.9962833622018551117250529383443035e-7` |
| `1000` | `1.1879375411452161946579555446887138` | `5.3852461092523267042219952610500539e-7` | `0.53852461092523267042219952610500539` | `2.8852461092523267042219952610500539e-7` |
| `2000` | `1.3440078260003367251495307467509845` | `1.5160662670416244357565390578729955e-7` | `0.60642650681664977430261562314919819` | `8.9106626704162443575653905787299549e-8` |
| `5000` | `1.5485738239457640734903380663464682` | `2.6797224818459480525570040725290535e-8` | `0.66993062046148701313925101813226337` | `1.6797224818459480525570040725290535e-8` |
| `10000` | `1.7028766831343941982947111923584631` | `7.0465686623400990046495219175233875e-9` | `0.70465686623400990046495219175233875` | `4.5465686623400990046495219175233875e-9` |
| `20000` | `1.8573688080533650862594433804397783` | `1.8297016537092515733142790645243493e-9` | `0.7318806614837006293257116258097397` | `1.2047016537092515733142790645243493e-9` |

All nine rows lie above `1/(4*k^2)`, and the sampled scaled profile
is strictly increasing. At `k=300`, the first-summand value agrees
with the repaired full-kernel Arb log gap to the displayed precision:

```text
full Arb: [4.05724319967749486814482971269169621E-6 +/- 2.03E-42]
n=1 mp:   4.0572431996774948681448297126916962e-6
```

This agreement is a cross-check, not an interval promotion of the mpmath
quadrature.

## Analytic Closure

The exact target is

```text
L_k^(1)>=1/(4*k^2), k>=319.
```

The paired interval theorem proves the seventh-order normalized remainder
floor `-79/1000` for every real mode `0.9264<=u<=5`, reaching
`t=V'(x(5))>1.5241916613e10`. The analytic ray certificate proves
the same paired remainder floor on `u>=5`, closing the half-line.
The later exact cumulant bridge sharpens this to the sufficient estimate
`kappa_3,t(2*log(U))>=-37/(50*t^2)` for every real `t>=318`.
The leading-saddle certificate proves caps 13/20, 1/100, and 1/1000
through fifth order; the compact and ray certificates close the remainder.

The first-summand dominance certificate gives

```text
|L_k-L_k^(1)|<=16/(k-1)^6.
```

Since `1/(4*k^2)>16/(k-1)^6` for every `k>=319`, the proved
dominant target closes the full-kernel tail and splices to the certified
finite collar through `k=318`.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The dominant n=1 integrand has exact strict-concavity geometry and an all-k saddle bracket. Sixty-digit saddle-centered quadrature at nine k values from 300 to 20000 finds L_k^(1)>1/(4*k^2), with k^2*L_k^(1) increasing from about 0.365 to 0.732 and formally tending to one. This remains finite floating evidence. The paired compact and ray theorems prove the analytic tail, and the certified full-kernel perturbation plus finite collar close the adjacent wall.
