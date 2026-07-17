# Jensen-Window PF Order-Four Localized Curvature Compact Certificate

Date: 2026-07-12

Status: rigorous compact interval certificate with an open analytic ray.
This is not a proof of complete order-four entry, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_localized_curvature_compact_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json
work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_h_tiles.jsonl
python work/rh_compute/scripts/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.py
```

## Certified Range

The cached Arb quadrature tiles and localized same-point inequality prove

```text
K_1(t)<=7/(2*t^2), 319<=t<=V'(2).
```

The lower mode endpoint has `t` ball `[318.99505050020666714217851899922620279230505439154 +/- 3.20E-48]` and
the handoff endpoint has `t` ball `[37850.322210211384816295110892561323868727705625771 +/- 3.22E-46]`.
The lower endpoint lies below 319, so the cover contains the full required start.

## Interval Construction

The deterministic cache contains `107452` adjacent mode tiles of width
`1/100000`. Each tile uses `200` paired
composite-Simpson panels on `|y|<=6`, the proved `V^(8)/a^4<=1/50000`
envelope, and outward-rounded Arb arithmetic through standardized moment order eight.
The expensive `H^(2),...,H^(8)` enclosures are cached once and reused.

Adaptive assembly accepted `1073` central blocks.
The least outward-rounded curvature margin is `1.14426650120038260530285332880335545546623246507456181669743E-10`;
the largest certified value of `t^2 U(t)` is `3.33723355756790554314260494045514876254074160915177414272085E+0`.

| mode interval | t ball | t^2 U upper | margin lower | cover tiles |
|---|---:|---:|---:|---:|
| `1159/1250..4641/5000` | `[3.2E+2 +/- 1.02]` | `2.39844033492698029750497377767239428573148533521454389019247E+0` | `1.07306486151606984577502451727104682892134069100945360183724E-5` | `388` |
| `747/625..5981/5000` | `[1.03E+3 +/- 6.71]` | `2.75524084335745274134492233869316795198934192949499418449836E+0` | `6.92961603676221631411987885253995641447762868909901636494613E-7` | `188` |
| `1829/1250..7321/5000` | `[3.4E+3 +/- 23.7]` | `2.99107519293676862219192691364207014762671809785681614692579E+0` | `4.34195440676075616044439759032864885901566099603725463632512E-8` | `128` |
| `1082/625..8661/5000` | `[1.14E+4 +/- 44.2]` | `3.15119995133655700957548028695420910470362206099252469003730E+0` | `2.68061966135575904580172937594460461936302501433289927885378E-9` | `108` |
| `9991/5000..2499/1250` | `[3.8E+4 +/- 4.53E+2]` | `3.33723355756790554314260494045514876254074160915177414272085E+0` | `1.14426650120038260530285332880335545546623246507456181669743E-10` | `104` |
| `2499/1250..2` | `[3.78E+4 +/- 84.7]` | `3.31827351870434774114528181168955445530573825506827225632098E+0` | `1.26846550560873917809914608906872361111667575422741488281799E-10` | `84` |

## Remaining Ray

The sole remaining part of the first-summand curvature theorem is

```text
K_1(t)<=7/(2*t^2) on the mode ray u>=2.
```

The ray requires analytic higher-cumulant saddle bounds. Finite scouts or
extrapolation from this compact certificate are not promoted to that theorem.

```text
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
outputs/formal_core.md
```
