# Jensen-Window PF Compound Order-Five Nested Curvature Compact Certificate

Date: 2026-07-13

Status: rigorous compact first-summand order-five curvature theorem
with one open analytic ray. This is not a proof of full order-five
entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py
```

## Interval Core

The two stable centered differences require a common three-unit
collar. On each central block, a single outward-rounded hull for
`H^(2),...,H^(8)` encloses

```text
B^(r)=integral_[-1,1](1-|s|)*H^(r+2)(t+s) ds
J^(r)=2*B^(r)-integral_[-1,1](1-|s|)*ell^(r+2)(t+s) ds
R^(r)=3*B^(r)-integral_[-1,1](1-|s|)*h^(r+2)(t+s) ds.
```

Truncated Arb Taylor-series arithmetic then evaluates
`ell=log(1-exp(-B))`, `h=2ell+log(1-exp(-J))`, and
`q=2h-ell+log(1-exp(-R))` through second order. These are finite
jet identities, not asymptotic series.

## Cached Cover

The source cache has `107452`
paired interval-Simpson tiles and SHA-256

```text
d721a22738543dd2f62181a31732b26666d13eb3f8c4f1c8c46ead3e84ada4cf
```

Adaptive assembly accepts `36` central blocks.
The first block straddles `t=320` and the last reaches
the exact mode `u=2`, so the certified range is

```text
320<=t<=V'(2)
```

Every block proves `J_1>0`, `R_1>0`, and a positive curvature
margin. The global diagnostics are

```text
weakest J lower: 4.09966568242201277102605766062470486046577941353033990939654E-5
weakest R lower: 6.14928302823389707801310674634727092894371032168610427734783E-5
largest t^2*q_1'' upper: 2.20266810828795017766530766838057878229544402251467398253894E+1
worst absolute margin: 3.13630306771174955427539380510717568375832632627546459216436E-8
```

Thus the compact theorem is

```text
q_1''(t)<=60/t^2 for every 320<=t<=V'(2).
```

## Remaining Ray

The only unproved part of the continuous first-summand target is

```text
q_1''(t)<=60/t^2 for every mode u>=2.
```

The exact-cumulant corridor cover on `2<=u<=20` and normalized-H
boxes on `u>=20` are the natural inputs for that final ray theorem.

| mode interval | t ball | t^2 q_1'' upper | margin lower | J lower | R lower |
|---:|---:|---:|---:|---:|---:|
| `11599/12500..5987/6250` | `[3E+2 +/- 65.1]` | `2.20266810828795017766530766838057878229544402251467398253894E+1` | `2.85025207342478162459272846833012496001586301802855615359602E-4` | `2.81733602512301863803341685049838034733768313383207848289866E-3` | `4.20636358306620988821211175766339821826087252157790200088841E-3` |
| `7487/6250..15349/12500` | `[1.1E+3 +/- 93.5]` | `1.61515405884962882332351923511350079499659118054884204687367E+1` | `3.07870967165756979520602085340039776571407586179301981686562E-5` | `1.05585171956506770805130361522794875407719738194261005138013E-3` | `1.58167405864647460651298678241051345710913125810103402065787E-3` |
| `18349/12500..4681/3125` | `[4E+3 +/- 5.24E+2]` | `1.65170551138358531420285723967286746875628140317935788067864E+1` | `2.74078400560736068923013856643734969592229551433392013722668E-6` | `3.52290112677541359124970683126711003359605388256286918390068E-4` | `5.28226758408282168053962610113712811381480117703871777898456E-4` |
| `5431/3125..22099/12500` | `[1E+4 +/- 3.40E+3]` | `1.73503590231677911480461647835865561400579038427555792045952E+1` | `2.37826253674774376827817518893691557244087926117473943476531E-7` | `1.11663159358678015952191624614409533572853245095679858654336E-4` | `1.67475142280039099963218452571861633234706524271639940414678E-4` |
| `6181/3125..2` | `[4E+4 +/- 5.73E+3]` | `1.50668024436921805086196165985870423964804030460251616539728E+1` | `3.13630306771174955427539380510717568375832632627546459216436E-8` | `4.09966568242201277102605766062470486046577941353033990939654E-5` | `6.14928302823389707801310674634727092894371032168610427734783E-5` |

```text
outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
