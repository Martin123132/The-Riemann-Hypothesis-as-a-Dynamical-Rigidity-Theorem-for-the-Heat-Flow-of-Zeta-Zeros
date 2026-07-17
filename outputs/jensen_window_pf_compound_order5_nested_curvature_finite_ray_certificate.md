# Jensen-Window PF Compound Order-Five Nested Curvature Finite-Ray Certificate

Date: 2026-07-13

Status: rigorous first-summand nested-curvature theorem on
`2<=u<=20` with one open asymptotic ray. This is not a proof of
complete order-five entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py
```

## Mode-Two Collar

The compact cache reaches slightly above `u=2`. One hundred additional
quadrature tiles extend the common `t+-3` derivative collar through
`u=2.001`. Direct potential-jet interval arithmetic proves

```text
abs(V^(8)/V''^4)<1/50000 on 2<=u<=2.002.
```

Thus the same paired Simpson remainder theorem applies. The collar gives

```text
q_1''(t)<=60/t^2 for every 2<=u<=2.001
largest scaled upper=7.93961101108720288769200812148118771867291019912296040070535E+0
margin lower=3.60155806340515842886026326469320324502469625630603569162001E-8
```

The extension cache has SHA-256

```text
3bd7b260926c48682930ad03f2233a2a79c05019fb8e9a04f4c581f3a916be19
```

## Exact-Corridor Cover

The proved seven exact cumulant corridors generate common
`H^(2),...,H^(8)` boxes on each central mode block. Every block has an
explicit outer pad exceeding three units in the continuous `t` variable.
The nested stable interval core then proves `J_1>0`, `R_1>0`, and
the curvature ceiling on all
`1850` blocks from `u=2.001` to `u=20`.

```text
largest t^2*q_1'' upper=11.9132395331211577148958456918358567151542484033635145131552
weakest J lower=2.78700264672900444459507581527921668289773350483977603777114E-37
weakest R lower=4.18050397042193599446874272391657414511023742379948266381521E-37
weakest absolute margin=1.04949542722598669876043610088376084619929188452655330686831E-72
```

Combining the collar and exact-corridor blocks proves

```text
q_1''(t)<=60/t^2 for every mode 2<=u<=20.
```

## Remaining Ray

The sole remaining part of the continuous first-summand theorem is

```text
q_1''(t)<=60/t^2 for every mode u>=20.
```

The existing normalized-H boxes `0<x_r<=1` for `2<=r<=8`, together
with `x_2>=97/100`, `x_3>=24/25`, and `1/t<=10^-30`, leave a very
large interval reserve for that asymptotic composition.

| mode interval | t ball | t^2 q_1'' upper | margin lower | J lower | R lower |
|---:|---:|---:|---:|---:|---:|
| `2001/1000..1003/500` | `[3.8E+4 +/- 8.78E+2]` | `1.02457266152346462652484290485820899290990359632666479908920E+1` | `3.29174749476935091938522193875726025136664465950669653142699E-8` | `3.93739880407802350057955518029163664337314350721982421867563E-5` | `5.90594608586855813723498915359914844538832457436158593465052E-5` |
| `2501/1000..2511/1000` | `[4E+5 +/- 5.19E+4]` | `1.19132395331211577148958456918358567151542484033635145131552E+1` | `3.63400566846094605894100113876886753127266269464693637174989E-10` | `4.42303797891198816158044888026247908254374191955098165540166E-6` | `6.63453717102081309729100956343164834233721114881885890642533E-6` |
| `6121/1000..6131/1000` | `[1.7E+12 +/- 4.71E+10]` | `1.01442222635154113133677602442953726110430524730677566306948E+1` | `1.67887111843854060397658530425544347015929262077537426692661E-23` | `1.06039724604777521594977104719938588850545743148940853818812E-12` | `1.59059586907114139053972781745869070828873767527967768916580E-12` |
| `10751/1000..10761/1000` | `[3.3E+20 +/- 9.36E+18]` | `9.68919269033680723244456827628978642964195879471181675078008E+0` | `4.50871567435833466319568364270535192134393761238020420146235E-40` | `5.67807692451224012512687182859484976532750462869699624143355E-21` | `8.51711538657114496137113325687143923554847207539494373548706E-21` |
| `15371/1000..15381/1000` | `[5E+28 +/- 1.36E+27]` | `9.51672400248364629473221824809338287910175486326844870205933E+0` | `1.96675922284159096447007189988921287472191511096117567232859E-56` | `3.79942586262598921752593808711203908316154655643415232850761E-29` | `5.69913879423285741399447900759962491942088613869887133400780E-29` |
| `19991/1000..20` | `[7E+36 +/- 2.87E+35]` | `9.12346931648153801605124455471636475911533577183637885624613E+0` | `1.04949542722598669876043610088376084619929188452655330686831E-72` | `2.78700264672900444459507581527921668289773350483977603777114E-37` | `4.18050397042193599446874272391657414511023742379948266381521E-37` |

```text
outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
