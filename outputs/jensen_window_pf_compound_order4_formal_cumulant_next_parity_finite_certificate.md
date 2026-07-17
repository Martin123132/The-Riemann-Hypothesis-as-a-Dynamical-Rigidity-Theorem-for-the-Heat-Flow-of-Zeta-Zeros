# Jensen-Window PF Order-Four Next-Parity Finite Certificate

Date: 2026-07-13

Status: rigorous finite interval theorem for the next-parity coefficient
functions. This is not a proof of the exact cumulant ray, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.json
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_chunks.jsonl
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.py
```

## Certified Bounds

A deterministic cover of `1800` adjacent
blocks of width `1/100` proves:

```text
-1/20 < C_2(u) < 0
-3/20 < C_3(u) < 0
-1/2 < C_4(u) < 0
0 < C_5(u) < 2
0 < C_6(u) < 16/5
0 < C_7(u) < 5
0 < C_8(u) < 61/10
```

Here `C_r` is the exact rational polynomial multiplying the first
omitted corridor-scaled power: `q^-3` for orders 2-4, `q^-2` for
orders 5-6, and `q^-1` for orders 7-8.

| order | certified lower | certified upper | lower margin | upper margin |
|---:|---:|---:|---:|---:|
| 2 | `-3.39429863998777575679748243861266758754525960720034340036451E-2` | `-2.58463410263116604598357789776504008619693367014130683855182E-2` | `1.60570135992127477302522473759582850620474039279965659963469E-2` | `2.58463410263116604598357789776504008619693367014130683855182E-2` |
| 3 | `-1.38180852002532073790585668489933319896647521290626811915571E-1` | `-1.24658520179043711117136803466266649465971826778600483516544E-2` | `1.18191479938299474023226185584065238533524787093731880844058E-2` | `1.24658520179043711117136803466266649465971826778600483516544E-2` |
| 4 | `-3.79255765437641569926259511302979952148790411595134710566324E-1` | `-3.03079558497495369642259370707919425805726772470307049473285E-2` | `1.20744234562358430073740488697020047851209588404865289433676E-1` | `3.03079558497495369642259370707919425805726772470307049473285E-2` |
| 5 | `1.70294120590486763540422951869562899998519832467537187705910E+0` | `1.99934731554569624376338980783983506651019266787515148680853E+0` | `1.70294120590486763540422951869562899998519832467537187705910E+0` | `6.52684454303756236610192160164933489807332124848513191482203E-4` |
| 6 | `2.56394571293543115820216732652909175940494946583078625013066E+0` | `3.11093222683365011620270397001106435369076540069418623219137E+0` | `2.56394571293543115820216732652909175940494946583078625013066E+0` | `8.90677731372460533405623263756543963092345993058137678081253E-2` |
| 7 | `3.17258059153290829781968974681191203652859812436964061218821E+0` | `4.89516142943693284165374604508445163223169807112034404094818E+0` | `3.17258059153290829781968974681191203652859812436964061218821E+0` | `1.04838570563067158346253954915548367768301928879655959051832E-1` |
| 8 | `3.72415667605322432810433840772647169096617023253853623103600E+0` | `6.03235588312468070770696521725972747199954878200262616271682E+0` | `3.72415667605322432810433840772647169096617023253853623103600E+0` | `6.76441168607673770646679309336319030004512179973738372829316E-2` |

## Taylor Gate

Direct interval substitution loses the cancellations in the 15-22-term
coefficient polynomials. Each block instead uses a centered sixth-order
Arb Taylor model. Exact rational-center coefficients enclose derivatives
through order six, while a full-block natural interval extension encloses
the seventh derivative remainder. Thus the range enclosure is rigorous,
not a point sample or floating-point fit.

The weakest outward-rounded margin is `0.000652684454303756236610192160164933489807332124848513191482203`
for order `5` on the
`minimum_upper_margin` side.

## Remaining Boundary

This closes only the explicit next-parity coefficient layer on the finite
interval. The asymptotic coefficient ray, all terms beyond epsilon eight,
and the exact-density central and two-tail estimates remain open.

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.md
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/formal_core.md
```
