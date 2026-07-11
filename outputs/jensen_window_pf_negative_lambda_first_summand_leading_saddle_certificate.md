# Jensen-Window PF Negative-Lambda First-Summand Leading-Saddle Certificate

Date: 2026-07-10

Status: interval leading-saddle theorem and open remainder target.
This is not a proof of the uniform remainder, cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: 8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows
```

## Log-Variable Potential

Put `x=2*log(u)` and `q=pi*exp(4u)`. The continuous first-summand
density is proportional to `exp(t*x-V(x))`, where, up to a constant,

```text
V(x)=100*u^2+q-5*u-log(2*q-3)-log(u).
```

At its mode `x_t`, write `a_t=V''(x_t)` and `b_t=V'''(x_t)`.
The exact original-u slope at `t=318,u=0.926` is positive:

```text
[0.35965164775903244596420837590841368229664609110202471 +/- 3.50E-54]
```

Since the original-u slope decreases in `u` and increases in `t`, every
x-density mode for `t>=318` lies on `u>=0.926`.

## Compact Certificate

Direct Arb evaluation of the exact rational-exponential formulas proves

```text
t(u)^2*V'''(u)/V''(u)^3<=13/20
```

on `40740/40740`
subintervals covering `0.926<=u<=5`. The same intervals also prove
the cubic Edgeworth correction is at most `1/100` and the fifth-order
correction is at most `1/1000`.

- Minimum outward-rounded curvature lower bound: `6.388903512954711914062499999999999999999999999999999999999999999999999E+2`.
- Minimum outward-rounded cap margin: `2.585626076543121598660945892333984374999999999999999999993627632355469E-2`.
- Minimum outward-rounded correction-cap margin: `5.536307266169015212917965912612872756575262075326505678724302012973085E-3`.
- Minimum outward-rounded fifth-correction margin: `9.048398819246355472550973349878408403626072605579984548297001637844215E-4`.

## Analytic Ray

On `u>=5`, put `r=2/(2q-3)`. The Arb endpoint gates and exact
monotonicity give `r<=1/100` and `q/u` increasing. Differentiating the
potential componentwise yields

```text
V'<=3*u*q,
V''>=(39/10)*u^2*q,
V'''<=12*u^3*q.
```

Therefore

```text
t^2*V'''/V''^3<=108000/(59319*u)<=108000/(59319*5)<13/20.
```

The final exact ray margin is `12561/43940`.
The derivative formulas through order seven and `q>=10^9` also give
the correction-cap margin `390971184173/39097152900000`.
The Bell-polynomial bounds through order seven give the fifth-correction
margin `9276816051497935508833/9276816051500400000000000`.

## Remaining Remainder

The cumulant bridge requires `t^2*kappa_3>=-37/50`. The leading saddle
term is certified below `13/20`, the cubic correction below `1/100`,
and the fifth-order correction below `1/1000`. It is therefore enough
to prove the seventh-order remainder bound

```text
scaled cumulant + leading + cubic correction + fifth correction
>=-79/1000, t>=318.
```

This signed uniform Laplace/Stein remainder is the sole open row in this
artifact. The recorded finite samples are:

| t | third-order remainder | fifth correction | seventh remainder | margin above -0.079 |
|---:|---:|---:|---:|---:|
| `318` | `-0.0000047900872336456278227269922572717659` | `0.0000048007568691244615813561506044316083` | `0.000000010669635478833758629158347159842314` | `0.07900001066963547883375862915834716` |
| `319` | `-0.0000047678369532683959850188022806467284` | `0.0000047784885418565753803741992695042559` | `0.000000010651588588179395355396988857527435` | `0.079000010651588588179395355396988858` |
| `400` | `-0.0000034420412636514221577766141243393302` | `0.0000034505922757872871601213531569897609` | `0.0000000085510121358650023447390326504307298` | `0.07900000855101213586500234473903265` |
| `700` | `-0.0000015569619706331105438661163589860474` | `0.0000015599826445992151024577034374113994` | `0.0000000030206739661045585915870784253520386` | `0.079000003020673966104558591587078425` |
| `1000` | `-0.00000091054494626537261565643491859796169` | `0.00000091182304300042701145611735851382065` | `0.0000000012780967350543957996824399158589598` | `0.079000001278096735054395799682439916` |
| `2000` | `-0.00000029469085446147090488299852334526449` | `0.00000029488640355029928141255426003622294` | `0.00000000019554908882837652955573669095844862` | `0.079000000195549088828376529555736691` |
| `5000` | `-0.000000058850333712340444976622768928994421` | `0.00000005886407793589575318498515966192756` | `1.3744223555308208362390732933138591e-11` | `0.079000000013744223555308208362390733` |
| `20000` | `-0.0000000045624608990187173726924209496064824` | `0.0000000045626916401450401892544559468853232` | `2.3074112632281656203499727884083508e-13` | `0.079000000000230741126322816562034997` |
| `100000` | `-0.00000000021898084301318233326642455403256557` | `0.00000000021898290763151665516268915321336152` | `2.0646183343218962645991807959534686e-15` | `0.079000000000002064618334321896264599` |

These rows inherit non-interval cumulant quadrature and are not a uniform
proof. They indicate that the proposed remainder floor has substantial room.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The leading saddle contribution is now certified globally: the leading cap is 13/20, the cubic correction cap is 1/100, and the fifth-order correction cap is 1/1000 for all t>=318. This reduces the cumulant wall to the seventh-order normalized remainder floor -79/1000. Nine finite samples clear that floor; the smallest sampled margin is 0.079000000000002064618334321896264599 at t=100000. The uniform remainder theorem remains open.
