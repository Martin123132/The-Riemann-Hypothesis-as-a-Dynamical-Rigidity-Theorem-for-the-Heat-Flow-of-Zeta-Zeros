# Jensen-Window PF Order-Four Gaussian Cumulant Ray Target

Date: 2026-07-12

Status: exact formal cumulant algebra, global formal-corridor theorem, and
conditional exact-density handoff. This is not a proof of the exact-density
cumulant corridor or continuum ray,
order-four entry, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_gaussian_cumulant_ray_target`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.py
```

## Exact Expansion

Put `epsilon=q^(-1/2)` and

```text
R_epsilon(y)=sum_(r=3)^8 L_r*epsilon^(r-2)*y^r/r!
```

Exact tilted-Gaussian power-series algebra through `epsilon^6` gives:

```text
kappa_2=1+epsilon^2*((2*L_3**2 - L_4)/2)+epsilon^4*((75*L_3**4 - 109*L_3**2*L_4 + 25*L_3*L_5 + 16*L_4**2 - 3*L_6)/24)+epsilon^6*((720*L_3**6 - 1745*L_3**4*L_4 + 503*L_3**3*L_5 + 951*L_3**2*L_4**2 - 103*L_3**2*L_6 - 313*L_3*L_4*L_5 + 14*L_3*L_7 - 66*L_4**3 + 25*L_4*L_6 + 15*L_5**2 - L_8)/48)
```

```text
kappa_3=epsilon^1*(-L_3)+epsilon^3*(-(8*L_3**3 - 7*L_3*L_4 + L_5)/2)+epsilon^5*(-(525*L_3**5 - 954*L_3**3*L_4 + 234*L_3**2*L_5 + 298*L_3*L_4**2 - 37*L_3*L_6 - 57*L_4*L_5 + 3*L_7)/24)
```

```text
kappa_4=epsilon^2*(3*L_3**2 - L_4)+epsilon^4*((48*L_3**4 - 59*L_3**2*L_4 + 11*L_3*L_5 + 7*L_4**2 - L_6)/2)+epsilon^6*((4725*L_3**6 - 10257*L_3**4*L_4 + 2592*L_3**3*L_5 + 4948*L_3**2*L_4**2 - 456*L_3**2*L_6 - 1406*L_3*L_4*L_5 + 52*L_3*L_7 - 298*L_4**3 + 94*L_4*L_6 + 57*L_5**2 - 3*L_8)/24)
```

```text
kappa_5=epsilon^3*(-15*L_3**3 + 10*L_3*L_4 - L_5)+epsilon^5*(-(384*L_3**5 - 605*L_3**3*L_4 + 125*L_3**2*L_5 + 160*L_3*L_4**2 - 16*L_3*L_6 - 25*L_4*L_5 + L_7)/2)
```

```text
kappa_6=epsilon^4*(105*L_3**4 - 105*L_3**2*L_4 + 15*L_3*L_5 + 10*L_4**2 - L_6)+epsilon^6*((3840*L_3**6 - 7365*L_3**4*L_4 + 1605*L_3**3*L_5 + 3095*L_3**2*L_4**2 - 237*L_3**2*L_6 - 745*L_3*L_4*L_5 + 22*L_3*L_7 - 160*L_4**3 + 41*L_4*L_6 + 25*L_5**2 - L_8)/2)
```

```text
kappa_7=epsilon^5*(-945*L_3**5 + 1260*L_3**3*L_4 - 210*L_3**2*L_5 - 280*L_3*L_4**2 + 21*L_3*L_6 + 35*L_4*L_5 - L_7)
```

```text
kappa_8=epsilon^6*(10395*L_3**6 - 17325*L_3**4*L_4 + 3150*L_3**3*L_5 + 6300*L_3**2*L_4**2 - 378*L_3**2*L_6 - 1260*L_3*L_4*L_5 + 28*L_3*L_7 - 280*L_4**3 + 56*L_4*L_6 + 35*L_5**2 - L_8)
```

At `L_3=...=L_8=1`, the leading terms are

```text
kappa_3~-epsilon, kappa_4~2*epsilon^2, kappa_5~-6*epsilon^3,
kappa_6~24*epsilon^4, kappa_7~-120*epsilon^5, kappa_8~720*epsilon^6.
```

## Candidate Corridor

The exact remainder theorem should prove

```text
2/5<=q*(kappa_2-1)<=4/5
1<=(-1)^3*kappa_3*q^(3/2-1)/1<=6/5
1<=(-1)^4*kappa_4*q^(4/2-1)/2<=27/20
1<=(-1)^5*kappa_5*q^(5/2-1)/6<=3/2
1<=(-1)^6*kappa_6*q^(6/2-1)/24<=17/10
1<=(-1)^7*kappa_7*q^(7/2-1)/120<=2
1<=(-1)^8*kappa_8*q^(8/2-1)/720<=5/2
```

These are theorem targets, not inferred from the formal series.

## Proved Formal Corridor

A `1800000`-block Arb certificate proves
all seven corridors for the exact epsilon-six formal polynomial on
`2<=u<=20`. Coefficient-positive analytic gates prove the same formal
corridors on every `u>=20`. Thus the formal model is closed for all
`u>=2`; only exact-minus-formal density errors remain.

## Conditional Collar Test

Assuming the corridor throughout each displayed mode collar, Arb boxes
cover the full `t+-2` derivative interval and give:

| u | central t | t^2 U upper | margin lower |
|---:|---:|---:|---:|
| `2` | `[37850.322210211384816295110892561324 +/- 1.32E-31]` | `3.41302921388726922813443269598269686597541191154695196523493E+0` | `6.07063657023823245605648341446970652867490597354437089805177E-11` |
| `5/2` | `[346604.16551917153573847482751638528 +/- 1.94E-30]` | `3.31360101108560140779560127274422914253548977481090957558863E+0` | `1.55158651583832356104323403295129732061609696339222988815317E-12` |
| `3` | `[3068741.5423332891027180035054560338 +/- 3.25E-29]` | `3.24973124054629564769390390736331501436994117603607551701168E+0` | `2.65757797490686725761044598198893708225808210272079511030523E-14` |
| `4` | `[223333897.74171562139160708003346697 +/- 2.54E-27]` | `3.17762683943435427900689759879005068582321250837335334353054E+0` | `6.46322962911638326154619636030675076415798969659142619932240E-18` |
| `5` | `[15241916613.768536123445468287824895 +/- 5.29E-26]` | `3.13757496136193659055582223171027832667951598844453811107410E+0` | `1.56005184413318624165544962642446145087844312913497145735256E-21` |
| `10` | `[14789692501169113422.216500789131780 +/- 4.94E-17]` | `3.06434436962938135359099927303431706838620575446388426690949E+0` | `1.99170513108481826117410246346200212181768327622493573090999E-39` |
| `20` | `[6.9625514316503258644917213381661481E+36 +/- 48.8]` | `3.03108040967912469979483304505520783669035696111910588605874E+0` | `9.67300794136377874241480509780005437057106007231226713009237E-75` |

This is a conditional finite compatibility check, not a continuum ray
certificate. It shows that the proposed corridor constants retain enough
margin for the localized order-four inequality.

## Remaining Theorem

Prove an explicit central-window remainder after subtracting the exact
`epsilon^6` polynomial, control both adaptive tails, and compose those
bounds with the exact normalized-potential geometry for every `u>=2`.
That theorem would close the sole remaining curvature ray.

```text
outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md
outputs/formal_core.md
```
