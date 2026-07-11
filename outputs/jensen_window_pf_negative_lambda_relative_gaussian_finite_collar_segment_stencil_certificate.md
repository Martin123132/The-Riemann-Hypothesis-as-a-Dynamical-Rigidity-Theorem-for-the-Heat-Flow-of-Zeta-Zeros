# Jensen-Window PF Negative-Lambda Relative-Gaussian Finite-Collar-Segment Stencil Certificate

Date: 2026-07-10

Status: bounded real-T fixed-k stencil certificate. This is not a proof
of the T>10000 ray, an all-k cone-entry theorem, RH, or `Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.json
```

Generator and checker:

```text
work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py
work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py
```

## Result

A two-regime exact-moment/Cauchy and incomplete-Gamma majorant proves the rational residual budgets uniformly for every real 1156<=T<=10000 and i=21..24. Composed with the Arb degree-40 perturbation ledger, this certifies the complete fixed-k=22 stencil system on the whole bounded collar segment. The unbounded ray T>10000, remaining k coverage, cone entry, RH, and Lambda <= 0 remain open.

The exact split is

```text
For 1156<=T<=15625/8 use x_*=8/25 and y_*=(64/625)T; for 15625/8<=T<=10000 use y_*=200 and x_*=sqrt(200/T).
```

Both regimes keep the compact core inside `x<=8/25<0.38`. The
compact Taylor polynomial is bounded by absolute full-Gamma moments,
the degree-384 remainder by a convergent Cauchy estimate, and the real
tails by endpoint-monotone finite-Phi majorants and upper
incomplete-Gamma moments. The global `n>=31` normalization correction
is then scaled at the worst finite endpoint `T=10000`.

## Uniform Residual Bounds

| Index | value scaled upper | fraction of 1/2 | derivative scaled upper | fraction of 9/1000 |
|---:|---:|---:|---:|---:|
| `F_21` | `6.757730132439229252683894686861241721083757336464753075017224062845976E-5` | `1.351546026487845850536778937372248344216751467292950615003444812569196E-4` | `1.248273417045932066294304752284783362926272786957385142815939671546479E-6` | `1.386970463384368962549227502538648181029191985508205714239932968384976E-4` |
| `F_22` | `1.346937511803755598508489832952951804540579863244450427407910463138745E-4` | `2.693875023607511197016979665905903609081159726488900854815820926277488E-4` | `2.489317755485167307932907581001963810556226986583532840862268477509947E-6` | `2.765908617205741453258786201113293122840252207315036489846964975011052E-4` |
| `F_23` | `2.625893479429003530173742599970272409704145957427375239669469031955057E-4` | `5.251786958858007060347485199940544819408291914854750479338938063910113E-4` | `4.855534549517396470369985090802091381171870135895810337592917020195589E-6` | `5.395038388352662744855538989780101534635411262106455930658796689106210E-4` |
| `F_24` | `5.014403071407194264651499379491189118436495197962007521139249792474134E-4` | `1.002880614281438852930299875898237823687299039592401504227849958494828E-3` | `9.277037504604688731663878157658373206814452285080325875479148766706465E-6` | `1.030781944956076525740430906406485911868272476120036208386572085189608E-3` |

All four rows therefore satisfy, for every real `1156<=T<=10000`,

```text
|R_i(1/T)|  <= (1/2) T^(-3)
|R_i'(1/T)| <= (9/1000) T^(-1).
```

## Stencil Composition

| Target | finite margin lower | perturbation upper | retained lower |
|---|---:|---:|---:|
| `normalizer` | `4.396279349275512583580856491807389999989460662793632652531137141158390E-1` | `3.236655688068670047095463507530522232789888658630038509677590150027123E-10` | `4.396279346038856895512186444711926492458938430003743993901098631480800E-1` |
| `B_product` | `3.720335729036512041561022096097329999989652433846591632890398864647758E+1` | `1.730103806508361218691061422759123140789837563390128776698100217100139E-3` | `3.720162718655861205439152989955054087675573450090252620012729054626048E+1` |
| `companion_product` | `1.722502237565501330236724575177899999993340738612960792171015944016988E+1` | `4.000000001941993413260239630015476153089032263105205690220702446380269E+0` | `1.322502237371301988910700612176352384684437512302440223148945699378961E+1` |
| `weighted_gap_derivative` | `2.742440242696800557239471524916489999987929124772433779765180231484331E+1` | `3.691040397770267727819495342225868845483873301706132333444232206394856E+0` | `2.373336202919773784457521990693903115439541794601820546420757010844846E+1` |

## Proof Boundary

Uniform fixed-k=22 residual and stencil certificate on the bounded real interval 1156<=T<=10000 only. It does not cover T>10000, does not cover all k, does not prove the full cone-entry or sign-regularity theorem, and does not prove RH or Lambda <= 0.

The numerical-integration and between-grid-T obligations are now closed
for this fixed k on the bounded segment. The next T obligation is the
unbounded ray `T>10000`; the separate all-k obligation also remains.

## Reproduction

```powershell
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian finite-collar-segment stencil certificate: 10 rows, 0 issues, 4 uniform residual rows, 2 T regimes, 4 positive perturbation margins, 0 open finite-segment stencil sources, 0 ready-to-apply rows
```
