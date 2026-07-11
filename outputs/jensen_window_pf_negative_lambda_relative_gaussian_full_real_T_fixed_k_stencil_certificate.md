# Jensen-Window PF Negative-Lambda Relative-Gaussian Full-Real-T Fixed-k Stencil Certificate

Date: 2026-07-10

Status: full real-T fixed-k stencil certificate. This is not a proof
of all-k cone entry, sign regularity, RH, or `Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.json
```

Generator and checker:

```text
work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py
work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py
```

## Result

Exact full-kernel evenness gives an order-42 residual zero. A full-kernel disk majorant, factored Cauchy bound, full real-tail majorant, and upper-Gamma hazard monotonicity certify the rational residual budgets on every T>=10000. Combined with the bounded segment theorem, the complete fixed-k=22 stencil system is now certified for every real T>=1156. Remaining k coverage, cone entry, sign regularity, RH, and Lambda <= 0 remain open.

The later interval counterexample
`outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md`
certifies `x_121<x_120` at `T=1156`. It leaves this fixed-`k=22` theorem
untouched, but proves that the theorem cannot be promoted to all-`k` cone
entry at its left endpoint.

The bounded segment certificate covers `1156<=T<=10000`. For the
unbounded ray, fix `x_*=1/5`, so the real-tail split starts at
`y_*=T/25>=400`.

## Full-Kernel Ray Bounds

- Full disk majorant on `|z|<=0.38`: `1.768508446736678439442251752233202656001885187906365988581239520991220E+6`.
- Full real value majorant at `x=1/5`: `1.333508477910873747986465182286269762759788914453547508634664424864289E-1`.
- Full real derivative majorant at `x=1/5`: `9.670711313562588006211702505501499440483482598101461249585525102725776E-1`.
- Value monotonicity margin: `1.896697212486137433573725989867344184211146590172698304273576869902013E+1`.
- Derivative monotonicity margin: `9.966972124861374335737259898673441842111465901726983042735768699020143E+0`.

Exact evenness and subtraction through degree 40 factor the compact
residual by `x^42`. After Gamma integration, the scaled compact value
and derivative bounds decay as `T^-18` and `T^-19`.

For the real tail,

```text
Gamma(s,y)<=y^s*exp(-y)/(y-s+1) for y>s-1; this makes every scaled tail term decrease once y-s+1 exceeds its positive T-power.
```

proves every scaled tail term decreases from `T=10000` onward.

## Ray Residual Rows

| Index | value scaled upper | fraction of 1/2 | derivative scaled upper | fraction of 9/1000 | hazard margin |
|---:|---:|---:|---:|---:|---:|
| `F_21` | `7.376768099446078747003213140698173794698095769516760877883051353223544E-17` | `1.475353619889215749400642628139634758939619153903352175576610270644710E-16` | `1.590103345880599196576248165883828573523811754762501789232568847250409E-19` | `1.766781495422887996195831295426476192804235283069446432480632052500455E-17` | `3.594999999999999999999999999999999999999999999999999999999999999999999E+2` |
| `F_22` | `1.458198345239341147663425853393825052440321256764941103767579918660469E-16` | `2.916396690478682295326851706787650104880642513529882207535159837320937E-16` | `3.143227544182579807185606839537800668593581375693317490343450046890343E-19` | `3.492475049091755341317340932819778520659534861881463878159388940989270E-17` | `3.584999999999999999999999999999999999999999999999999999999999999999999E+2` |
| `F_23` | `2.819183467462726218815956649894728434717954429745552800617321176076905E-16` | `5.638366934925452437631913299789456869435908859491105601234642352153810E-16` | `6.076906585419654293892173223106414625947590659673747147997336757321328E-19` | `6.752118428244060326546859136784905139941767399637496831108151952579254E-17` | `3.574999999999999999999999999999999999999999999999999999999999999999999E+2` |
| `F_24` | `5.338453800088992201587662592353847461487190303135195728828544354698820E-16` | `1.067690760017798440317532518470769492297438060627039145765708870939765E-15` | `1.150733374685849430120007269907384897253905465342475523769708449790635E-18` | `1.278592638539832700133341411008205441393228294824972804188564944211817E-16` | `3.564999999999999999999999999999999999999999999999999999999999999999999E+2` |

## Full T Composition

```text
[1156,10000] union [10000,infinity) = [1156,infinity).
```

The exact rational residual budgets therefore preserve all four
normalizers, the B product, companion product, and weighted-gap
derivative for fixed `k=22` at every real `T>=1156`.

## Proof Boundary

Uniform full real-T fixed-k=22 residual and stencil theorem for every T>=1156. It does not cover all k, does not by itself prove cone entry or sign regularity, and does not prove RH or Lambda <= 0.

The missing cone-entry work has shifted from real-T control to k
coverage: this theorem controls only the window `F_21,...,F_24`.

## Reproduction

```powershell
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian full-real-T fixed-k stencil certificate: 10 rows, 0 issues, 4 ray residual rows, 4 full-kernel n-tail channels, 4 positive perturbation margins, 0 open full-T fixed-k stencil sources, 0 ready-to-apply rows
```
