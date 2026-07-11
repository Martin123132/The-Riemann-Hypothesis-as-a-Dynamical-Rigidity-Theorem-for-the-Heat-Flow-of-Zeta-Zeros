# Jensen-Window PF Negative-Lambda First-Summand Cumulant Bridge

Date: 2026-07-10

Status: exact cumulant reduction with analytic hypothesis discharged; not a proof of the remaining all-order bridge.
This proves the lambda=-100 adjacent wall, not the remaining all-order bridge, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_first_summand_cumulant_bridge`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda first-summand cumulant bridge: 8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, 0 open requirements, 3 ready-to-apply rows
```

## Exact Decomposition

For real `t>=318`, define

```text
M_t^(1)=2*integral_0^infinity u^(2t)*exp(-100*u^2)*phi_1(u)du,
F(t)=log(M_t^(1)).
```

Under the probability measure proportional to the integrand, differentiation
under the integral gives

```text
F'''(t)=kappa_3,t(2*log(U)).
```

Three applications of the fundamental theorem of calculus give

```text
Delta^3 F(k-1)=integral_[0,1]^3 F'''(k-1+r+s+v)drdsdv.
```

The Gamma normalization is exact:

```text
-Delta^3 log Gamma(k-1/2)=-log(1-1/(k+1/2)^2).
```

## Sufficient Cumulant Bound

The paired compact and ray certificates prove the single continuous estimate

```text
kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318.
```

Indeed, `-log(1-z)>=z` and the cubic-difference integral imply

```text
L_k^(1)>=1/(k+1/2)^2-37/(50*(k-1)^2).
```

After subtracting `1/(4*k^2)`, the exact rational numerator is

```text
4*k^4-996*k^3+401*k^2-50*k-25.
```

With `n=k-319`, this becomes

```text
4*n^4+4108*n^3+1489493*n^2+215582064*n+9130082706.
```

Every coefficient is positive for `n>=0`. Therefore the cumulant estimate
implies `L_k^(1)>=1/(4*k^2)` for every integer `k>=319`.

## High-Precision Scout

The following 60-digit saddle-centered calculations use finite upper
truncation and are not interval enclosures.

| t | saddle | t^2 kappa_3 | margin above -37/50 |
|---:|---:|---:|---:|
| `318` | `0.92612054479248795264653541275789266` | `-0.62398825392382784408826296832121367` | `0.11601174607617215591173703167878633` |
| `319` | `0.92684312642253437374664751576064431` | `-0.62345680091994747111173551642919344` | `0.11654319908005252888826448357080656` |
| `400` | `0.97888094701689181843228254728997258` | `-0.58619255414281535730902534503577599` | `0.15380744585718464269097465496422401` |
| `700` | `1.1069452994078673013038750783415861` | `-0.50456881698434643552227489961979525` | `0.23543118301565356447772510038020475` |
| `1000` | `1.1879375411452161946579555446887138` | `-0.46099339683943611590699232273178012` | `0.27900660316056388409300767726821988` |
| `2000` | `1.3440078260003367251495307467509845` | `-0.39329116138832195333736837391724843` | `0.34670883861167804666263162608275157` |
| `5000` | `1.5485738239457640734903380663464682` | `-0.3299410842432296615175601567087721` | `0.4100589157567703384824398432912279` |
| `20000` | `1.8573688080533650862594433804397783` | `-0.26808362081845586697830696908227599` | `0.47191637918154413302169303091772401` |
| `100000` | `2.2183616301357350926354408275557385` | `-0.22253562884308290948889082506809256` | `0.51746437115691709051110917493190744` |

All nine recorded samples exceed the proposed floor, and the sampled
scaled cumulant is strictly increasing. This is theorem-search evidence only.

## Proof Boundary

The exact algebra reduces the existing first-summand wall target to a
special-kernel skewness estimate. Generic strict log-concavity is not enough
to determine third-cumulant sign. The paired remainder theorem supplies
the required explicit `q=pi*exp(4u)` saddle geometry.

With the cumulant estimate proved, the existing dominance certificate
and collar extension give

```text
L_k>=1/(4*k^2)-16/(k-1)^6>0, k>=319,
```

spliced to the certified finite adjacent wall through `k=318`.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md
outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The first-summand wall is exactly reduced to the uniform bound kappa_3,t(2 log U)>=-37/(50 t^2) for t>=318. The Gamma contribution and rational transfer are exact; nine finite high-precision samples clear the floor, with the smallest sampled margin 0.11601174607617215591173703167878633 at t=318. The uniform cumulant estimate is discharged by the paired ray certificate.
