# Jensen-Window PF Negative-Lambda Relative-Gaussian Recorded-Grid Stencil Composition Certificate

Date: 2026-07-10

Status: recorded-grid fixed-k stencil composition certificate. This is not a proof
of a real-T collar, an all-k cone-entry theorem, RH, or `Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.json
```

Generator and checker:

```text
work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py
work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py
```

## Result

The 20 direct expectation balls satisfy the exact rational residual budgets A=1/2 and B=9/1000. An independent Arb perturbation ledger preserves all degree-40 fixed-k=22 normalizer, B-product, companion-product, and weighted-gap derivative margins. Therefore the complete fixed-k stencil system is certified at all five recorded T values. The real-T collar, remaining k coverage, cone entry, RH, and Lambda <= 0 remain open.

The exact residual identities are

```text
E[value core]      = R_i(u)
E[derivative core] = u R_i'(u).
```

The composition uses the exact rational sufficient budgets

```text
|R_i(u)|  <= (1/2) u^3
|R_i'(u)| <= (9/1000) u.
```

The older floating threshold decimals are recorded for comparison but are
not proof inputs.

## Perturbation Ledger

| Target | finite margin lower | perturbation upper | retained lower |
|---|---:|---:|---:|
| `normalizer` | `4.396279349275512583580856491807389999989460662793632652531137141158390E-1` | `3.236655688068670047095463507530522232789888658630038509677590150027123E-10` | `4.396279346038856895512186444711926492458938430003743993901098631480800E-1` |
| `B_product` | `3.720335729036512041561022096097329999989652433846591632890398864647758E+1` | `1.730103806508361218691061422759123140789837563390128776698100217100139E-3` | `3.720162718655861205439152989955054087675573450090252620012729054626048E+1` |
| `companion_product` | `1.722502237565501330236724575177899999993340738612960792171015944016988E+1` | `4.000000001941993413260239630015476153089032263105205690220702446380269E+0` | `1.322502237371301988910700612176352384684437512302440223148945699378961E+1` |
| `weighted_gap_derivative` | `2.742440242696800557239471524916489999987929124772433779765180231484331E+1` | `3.691040397770267727819495342225868845483873301706132333444232206394856E+0` | `2.373336202919773784457521990693903115439541794601820546420757010844846E+1` |

## Recorded T Systems

| T | max |R|/u^3 | max |R'|/u | normalizers | B | companion | weighted gap |
|---:|---:|---:|---:|---:|---:|---:|
| `1156` | `2.824904766903341294024455616699704095022017472896465328791038018607652E-4` | `5.078764204556721600091363717104104023901747777536286321264924481511117E-6` | `True` | `True` | `True` | `True` |
| `1500` | `2.733853646250007555792814628156111946056257657518967496557138474599925E-6` | `3.795342750000097379509041206678641305893862067932786885648965835571290E-8` | `True` | `True` | `True` | `True` |
| `2000` | `1.612358888000022189044739478304132617779653941295081587270487943897025E-8` | `1.681837720000188846348459460411259598481593935392908178982906974852086E-10` | `True` | `True` | `True` | `True` |
| `5000` | `1.208568929584068358253054880163250000110886273803224166719542829285257E-15` | `5.061340683167061023023142200535175000032972700245449492257149638686311E-18` | `True` | `True` | `True` | `True` |
| `10000` | `4.754463729209350530609360800000176159143010136612210436986238842757647E-21` | `9.969502254878374396517988980000030541130693227476141061665009067800048E-24` | `True` | `True` | `True` | `True` |

Maximum fraction of the value budget: `5.649809533806682588048911233399408190044034945792930657582076037215303E-4`.
Maximum fraction of the derivative budget: `5.643071338396357333434848574560115582113053086151429245849916090567908E-4`.

## Proof Boundary

Complete fixed-k=22 stencil composition certificate at T=1156,1500,2000,5000,10000 only. It does not certify an interval in T, does not cover all k, does not prove the full cone-entry theorem or sign-regularity bridge, and does not prove RH or Lambda <= 0.

This closes the fixed-k stencil composition source at the five recorded
T values. It does not interpolate between them and does not extend the
single k-window to all k.

## Reproduction

```powershell
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian recorded-grid stencil composition certificate: 8 rows, 0 issues, 20 residual rows, 4 positive perturbation margins, 5 certified T systems, 0 open recorded-grid stencil sources, 0 ready-to-apply rows
```
