# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Collar Scan

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof
of a real-`T` collar theorem, uniform Taylor-tail theorem,
scaled-curvature monotonicity, cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan`.

Proof boundary: this artifact scans an integer `T` collar for `k=22`
and degree-16 continuation. It does not prove an analytic collar.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: 7 rows, 0 issues, 1301 scan rows, 1045 continuation-positive rows, 718 half-safety rows, 0 ready-to-apply rows
```

## Collar Thresholds

```text
scan rows: 1301
baseline-positive rows: 1283
continuation-positive rows: 1045
half-safety rows: 718
pointwise budget successes: 0
pointwise budget failures: 1283
first baseline-positive T: 918
first continuation-positive T: 1156
first continuation q/T: 1.946366782006920415E-2
first half-safety T: 1483
first half-safety q/T: 1.517194875252865813E-2
worst pointwise over-budget: 2.638097770251728679E+5 at T=918
worst stencil abs over half-margin: 7.209437565352806308E+2 at T=918
```

Selected scan rows:

```text
T=900: baseline=False, continuation=False, half_safety=False, companion_M8=-2.585087157631684242E-6, weighted_gap_M8=5.383536975387089451E-4
T=950: baseline=True, continuation=False, half_safety=False, companion_M8=-1.355197819587097835E-6, weighted_gap_M8=4.352810148466106944E-4
T=1000: baseline=True, continuation=False, half_safety=False, companion_M8=-6.930385931843005222E-7, weighted_gap_M8=3.672349168521610482E-4
T=1100: baseline=True, continuation=False, half_safety=False, companion_M8=-1.238943873083390585E-7, weighted_gap_M8=2.831230178557149462E-4
T=1153: baseline=True, continuation=False, half_safety=False, companion_M8=-3.950543371415399732E-9, weighted_gap_M8=2.533346006074945900E-4
T=1154: baseline=True, continuation=False, half_safety=False, companion_M8=-2.319815056389729095E-9, weighted_gap_M8=2.528349616498292444E-4
T=1155: baseline=True, continuation=False, half_safety=False, companion_M8=-7.079873820087283923E-10, weighted_gap_M8=2.523373321469098035E-4
T=1156: baseline=True, continuation=True, half_safety=False, companion_M8=8.851471227548859301E-10, weighted_gap_M8=2.518416978813492713E-4
T=1157: baseline=True, continuation=True, half_safety=False, companion_M8=2.459793526702945523E-9, weighted_gap_M8=2.513480447752841385E-4
T=1158: baseline=True, continuation=True, half_safety=False, companion_M8=4.016154526233547062E-9, weighted_gap_M8=2.508563588887468845E-4
T=1200: baseline=True, continuation=True, half_safety=False, companion_M8=5.523957637225961793E-8, weighted_gap_M8=2.318158543664976342E-4
T=1400: baseline=True, continuation=True, half_safety=False, companion_M8=1.167906699005006557E-7, weighted_gap_M8=1.690333664807173359E-4
T=1480: baseline=True, continuation=True, half_safety=False, companion_M8=1.126814727964722694E-7, weighted_gap_M8=1.516368261888900678E-4
T=1481: baseline=True, continuation=True, half_safety=False, companion_M8=1.125997379393203429E-7, weighted_gap_M8=1.514380843190536180E-4
T=1482: baseline=True, continuation=True, half_safety=False, companion_M8=1.125174810934742418E-7, weighted_gap_M8=1.512397628351518590E-4
T=1483: baseline=True, continuation=True, half_safety=True, companion_M8=1.124347089699779384E-7, weighted_gap_M8=1.510418603694484371E-4
T=1484: baseline=True, continuation=True, half_safety=True, companion_M8=1.123514282141585092E-7, weighted_gap_M8=1.508443755609399430E-4
T=1485: baseline=True, continuation=True, half_safety=True, companion_M8=1.122676454062708496E-7, weighted_gap_M8=1.506473070553093278E-4
T=1600: baseline=True, continuation=True, half_safety=True, companion_M8=1.005934937316762764E-7, weighted_gap_M8=1.304529457317573060E-4
T=1800: baseline=True, continuation=True, half_safety=True, companion_M8=7.920154197169507982E-8, weighted_gap_M8=1.041170141427360220E-4
T=2000: baseline=True, continuation=True, half_safety=True, companion_M8=6.191220017299330373E-8, weighted_gap_M8=8.513530052272832770E-5
T=2200: baseline=True, continuation=True, half_safety=True, companion_M8=4.889720300717582147E-8, weighted_gap_M8=7.095007736307137735E-5
```

Interpretation:

The finite collar signal is now explicit: continuation positivity starts
at `T=1156`, while the stricter half-safety condition starts at `T=1483`.
The pointwise budget fails on every baseline-positive row in the scan, so
the live theorem-search target remains a direct signed stencil-tail collar.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The degree-16 collar scan on T=900..2200 makes the large-T signal concrete: the M=7 baseline first has positive stencils at T=918, the M=8 continuation first preserves signs at T=1156 (q/T=1.946366782006920415E-2), and half-safety first holds at T=1483 (q/T=1.517194875252865813E-2). The crude pointwise budget fails on all baseline-positive rows, so the live route remains a direct signed stencil-tail collar theorem.
