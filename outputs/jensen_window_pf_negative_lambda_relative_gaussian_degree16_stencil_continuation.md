# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Stencil Continuation

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof
of a uniform Taylor-tail theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation`.

Proof boundary: this artifact computes the degree-16 coefficient
and stress-tests the current finite baselines one Taylor step further.
It does not prove an infinite Taylor-tail estimate.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: 7 rows, 0 issues, 4 tested continuation rows, 3 stencil-sign-preserving rows, 1 stencil-sign-failure rows, 0 ready-to-apply rows
```

## Degree-16 Coefficient

```text
degree: 16
sign: negative
ratio c16/c0: [-237753170.9494156506565064820654116030884 +/- 2.39e-32]
tail radius: 1583608363045672366325201511766912819440207476015802478781126089830842341935593456e-8958
```

## Continuation Summary

```text
tested continuation rows: 4
pointwise budget failures: 4
stencil-sign-preserving rows: 3
stencil-sign-failure rows: 1
half-safety stencil failures: 3
degree-14 baseline survivors: 1
degree-14 baseline failures: 1
worst pointwise over-budget factor: 3.010798908654295615E+3
worst stencil abs over half margin: 3.998025843926772743E+0
```

Continuation rows:

```text
nlrgts_M5_T2000: M->6, T=2000, pointwise over-budget=2.028701787378877343E+3, stencil signs preserved=True, failing stencils=[]
  B: increment 1.115522375174383601E-6, abs/half-margin 5.180086598431442089E-2, margin after increment 4.418516431946812336E-5
  companion: increment -2.257207596764649147E-7, abs/half-margin 1.812650007602512478E+0, margin after increment 2.332981128843164882E-8
  weighted_gap: increment 1.283992045514261819E-5, abs/half-margin 3.450019208573621088E-1, margin after increment 8.727382750837995921E-5
nlrgts_M6_T2000: M->7, T=2000, pointwise over-budget=3.010798908654295615E+3, stencil signs preserved=True, failing stencils=[]
  B: increment -1.906159608985152175E-7, abs/half-margin 8.628052597940853252E-3, margin after increment 4.399454835856960814E-5
  companion: increment 4.663659423254214610E-8, abs/half-margin 3.998025843926772743E+0, margin after increment 6.996640552097379494E-8
  weighted_gap: increment -2.573151850726511302E-6, abs/half-margin 5.896731985266578194E-2, margin after increment 8.470067565765344792E-5
nlrgts_M7_T1000: M->8, T=1000, pointwise over-budget=1.139255507841004152E+3, stencil signs preserved=False, failing stencils=['companion']
  B: increment 1.287046862588277294E-5, abs/half-margin 1.666505251701368401E-1, margin after increment 1.673310514862494618E-4
  companion: increment -4.145075323704713633E-6, abs/half-margin 2.401524460650696046E+0, margin after increment -6.930385931843005219E-7
  weighted_gap: increment 2.205594774658870866E-4, abs/half-margin 3.007449350603783080E+0, margin after increment 3.672349168521610482E-4
nlrgts_M7_T2000: M->8, T=2000, pointwise over-budget=1.248008080287443688E+2, stencil signs preserved=True, failing stencils=[]
  B: increment 2.803860685989834727E-8, abs/half-margin 1.274640059098902619E-3, margin after increment 4.402258696542950649E-5
  companion: increment -8.054205347980491209E-9, abs/half-margin 2.302306453506772200E-1, margin after increment 6.191220017299330373E-8
  weighted_gap: increment 4.346248650748797813E-7, abs/half-margin 1.026260680213612010E-2, margin after increment 8.513530052272832770E-5
```

Interpretation:

The degree-16 continuation removes the previous missing-coefficient
frontier. The `T=1000` degree-14 baseline fails through the companion
stencil, while the `T=2000` degree-14 baseline survives. This is a
finite signal for an explicit large-`T` or `q/T` collar in any direct
signed stencil-tail theorem.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md
outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md
outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
```

Summary:

Adding the degree-16 Taylor coefficient removes the previous next-increment frontier for the degree-14 positive baselines. All four current positive baselines still fail the crude pointwise tail budget, while three preserve the structured stencil signs. The degree-14 T=1000 baseline fails through the companion stencil, and the degree-14 T=2000 baseline survives. This sharpens the direct signed stencil-tail route toward an explicit large-T/q-over-T collar.
