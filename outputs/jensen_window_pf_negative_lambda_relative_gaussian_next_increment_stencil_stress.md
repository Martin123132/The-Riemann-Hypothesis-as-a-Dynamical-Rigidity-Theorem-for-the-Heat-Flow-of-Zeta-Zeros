# Jensen-Window PF Negative-Lambda Relative-Gaussian Next-Increment Stencil Stress

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof
of a uniform Taylor-tail theorem, scaled-curvature monotonicity,
cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress`.

Proof boundary: this artifact compares the known next Taylor
increment with the current pointwise and structured stencil budgets.
It does not prove an infinite Taylor-tail estimate.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: 8 rows, 0 issues, 2 tested next-increment rows, 2 pointwise budget failures, 2 stencil-sign-preserving rows, 0 ready-to-apply rows
```

## Stress Summary

```text
tested next-increment rows: 2
missing next-increment rows: 2
pointwise budget failures: 2
stencil-sign-preserving rows: 2
half-safety stencil failures: 2
worst pointwise over-budget factor: 3.010798908654295615E+3
worst stencil abs over half margin: 3.998025843926772743E+0
```

Tested rows:

```text
nlrgts_M5_T2000: M->6, max |theta| / allowed rho = 2.028701787378877343E+3, stencil signs preserved = True
  B: increment 1.115522375174383601E-6, abs/half-margin 5.180086598431442089E-2, margin after increment 4.418516431946812336E-5
  companion: increment -2.257207596764649147E-7, abs/half-margin 1.812650007602512478E+0, margin after increment 2.332981128843164882E-8
  weighted_gap: increment 1.283992045514261819E-5, abs/half-margin 3.450019208573621088E-1, margin after increment 8.727382750837995921E-5
nlrgts_M6_T2000: M->7, max |theta| / allowed rho = 3.010798908654295615E+3, stencil signs preserved = True
  B: increment -1.906159608985152175E-7, abs/half-margin 8.628052597940853252E-3, margin after increment 4.399454835856960814E-5
  companion: increment 4.663659423254214610E-8, abs/half-margin 3.998025843926772743E+0, margin after increment 6.996640552097379494E-8
  weighted_gap: increment -2.573151850726511302E-6, abs/half-margin 5.896731985266578194E-2, margin after increment 8.470067565765344792E-5
```

Rows awaiting a degree-16 coefficient:

```text
nlrgts_M7_T1000: M=7, T=1000, reason=degree-16 Taylor ratio is not available in the current degree<=14 coefficient scaffold
nlrgts_M7_T2000: M=7, T=2000, reason=degree-16 Taylor ratio is not available in the current degree<=14 coefficient scaffold
```

Interpretation:

The crude pointwise triangle envelope is too conservative for the
known next increment at the current positive baselines. The structured
stencil increments still preserve the tested finite signs, so the live
route is a direct signed stencil-tail theorem rather than a purely
pointwise epsilon_j theorem at this truncation order.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md
outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
```

Summary:

The known next Taylor increment defeats the crude pointwise half-safety envelope on both tested positive baselines, with worst over-budget factor 3.010798908654295615E+3. However, the structured B, companion, and weighted-gap next-increment stencils preserve all tested finite signs, while the companion stencil exceeds the half-safety margin. This points away from a crude pointwise triangle proof and toward a direct signed stencil-tail theorem.
