# Jensen-Window PF Negative-Lambda Half-Width Tail Target

Date: 2026-07-06

Status: finite-rejected target. This is not a proof of a replacement
scaled-defect tail theorem, zeta cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_half_width_tail_target`.

Proof boundary: this artifact records that the fixed half-width
scaled-defect route is rejected by finite stress. It does not prove
a corrected route, does not prove the separate monotone defect bridge,
and does not establish `Lambda <= 0`.

Machine-readable target:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_half_width_tail_target.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_half_width_tail_target.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py
```

Current result:

```text
validated Jensen-window PF negative-lambda half-width tail target: 9 rows, 0 issues, 0 ready-to-apply rows, 0 live routes, 430 half-width rows, 17 half-width failures
```

## Target Statement

Let:

```text
d_k = 1 - x_k
s_k = ((2*k+1)/2) * d_k
```

The fixed half-width tail target was:

```text
0 <= s_k <= 1/2 for all k >= 150
d_(k+1) <= d_k for all k >= 149 remains separately required
```

Finite k150 stress rejects this fixed half-width target; it is retained only as a documented failed shortcut.

The finite prefix currently supplies:

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_150
finite contractions: x_1..x_149
next finite coefficient needed: A_151
```

## Finite Evidence

```text
half-width rows: 430 / 447
half-width failure rows: 17
one-third failure rows: 268
minimum half-width margin: -1.107359370078869990E-2 at lambda=-25.0, k=149
maximum scaled defect: 5.110735937007886999E-1 at lambda=-25.0, k=149
minimum one-third margin: -1.777402603674553666E-1 at lambda=-25.0, k=149
```

First half-width failures inferred from the increasing scaled-defect frontier:

```text
lambda=-25.0: first failure k=133, failure rows=17, max s=5.110735937007886999E-1 at k=149
```

## Live Routes

```text
none for the fixed half-width threshold; repair routes must target the exact cone or an adaptive scaled-defect buffer
```

Rejected shortcuts:

```text
one-third-width buffer
fixed half-width buffer
scaled-defect nonincrease
endpoint PF/RH/Newman assumptions
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
```

Summary:

The k150 scaled-defect frontier finite-rejects the fixed half-width bound s_k<=1/2: only 430/447 checked rows pass, with 17 half-width failures and minimum margin -1.107359370078869990E-2 at lambda=-25.0, k=149. The route must be replaced by an exact-cone or adaptive scaled-defect tail target plus the separate monotone defect bridge.
