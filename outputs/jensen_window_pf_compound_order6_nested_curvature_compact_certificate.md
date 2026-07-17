# Jensen-Window PF Compound Order-Six Nested Curvature Compact Certificate

Date: 2026-07-13

Status: rigorous first-summand order-six curvature theorem on
`321<=t<=V'(2)`. This is not a proof of the remaining saddle ray,
order-six entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py
```

## Extended Derivative Cache

The existing 107452-tile compact cache supplies `H^(2),...,H^(8)`.
An aligned extension supplies only the two genuinely new columns
`H^(9),H^(10)`. Both files are hashed in the JSON artifact.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_h_tiles.jsonl
work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_compact_h9_h10_tiles.jsonl
work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_compact_right_collar.json
aligned rows=107452
```

## Certified Cover

Every adaptive block uses a common `t+-4` derivative collar and
outward-rounded stable power-series arithmetic. It proves

```text
J_1(t)>0, R_1(t)>0, S_1(t)>0,
p_1''(t)<=200/t^2 for 321<=t<=V'(2).
```

Accepted blocks: `38`.
Largest scaled curvature upper: `1.97869982971746450359369590541639821427159936246516490973153E+2`.
Weakest curvature margin lower: `1.11892969209662504790677781810430993232578460578084967353093E-7`.
Weakest third-gap lower: `8.19820850331485930668672604808464431251440082120441608273389E-5`.

## Remaining Ray

The unresolved continuous region begins at saddle mode `u=2` and is
handled separately by finite-ray and asymptotic-ray certificates.

```text
outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md
outputs/formal_core.md
```
