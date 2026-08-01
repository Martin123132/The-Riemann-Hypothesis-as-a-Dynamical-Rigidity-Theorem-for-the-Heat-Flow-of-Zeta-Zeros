# Order-Eleven Compact Block Handoff

Date: 2026-07-22

Status: complete and independently validated compact quarter-block theorem,
now composed through the global first-summand, endpoint, heat, and fixed-order
lambda-zero steps. This is a fixed-order-eleven result, not PF-infinity or RH.
It is not a proof of RH or `Lambda <= 0`.

## Current Readiness

The adaptive point source is complete and independently validated at
`8085/8085` rows on `5692..38028`. Its final SHA-256 is
`63dead138ec03af796ffd37a53307ac21895c31feec960b819df7d3013e2d5bb`.
That final hash is frozen into the canonical compact run contract.

The driver and independent checker are implemented. Named pilots at segments
0, 1207, and 2019 validate 192 quarter blocks, including both true endpoints
and the mixed-profile transition segment. Their largest scaled upper is
`2122.86516691011593226179209689464539687227984653146836129212`, and
their smallest margin is
`3877.13483308988406773820790310535460312772015346853163870788`.

The canonical append-only cache validates all `2020/2020` contiguous segments
and 129,280 quarter blocks from `t=5700` through `t=38020`. Its final SHA-256 is
`70062a5bf5c906ef894a559db1cc8ebab3379d69ffcaf7f6f88561986b7ae9a2`
and its byte length is `2361570130`. Its largest
scaled upper is
`2122.86516691011593226179209689464539687227984653146836129212`, and
its smallest margin is
`3877.13483308988406773820790310535460312772015346853163870788`.
Both extrema occur in final segment 2019, block 63, at anchor `152079/4`, on
the quarter block ending at `t=38020`. Canonical transition segment 1207 is
identical to the independent pilot record, with structured-record SHA-256
`fce51d9b8ba0d87554dda72d095dc61ebec4286a0c35e303bc8ee83393c9e0bd`.
It records both source profiles and has largest scaled upper
`789.049774390264452939526925992731550707469231213217022028614`.
The first wholly upper-profile segment, index 1208, has largest scaled upper
`48.0555004956989427715303537194518087003397446450089055323501`.
Across the fully upper-profile region, the largest scaled upper is the global
value `2122.86516691011593226179209689464539687227984653146836129212`
in final segment 2019, block 63, at anchor `152079/4`. The final bounded
session added 115 segments and independently validated the complete cache.
Daytime runs use at most two below-normal workers; unattended-night mode
adapts between zero and four after a baseline CPU sample and targets at most
75 percent total CPU.

The compact theorem promoter and independent checker pass on all 2,020
segments. The global four-range composer then joins the lower, compact,
finite-saddle, and asymptotic-saddle certificates to prove
`y_1''(t)<=6000/t^2` for every real `t>=1252`. Its circularity counters remain
at zero full-kernel, zero heat-forward, and zero RH claims. The separately
checked endpoint and heat composition proves `Q_(11,n)(0)>0` for every
integer `n>=0` and transfers all completed lower layers to arbitrary columns
through fixed order eleven.

## Exact Geometry

```text
target interval:          5700..38020
root segment width:       16
root segments:            2020
half-cell width:          1/2
half-cells:               64640
quarter-block width:      1/4
quarter blocks:           129280
point-anchor lattice:     5692..38028, step 4
maximum propagation:      2
lower source profile:     5692..25016, p896/b200/n96
upper source profile:     25020..38028, p896/b200/n108
```

The root segmentation closes exactly. Segment 0 is `[5700,5716]`, segment
1207 is `[25012,25028]` and crosses the source-profile boundary, and segment
2019 is `[38004,38020]`.

## Validated Inputs

The adaptive H0-H23 source is complete at 8,085 rows and passes its checker
with `--require-complete`. Its final cache hash is bound by the compact run
contract.

The inherited H2-H24 source is already complete at 32,336 unit tiles on
`5692..38028`; its current cache SHA-256 is
`2040c13cdb9610602c25099962a69e7588662a10b0d7dbc627531dab853d35e2`.
The checked order-eleven boundary extension supplies the two additional tiles
`5691..5692` and `38028..38029`.

The computational cores to reuse unchanged are:

```text
work/rh_compute/scripts/jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core.py
work/rh_compute/scripts/jensen_window_pf_compound_order11_shifted_taylor_model_core.py
```

The complete lower-bridge segment driver supplies the canonical streaming,
run-contract, segment-cache, and checkpoint architecture. The compact
step-four geographic and endpoint pilots supply the source lookup,
propagation-distance-two geometry, quarter-block evaluator, and one-sided
endpoint behavior.

## Canonical Driver Requirements

1. Refuse a partial adaptive source and verify its manifest, cache hash,
   source contract, profile counts, contract provenance, and final row count.
2. Verify the complete H2-H24 unit cache and both boundary-extension tiles.
3. Enumerate all 2,020 width-16 segments deterministically and append one
   complete segment record at a time.
4. Within each segment, enumerate 32 half-cells and exactly 64 adjacent quarter
   blocks, with no overlap or gap.
5. Select the nearest admissible step-four anchor deterministically by
   `(distance, anchor)` and reject propagation distance greater than 2.
6. Preserve each source row's recorded profile. In particular, do not infer
   the profile from a block center across segment 1207.
7. Use one-sided blocks at the true endpoints and the checked H24 boundary
   tiles whenever the stencil reaches `5691` or `38029`.
8. Require every block to return `passed=true`, record its scaled curvature
   upper bound and positive margin, and fail the segment atomically otherwise.
9. Bind the source caches, manifests, boundary artifact, generator, checker,
   propagation core, Taylor-model core, and driver by SHA-256 in an immutable
   run contract.
10. Support bounded runtime, append-only resume, four-worker operation, and
    completion of only the in-flight segment batch after the runtime ceiling.

## Independent Checker Gates

The checker must stream the segment cache and independently verify the run
contract hashes, segment order, exact domain endpoints, block count, quarter
adjacency, expansion anchors, source profiles, propagation distances, model
degrees, stable-stage records, positive margins, and segment extrema. A
`--require-complete` mode must demand all 2,020 segments and all 129,280 blocks.

Before a full run, exact pilot mode must pass segment 0, transition segment
1207, and final segment 2019. Deliberate corruptions of segment order, profile
selection, source hash, propagation distance, and quarter adjacency must each
be rejected.

All three named pilots pass, and all five listed corruptions are rejected by
the independent checker.

## Reproduction

```powershell
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py --require-complete --rebuild-index 4831 --rebuild-index 4832
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_adaptive_h23_segments.py --workers 3 --pilot-segment 0 --pilot-segment 1207 --pilot-segment 2019
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_h23_segments.py --pilot work/rh_compute/results/jensen_window_pf_compound_order11_compact_adaptive_h23_segment_pilot.json
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_adaptive_h23_segments.py --workers 2 --runtime-seconds 4800
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_h23_segments.py --require-complete
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.py
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_adaptive_h23_certificate.py
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_first_summand_curvature_certificate.py
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_first_summand_curvature_certificate.py
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_m100_entry_certificate.py
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_m100_entry_certificate.py
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_lambda0_completion_certificate.py
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_lambda0_completion_certificate.py
python -u work/rh_compute/scripts/check_jensen_window_pf_compound_order11_completion_composition_gates.py
```

## Proof Boundary

The block cache proves only the compact first-summand premise. The checked
downstream artifacts now complete its composition through fixed order eleven:
`Q_(m,n)(0)>0` for `1<=m<=11`, every `n>=0`, and the corresponding signed
consecutive-row arbitrary-column minors through order eleven. No order above
eleven is claimed. The all-order Jensen/PF bridge, PF-infinity, `Lambda<=0`,
RH, and the Clay-prize conclusion remain open.
