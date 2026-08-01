# Order-Twelve Sparse-H23 Physical Lower Bridge

Date: 2026-08-01

Status: **incomplete interval certificate; a contiguous 125-segment lower
prefix is validated, while the remaining 138 segments are open**.
This is not a proof of order twelve, RH, or `Lambda <= 0`.

## Target

The remaining physical lower obligation is

```text
v_1''(t) <= 8000/t^2 for every real 1503 <= t <= 5700.
```

The deterministic cover has `263` root segments, `8,394` half-cells, and
`16,788` quarter blocks of width `1/4`.

## Exact Jet Transport

For `0<=d<=16`, an exact `H0`--`H23` jet at `a` is transported to `x` by

```text
H^(d)(x)/d! in
  sum_(m=d)^23 H^(m)(a)/m! binom(m,d)(x-a)^(m-d)
  + [-E_d,E_d],

E_d = M_24 |x-a|^(24-d) / ((24-d)! d!),
|x-a| <= 2.
```

The dedicated checker passes seven exact-polynomial translations, seven
exponential H24 remainder enclosures, and a strict rejection at distance
`201/100`.

## Bound Sources

```text
inherited exact H0-H23 anchors:
  rows=2233, t=1244..5708, step=2
  SHA-256=a6406ef805825962aa3dba7ea883a6fc8f741dcf8cd36ed50629a271da1ae90c

broad H2-H24 wall:
  rows=42180, t=1491..5709, step=1/10
  SHA-256=0c2b123ca508280a0d8ba5a325a020da80ffc796a1d0c41fb87913b0f59cf5fd

compact H2-H24 wall:
  rows=32336, t=5692..38028, step=1
  SHA-256=2040c13cdb9610602c25099962a69e7588662a10b0d7dbc627531dab853d35e2
```

At the extreme stencil target `t=5710`, the nearest exact anchor is `5708`.
The broad H24 wall covers through `5709`, and the compact unit tile covers
`5709..5710`; the preflight finds no gap in the mixed collar.

The independent preflight validates all `263` segment specifications,
`16,788` quarter-block specifications, and `8,435` distinct stencil targets.
The separate algebra checker also validates three recurrence products and
`900` `D0`--`D9` value/first/second-derivative interval enclosures.

## Validated Prefix

Named segment pilots `0`, `131`, and `262` pass all `148` quarter blocks. Their
largest scaled upper is `67.7240880954849832`, including the final mixed-wall
collar.

The canonical cache validates segments `0..124`, all `8,000` quarter blocks,
and therefore every real point in

```text
1503 <= t <= 3503.
```

Its largest scaled upper is `49.5710796974064546`; the smallest margin lower
is `7950.42892030259355`. The next deterministic segment is index `125` on
`3503..3519`.

## Promotion Gates

The segment driver is append-only and deterministic. The promotion generator
is fail-closed unless all `263` segment rows and all `16,788` quarter blocks
pass the independent checker with positive margin.

```powershell
python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h0_h23_propagation.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py --preflight-only

python work/rh_compute/scripts/jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py `
  --pilot-segment 0 --pilot-segment 131 --pilot-segment 262 --workers 2

python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py `
  --pilot work/rh_compute/results/jensen_window_pf_compound_order12_sparse_h23_lower_bridge_pilot.json

python work/rh_compute/scripts/jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py `
  --workers 1 --runtime-seconds 7200

python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py `
  --require-complete
```

The repository publishes the validated prefix as a lossless gzip archive. Its
archive wrapper restores a temporary JSONL stream and runs the same independent
segment checker:

```powershell
python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_archive.py
```

Only after the complete checker passes may the promotion generator and its
checker be run.

```powershell
python work/rh_compute/scripts/jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate.py
```

## Boundary

The source manifests, exact transport lemma, deterministic geometry, and
mixed-wall coverage are ready, and the first `125` canonical lower segments are
proved. Segments `125..262`, the physical compact interval, global first-summand
composition, fixed order twelve, PF-infinity, RH, and `Lambda <= 0` remain
separate and open.
