# Order-Eleven Sparse-H23 Lower-Bridge Progress

Date: 2026-07-17

Status: rigorous full-source and canonical segment-prefix progress toward an
open continuum target. This is not a proof of the full lower bridge, order eleven,
PF-infinity, `Lambda<=0`, RH, or a Clay-prize result.

## Target

The remaining order-eleven endpoint premise is

```text
y_1''(t) <= 6000/t^2 for every real t >= 1252.
```

The localized lower bridge is `1252 <= t <= 5700`. It contains 8,896
half-cells, 17,792 quarter blocks, and 279 resumable root segments. None of
that full interval coverage is claimed here.

## Repaired Route

The direct interval recurrence lost about 26 decimal orders to dependency
wrapping. A finite H0-H8 stencil repaired low blocks but lost point
cancellation at larger `t`. The active route instead uses:

1. Exact 896-bit Arb H0-H23 jets on the even lattice `1244,1246,...,5708`.
2. The localized H24 wall to propagate normalized H0-H16 coefficients by at
   most one unit to every required half-grid target.
3. Common-variable Taylor models of degrees `(16,15,14)` for `H,H',H''` on
   each quarter block.
4. Coefficient-wise shifted recurrences and explicit product, stable-log,
   analytic-tail, and input-remainder bounds through the eighth nested stage.

The canonical driver binds the full source manifest, H caches, generator,
propagation core, Taylor-model core, and its own source by SHA-256. Partial
sources are accepted only by named pilot mode and cannot populate the
canonical segment cache.

## Source Certificate

```text
required exact anchors: 2233
validated anchors:      2233
anchor range:           1244..5708, step 2
cache SHA-256:           a6406ef805825962aa3dba7ea883a6fc8f741dcf8cd36ed50629a271da1ae90c
manifest SHA-256:        1df0cd867f6f7dc6be669795404081fb251eba732d9cd34c5d0f036926cfa8f1
row contract:            order11-lower-sparse-point-h0-h23-step2-p896-b180-n80-w15-t30-v1
```

Every row passed independent identity, interval-parse, mode-bracket,
target-ball, manifest-hash, and propagation-geometry checks. Rows 0, 1116,
and 2232 were rebuilt from quadrature and matched the cache exactly. The
source-independent gates also pass five polynomial plus five H24 propagation
fixtures, nine forced-truncation products, and 81 stable-log enclosures.

## Canonical Segment Checkpoint

```text
validated segments:     32/279
covered interval:       1252..1756
validated half-cells:   1008/8896
validated quarter blocks: 2016/17792
maximum scaled upper:   38.4152805598063436179290619863693617309716426022555744061374
minimum margin:         5961.58471944019365638207093801363063826902835739774442559385
segment cache SHA-256:  57db8dded767aeef5f79f1d8c0cf997b1b76f6d7a56deebbe56f4f5df07adf78
run contract SHA-256:   b3e2f8fb4edea88815601e85b51ae6ed3c0cc95b7eb3336d9dfdcf6d38d729a4
```

The streaming independent checker verifies the immutable run contract, exact
quarter-block adjacency, expansion anchors, model degrees, all eight stable
stages, positive margins, and segment extrema. This is a resumable prefix,
not coverage of `1756<t<=5700`.

## Cell Pilots

Each row below certifies both quarter blocks of one half-cell over hashed
inputs. The numbers are outward upper bounds for `t^2 y_1''(t)`.

| Half-cell | Quarter 1 | Quarter 2 | Source |
|---|---:|---:|---|
| `1252..1252.5` | 36.9763558314 | 36.9773107178 | canonical prefix |
| `1500..1500.5` | 37.7880503817 | 37.7887524981 | canonical prefix |
| `2200..2200.5` | 39.1998972625 | 39.2002629161 | canonical prefix |
| `3000..3000.5` | 40.0925956822 | 40.0928126139 | canonical full source |

All eight bounds are below 6000 and pass the independent pilot-artifact
checker over canonical source records.

## Open Work

1. Complete the remaining 247 canonical segments and independently check all 17,792 quarter
   blocks, extrema, continuity of coverage, and source hashes.
2. Complete the finite-saddle/asymptotic splice beyond `t=5700`.
3. Only then promote the continuum curvature premise and compose the existing
   order-eleven endpoint and heat-flow reductions.

Bounded segment resume command:

```powershell
python -u work/rh_compute/scripts/jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments.py --workers 4 --runtime-seconds 9000
```
