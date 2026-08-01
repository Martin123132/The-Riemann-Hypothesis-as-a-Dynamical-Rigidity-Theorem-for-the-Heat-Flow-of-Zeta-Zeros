# Order-Eleven Sparse-H23 Lower-Bridge Progress

Date: 2026-07-22

Status: rigorous complete finite lower-bridge theorem, now composed with the
compact and saddle certificates into the global first-summand theorem and the
fixed-order-eleven lambda-zero result. This is not PF-infinity, `Lambda<=0`,
RH, or a Clay-prize result.
It is not a proof of RH or `Lambda <= 0`.

## Target

The order-eleven endpoint premise discharged by the completed chain is

```text
y_1''(t) <= 6000/t^2 for every real t >= 1252.
```

The localized lower bridge is `1252 <= t <= 5700`. It contains 8,896
half-cells, 17,792 quarter blocks, and 279 root segments. The complete interval
coverage is certified below.

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

## Canonical Segment Certificate

```text
validated segments:     279/279
covered interval:       1252..5700
validated half-cells:   8896/8896
validated quarter blocks: 17792/17792
maximum scaled upper:   41.4308043260243735131882308377712013897305561203006129174507
minimum margin:         5958.56919567397562648681176916222879861026944387969938708254
segment cache SHA-256:  7b1b5fc81dff51f01debc57e04fdf88e6e421ed7d107f49063f829172fd2792e
run contract SHA-256:   b3e2f8fb4edea88815601e85b51ae6ed3c0cc95b7eb3336d9dfdcf6d38d729a4
```

The streaming independent checker verifies the immutable run contract, exact
quarter-block adjacency, expansion anchors, model degrees, all eight stable
stages, positive margins, and segment extrema. The `--require-complete` gate
passes over every segment and quarter block.

## Cell Pilots

Each row below certifies both quarter blocks of one half-cell over hashed
inputs. The numbers are outward upper bounds for `t^2 y_1''(t)`.

| Half-cell | Quarter 1 | Quarter 2 | Source |
|---|---:|---:|---|
| `1252..1252.5` | 36.9763558314 | 36.9773107178 | canonical certificate |
| `1500..1500.5` | 37.7880503817 | 37.7887524981 | canonical certificate |
| `2200..2200.5` | 39.1998972625 | 39.2002629161 | canonical certificate |
| `3000..3000.5` | 40.0925956822 | 40.0928126139 | canonical certificate |

All eight bounds are below 6000 and pass the independent pilot-artifact
checker over canonical source records.

## Completed Handoff

The finite saddle range `2001/1000<=u<=20` and the asymptotic range `u>=20`
are now rigorous and independently checked. Only the compact real-t handoff

```text
5700 <= t <= V'(2001/1000) < 38020
```

was the final interval needed for the global continuum curvature premise. The
current compact pilots support an adaptive step-four exact H0-H23 lattice.
The p896/b200/n96 contract passes representative neighborhoods through
`t=25004`; p896/b200/n108 passes at `t=30004`, `35004`, and `37800` with
largest tested upper `1949.4594`. Both endpoint collars through `t=5700` and
`t=38020` pass, including the two boundary H24 tiles. The active full compact
source cache validates all `8085/8085` exact anchors through `t=38028`; its
final SHA-256 is
`63dead138ec03af796ffd37a53307ac21895c31feec960b819df7d3013e2d5bb`.
The profile split is exactly 4,832 lower plus 3,253 upper rows. Independent
exact rebuilds match on both sides of the transition,
at `t=25016` and `t=25020`. The checker also binds the profile counts, next
task, and explicit inner-kernel versus outer-source contract provenance.
The canonical compact driver and checker validate all `2020/2020` contiguous
segments and 129,280 quarter blocks over `5700..38020`. Canonical transition
segment 1207 exactly matches its independent pilot record, and all later rows
use the upper p896/b200/n108 profile. The largest scaled upper in the validated
cache is `2122.8651669101`, leaving margin `3877.1348330899`. The final cache
SHA-256 is `70062a5bf5c906ef894a559db1cc8ebab3379d69ffcaf7f6f88561986b7ae9a2`.
The compact promoter and global four-source composer pass, proving
`y_1''(t)<=6000/t^2` for every real `t>=1252`. The endpoint and heat-flow
composition then proves `Q_(11,n)(0)>0` for every integer `n>=0` and the
arbitrary-column signed-minor signs through fixed order eleven.

The remaining frontier is no longer this lower or compact bridge. It is order
twelve and, ultimately, an all-order mechanism. No fixed-order certificate is
being promoted to PF-infinity, `Lambda<=0`, or RH.

The implemented geometry, pilot evidence, corruption gates, and resume commands
are specified in
`outputs/jensen_window_pf_compound_order11_compact_block_handoff.md`.

Complete-certificate reproduction:

```powershell
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments.py --require-complete
python work/rh_compute/scripts/jensen_window_pf_compound_order11_sparse_h23_lower_bridge_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_sparse_h23_lower_bridge_certificate.py
```
