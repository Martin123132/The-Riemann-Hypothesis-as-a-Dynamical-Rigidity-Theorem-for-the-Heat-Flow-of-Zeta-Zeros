# Complete displayed lower-fold height window

Date: 2026-08-11

Status: rigorous canonical lower-fold enclosure on all displayed source
heights; not a proof of the complete source formula or a selector-fold join

The saved evaluator output requests the exact 15-point roster

```text
t-0.07, t-0.06, ..., t, ..., t+0.06, t+0.07,
t=10^10.                                                (DW1)
```

The roster is parsed from the pinned source output.  Its printed Hardy values
are not used as proof data.

Across the same height cell, the exact inverse roots remain in

```text
39852 < N_-(C;t) < 39853,
39936 < N_+(C;t) < 39937,
39894 < sqrt(t/(2pi)) < 39895,                         (DW2)
```

with minimum integer-boundary margin
`[0.1084806914476379886815640987129354315759906196584169566915471929534286 +/- 8.75e-64]`.  Hence the source turning
roster stays exactly `39853..39936` with 84 modes throughout the window.

Throughout the one selector-stable interval `|t-10^10|<=0.07`, the exact
fixed-selector identity remains

```text
G_tt=(-lambda G-i G')/beta^2.                          (DW3)
```

An Arb cover with 2560 intervals proves
`|Ai(x)|<0.54` on the complete Airy argument range.  A differently partitioned
independent checker repeats the cover.  Therefore

```text
|S''(t)|<0.146,                                         (DW4)
|S(t)-S(t0)-(t-t0)S'(t0)|<0.000358,                   (DW5)
|S(t)-S(t0)|<0.009494.                                 (DW6)
```

All 15 displayed heights lie in this single certified cell.
This proves height variation only for the canonical grouped lower-fold
transform `S`; it does not certify the other source blocks or their final
assembly.

Proof boundary: complete displayed-height window for the canonical lower-fold
transform only.  No adjacent-selector crossing, complete source evaluator
error theorem, lower-interior join, `T_upper`, `Lambda<=0`, RH, or prize-level
conclusion is proved.
