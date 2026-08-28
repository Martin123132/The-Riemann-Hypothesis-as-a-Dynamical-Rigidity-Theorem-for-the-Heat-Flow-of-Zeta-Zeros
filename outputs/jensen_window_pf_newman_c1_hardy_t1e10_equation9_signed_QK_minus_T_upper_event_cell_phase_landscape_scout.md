# Upper-event-cell phase landscape scout

Date: 2026-08-28

Status: cache-backed route diagnostic; not an interval theorem or RH proof

## What the certified local boxes say

The complete derivative has the sign sequence

```text
offset 5: positive, [[0.002059354425517174407213171989373157006318752808162511360 +/- 1.44e-58], [0.002098220163707626655304050118378589135225002808162511360 +/- 1.44e-58]]
offset 10: negative, [[-0.002310835679387249161159873813700601904288334268755159753 +/- 8.13e-59], [-0.002281241963706573347231627315592691748038334268755159753 +/- 8.13e-59]]
offset 15: negative, [[-0.0003586256033402834319658113224915429506103439071384854402 +/- 3.31e-59], [-0.0003203218536926443804203999741486474427978439071384854402 +/- 3.31e-59]]
offset 20: positive, [[0.0004876150181871707903666913722968057518790564860630687427 +/- 9.00e-60], [0.004607758708857545571260017177960868251879056486063068743 +/- 2.92e-58]]
```

Thus a single positive-derivative bridge from offset 20 toward offset 5 is impossible.
Continuity gives at least one stationary point in `(5,10)` and at least one in
`(15,20)`, without proving uniqueness.

All 20 integer-height boxes, all 13 candidate-trough boxes, and all 13
phase-predicted candidate-peak boxes are strictly negative.  The closest sampled
candidate peak to zero is

```text
offset 12.75: upper <= -0.001981574179808376356959342956542968750000 < 0.
```

The observed spacing is approximately `1.5` height units.  This is route-sizing
evidence, not a proved period or a certified list of extrema.

## Next certificate

Build a phase-normalized interval model on the 13 candidate cells covering offsets
`[0,19.5]`, isolate every derivative zero or certify a Taylor/Chebyshev upper
envelope, and handle the remaining `[19.5,19.93]` tail before the existing
offset-20 sign island.

## Proof boundary

Every listed local Arb box is rigorous, and the derivative signs rule out a single positive-derivative bridge.  The sampled candidate peaks are not certified extrema; their negativity does not fill any interval, count all stationary points, or imply RH.
