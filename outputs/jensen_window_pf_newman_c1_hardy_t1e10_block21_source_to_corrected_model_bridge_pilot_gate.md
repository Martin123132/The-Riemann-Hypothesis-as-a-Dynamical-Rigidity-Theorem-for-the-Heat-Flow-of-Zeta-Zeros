# Hardy t=1e10 block-21 signed bridge pilot

Date: 2026-08-09

Status: representative finite signed-bridge gate closed; not a proof of a coefficient neighborhood, arbitrary height, or RH

This pilot tests chains `432, 443, 448, 489`.  They are selected before exact contour
evaluation to cover the largest output-15 block-21 contribution and every
occupied `(W3/W4, sign(Phi3))` branch.  All have parent length 119, so none is
a relabelled block-20 length-104 call.

For each witness, exact source algebra is transported in four signed stages:

```text
Dsrc = Psrc - MC - qsrc
D1   = Dsrc - (I69src - MC - t5src)
D2   = D1 + qextra = Psrc - I69src - Qsrc,drop
Dcorr = D2 + (Ppi-Psrc) - (I69pi-I69src) - (Qpaper-Qsrc,drop)
      = Ppi - I69pi - Qpaper.
```

An independent sign-aware endpoint-contour calculation supplies
`Qexact = Ppi-I69pi`.  Thus `Dcorr=Qexact-Qpaper` closes without a fitted
correction or cancellation between calls.

Witnesses closed: `4/4`

Maximum formula-to-independent-target gap: `3.57018059931571603293765093667389010079205036163330E-14`

Maximum signed bridge identity gap: `3.57018060460967195327702805585090573003981262445450E-14`

Exact corrections inside the earlier retained majorant:
`4/4`

The result is a representative finite bridge, not a proof for all 212
block-21 pivots, blocks 22--28, a coefficient neighborhood, arbitrary height,
or RH.  Its next valid use is to promote the length-generic machinery to the
complete block-21 roster with an append-only cache.
