# Hardy t=1e10 complete later-recursive signed coefficient-cell bridge

Date: 2026-08-09

Status: complete finite corrected-model cell column at one saved height; not a proof of the complete recurrence, height uniformity, or RH

The two-witness pilot is scaled to all `1040` recursive
calls in blocks 21--28.  Every row has a nonzero-radius selector-stable box,
overlapping 70/110-digit Arb enclosures of the signed analytic correction
`Qexact-Qpaper`, and an enclosure containing the saved exact-point correction.

The calculation keeps

```text
Dcorr = Ppi - I69pi - Qpaper = Qexact - Qpaper
```

as one complex interval before taking its magnitude.  It uses no signed
cancellation between different calls.

Cells closed: `1040/1040`

Maximum cell/point inflation: `1.03235836575090763952065052440786589634990265841376E+0`

Worst all-recursive cell complete output: `4.48798733005244497141947844440854195181195374505537E-4`

Minimum margin below `0.005`: `4.55120126699475550285805215555914580481880462549446E-3`

Boundary: this closes the analytic corrected-model coefficient-cell column at
the finite saved height only.  Source binary128 rounding over boxes, the
full result is not a complete recurrence.  The
complete cubic recurrence tail in blocks 21--28, the outer Hardy
representation/remainder, and height-uniform selector and accumulation
theorems remain open.  It does not prove Lambda<=0, PF-infinity, RH, or a
prize-level conclusion.
