# Hardy t=1e10 later-block signed coefficient-cell bridge pilot

Date: 2026-08-09

Status: rigorous two-witness corrected-model coefficient-cell pilot; not a proof of the full recurrence, height uniformity, or RH

The worst exact/retained-ratio witnesses from block 21 and the later recursive
blocks are promoted from exact coefficient points to nonzero-radius,
selector-stable boxes.  The source bridge is kept signed algebraically:

```text
Dsrc = Psrc - MC - qsrc
D1   = Dsrc - (I69src - MC - t5src)
D2   = D1 + qextra = Psrc - I69src - Qsrc,drop
Dcorr = D2 + (Ppi-Psrc) - (I69pi-I69src) - (Qpaper-Qsrc,drop)
      = Ppi - I69pi - Qpaper = Qexact - Qpaper.
```

The last expression is evaluated as one complex Arb interval on each cell;
the exact contour and paper terms are not separately converted to absolute
majorants before subtraction.  Both 70- and 110-digit enclosures overlap the
saved exact-point correction.

Witness cells closed: `2/2`

Maximum cell/point correction inflation: `1.02786775833167184220076791455195161607238751829546E+0`

Maximum transported two-cell inflation: `4.12048229426223632554000851598731380612357170791402E-8`

Remaining finite exact-point margin: `4.55601112596361510673046545569964374799791938123474E-3`

Boundary: these are analytic corrected-model cells.  The theorem does not yet
transport source binary128 rounding, the complete cubic recurrence tail, all
1,040 later recursive calls, the outer Hardy representation/remainder, or any
height-uniform conclusion.  Those remain explicit obligations.

This is not a full recurrence, height-uniform theorem, or RH proof.
