# Extreme-detuning full-line pilot

Date: 2026-08-12

Status: rigorous two-endpoint diagnostic on the common event contour; not a
proof of the 399-event propagation theorem

The event-parameter geometry gate proves that `R=14`, `Im z=1`, and
`|Re z|<=21` contain all exact and canonical saddles and retain positive
lifted-contour margins throughout the 399-event range.  This pilot now runs
the signed retained-Fresnel/logistic minus Airy integral at the two extreme
event detunings, before absolute values:

```text
event 398 (mode 39695, 4m-C=-797): normalized=[0.012353512577189333317854926815428442954392096453375415876507759094238281250000000 +/- 3.01e-14], physical=[0.0055610344675477890067053845932687136666139338103675981983542442321777343750000000 +/- 1.36e-14]
event 397 (mode 40093, 4m-C=795): normalized=[0.012798667260632164341010487147026752907663649239111691713333129882812500000000000 +/- 4.67e-13], physical=[0.0057614244799064551697320778027300303136826187255792319774627685546875000000000000 +/- 2.10e-13]
```

Prototype target survives both extremes: `False`.
Prototype target is rigorously excluded at both extremes:
`True`.

This diagnostic tests quantitative propagation after contour geometry.
It does not prove a uniform bound between the endpoints, a
paired-event cancellation theorem, all 399 event remainders, complete
`T_upper`, `Lambda<=0`, RH, or a prize-level conclusion.
