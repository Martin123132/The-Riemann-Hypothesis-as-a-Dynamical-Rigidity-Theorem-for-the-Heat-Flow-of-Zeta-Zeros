# Hardy block-20 endpoint-linearization decomposition full gate

Date: 2026-08-07

Status: rigorous all-374 finite-input decomposition; not a proof of a height-uniform remainder theorem or RH

The endpoint-linearization pilot is replayed at 70 and 110 decimal digits on
all 374 saved recursive calls.  Each chain reuses the certified full-ray
cutoff and is fsynced atomically to an append-only cache.

On every call the exact identities

```text
endpoint_half+a_endpoint+B_linear+C_linear = paper_half+W5,
Q_exact-Q_paper = R_b+R_c+zero_mode-W2-(W3 or W4)
```

enclose zero.  Aggregate bounds are

```text
maximum generic/W5 identity gap       <= 1.82388016188312862703883057071963467016256349900634E-106
maximum correction decomposition gap  <= 5.23698368491970890712533526212268952804151922464371E-14
maximum component triangle sum        <= 2.12607629404946580399332176464137723960672396046398E+2
maximum triangle/correction ratio      <= 1.38940364315770530666447663559917886409495020516720E+5
```

The worst triangle ratio is witnessed by chain
`36`.  It
measures cancellation already lost by a fully componentwise absolute bound;
it is a diagnostic for designing the next analytic majorant, not a fitted
constant and not a recurrence budget.

This gate rigorously localizes the complete finite correction into generic
endpoint-linearization remainders and the exceptional W2--W4 replacements.
It supplies neither uniform-height piece bounds nor recurrence propagation or
outer Hardy control.  It is not a proof of `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
