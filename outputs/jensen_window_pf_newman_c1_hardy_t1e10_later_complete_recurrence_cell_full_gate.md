# Hardy t=1e10 complete later-cell recurrence-tail transport

Date: 2026-08-09

Status: complete later-cell tail and source multiply-add rounding columns; not a proof of source q, the outer Hardy representation, or RH

The accepted nonzero coefficient cells for all `1040`
recursive calls in blocks 21--28 preserve the exact factored cubic Legendre
remainder at 70 and 110 decimal digits.  Every saved point tail is contained
in its cell enclosure.  The common transformed-child phase is removed at unit
modulus before absolute values.

The pinned source environment uses radix-2, 113-bit, round-to-nearest
ties-to-even arithmetic.  Parent, reconstructed-child, anchored-multiplier,
child-state, and `tpm=-2*pi` rounding preimages are contained for every call,
and every emitted source complex multiply-add lies inside its operation bound.

```text
later cells                                      1040
precision overlaps                               1040 / 1040
maximum local complete tail                      1.63664601327194522170888907502650825178431619309559E-3
maximum cell/point tail inflation                1.00000001330007527320867365347722463515343043039605E+0
maximum transported complete tail                7.36544950300080392703004465341499402300333384753994E-4
maximum transported source multiply-add rounding 1.82398763256928169651187850850036596334521037598653E-32
maximum partial complete output                  1.18534368330532488984495230980059347380722157622465E-3
minimum remaining 0.005 margin                   3.81465631669467511015504769019940652619277842377535E-3
outputs below 0.005                              15 / 15
```

No cancellation between recursive calls is used.

## Boundary

This closes the analytic complete-tail and source recurrence multiply-add
rounding columns on the finite saved later cells.  It does not yet enclose
source `q` or special-function variation over those cells.  The outer Hardy
representation/remainder, height uniformity, `Lambda<=0`, PF-infinity, RH,
and a prize-level conclusion remain open.
