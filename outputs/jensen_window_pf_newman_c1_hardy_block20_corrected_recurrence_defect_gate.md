# Hardy block-20 corrected recurrence-defect gate

Date: 2026-08-06

Status: rigorous finite exact-point recurrence defects after the PSI/ERF replacement; not a proof or uniform recurrence theorem

## Comparison

For each of the 374 `MIT=2` calls, Arb independently sums the exact saved
length-104 parent and length-one-or-two child kernels.  It then compares

```text
D_source = parent - A(source child, source q),
D_corr   = parent - A(source child, source q + delta q_special),
```

where `A` includes the exact saved multiplier, conjugation, and subtract-one
branches.  The correction `delta q_special` is imported from the dyadic Arb
endpoints of Section 11.249, so no decimal re-rounding is used.  Evaluations at
180 and 260 decimal digits overlap for every source and corrected defect.

## Result

```text
maximum independent source/logged defect gap <= 5.26109688231931655659882087643541965280697931813432E-31
minimum source defect magnitude              >= 4.66665048371366388692275634829337527000039622904166E-3
minimum corrected defect magnitude           >= 4.66665048362406825963434178913520943469421998860976E-3
maximum source defect magnitude              <= 2.05015438240018819216760444317311729323888378265636E-1
maximum corrected defect magnitude           <= 2.05015438240021097313750686746075367107459698161638E-1
maximum special-function defect change       <= 7.15521887106601655049865376374079065434789044715391E-4
maximum relative change budget               <= 9.21494747578374093750752321024263389745875540642765E-3
rigorously improved calls                      = 232
rigorously worsened calls                      = 142
interval-indeterminate comparisons             = 0
corrected defects excluding zero               = 374 / 374
```

The special-function replacement is therefore propagated through the actual
one-step recurrence, but it does not remove the finite discrepancy.  The
remaining nonzero defect is an explicit target for the intrinsic real `erfc`,
endpoint-saddle, Euler--Maclaurin, and other mathematical `t5/q` approximation
terms rather than an untracked numerical effect.

This is a finite low-height diagnostic and error certificate.  It does not
derive a sign or uniform bound for the remaining defect, transport it to large
heights, combine the 424 calls with an outer Hardy remainder, prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
