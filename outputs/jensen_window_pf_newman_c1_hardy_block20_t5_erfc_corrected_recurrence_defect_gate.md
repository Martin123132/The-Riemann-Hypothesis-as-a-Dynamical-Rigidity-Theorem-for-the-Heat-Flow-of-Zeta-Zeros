# Hardy block-20 t5-erfc-corrected recurrence-defect gate

Date: 2026-08-06

Status: rigorous finite recurrence propagation of the admitted PSI/complex-ERF and intrinsic real-erfc replacements; not a proof or saddle remainder theorem

## Comparison

For each of the 374 recursive calls, Arb independently resums the exact saved
length-104 parent and length-one-or-two child.  The three models are

```text
D_source  = parent - A(child, q_source),
D_special = parent - A(child, q_source + delta_q_PSI/complex-ERF),
D_full    = parent - A(child, q_source + delta_q_PSI/complex-ERF
                                  + delta_q_intrinsic-real-erfc).
```

`A` includes the exact saved multiplier, conjugation, and subtract-one branch.
All imported correction balls use their stored dyadic endpoints.  Evaluations
at 180 and 260 decimal digits overlap, and every recomputed `D_special`
overlaps the previously admitted corrected-defect certificate.

## Result

```text
maximum intrinsic-real-erfc defect increment <= 1.97694940650258454908778340217708987471483470054228E-34
maximum increment / prior-defect ratio        <= 6.87419348951299649624507272329385348188163150311951E-33
minimum surviving full-defect magnitude       >= 4.66665048362406825963434180168989386723463716839527E-3
maximum surviving full-defect magnitude       <= 2.05015438240021097313750686731496617932304852838501E-1
minimum defect / max erfc-increment separation >= 2.36053106279529230392032788934410009116999754008535E+31
rigorously improved by real-erfc replacement    = 192
rigorously worsened by real-erfc replacement    = 182
interval-indeterminate magnitude comparisons    = 0
full corrected defects excluding zero           = 374 / 374
```

Thus the source intrinsic real `erfc` accuracy is not the observed recurrence
wall: its largest possible propagated change is over thirty orders of magnitude
below even the smallest surviving defect.  The surviving finite discrepancy is
now quantitatively assigned to the still-unproved endpoint-saddle, finite-`ip`,
Euler--Maclaurin, and remaining `t5/q` approximation contract (plus separately
unbounded arithmetic outside the admitted replacements).

This finite diagnostic does not establish a formula or uniform bound for that
remainder, transport it in height, control block or outer Hardy accumulation,
establish `Lambda<=0`, PF-infinity, RH, or a prize-level theorem.
