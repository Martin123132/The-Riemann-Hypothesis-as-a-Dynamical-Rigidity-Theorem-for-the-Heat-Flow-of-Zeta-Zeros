# Rigorous Arb source-q dependency cut

Date: 2026-08-09

Status: finite saved-height source-q dependency cut closed; not a proof of the outer Hardy theorem, height uniformity, or RH

The corrected-model cell evaluator was audited function by function.  Its
only access to the chain record is through the exact coefficient `levels`.
It does not read emitted `q_terms`, recurrence `steps`, `qq_hex`, `t1`--`t5`,
saved PSI/ERF payloads, source recurrence states, or source multipliers.
Mathematical pi enters through `arb.pi()`, while `erf` and `digamma` are
evaluated independently by Arb/FLINT rather than imported from the source.

As an executable falsification check, chains 626, 653, 3406, 3423
were replayed at 110 decimal digits with a dictionary whose only
top-level key was `levels`.  All 4 replayed
correction balls overlap the saved independent cell certificates.

For all 1040 later recursive cells, a source-q-free projection of
the accepted cell, Arb correction, and analytic complete-tail fields has SHA-256
`e95ac19315fad25140e1e97e4e56cdc8044d7d8c45a861ab004f809a54b18e33`.  Exact nonnegative outer transport was
then rebuilt directly as

```text
block-20 analytic correction
  + direct nonrecursive kernel column
  + later Arb cell correction
  + later analytic complete tail.
```

No source-q point column or signed cross-call cancellation occurs in this
recombination.  The largest complete output is
`1.18534368330532563659717863467183140923914676684614E-3` and the minimum
margin below `0.005` is `3.81465631669467436340282136532816859076085323315386E-3`.
All 15 outputs close.

This proves a dependency classification, not RH.  Emitted source `q`, PSI/ERF
values, and recurrence states are not mathematical antecedents of this finite
corrected-model Arb bound.  They remain relevant only to the separate task of
reproducing or proving the legacy executable itself.  The direct-kernel column
still audits nonrecursive source arithmetic, and the outer Hardy representation
and remainder, coefficient coverage beyond the saved height, a height-uniform
route theorem, `Lambda<=0`, PF-infinity, RH, and a prize-level proof remain open.
