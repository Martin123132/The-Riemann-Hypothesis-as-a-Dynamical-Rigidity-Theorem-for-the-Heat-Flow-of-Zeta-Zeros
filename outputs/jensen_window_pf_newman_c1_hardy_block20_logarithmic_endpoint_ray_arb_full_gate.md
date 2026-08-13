# Hardy block-20 logarithmic endpoint-ray full Arb gate

Date: 2026-08-07

Status: rigorous all-374 finite-input logarithmic endpoint-ray enclosure; not a proof of a height-uniform recurrence or RH

## Scope

The four-branch pilot is extended to every recursive block-20 call.  Each
chain is evaluated atomically at 70 and 110 decimal digits and fsynced to an
append-only JSONL cache before the next chain begins.  The run remains serial
with one below-normal RH worker.

For each call the ray cutoff is the smallest power-of-two multiple of 64 for
which the exact polynomial-times-exponential tail bound is below `1e-30`:

```text
{"1024": 1, "128": 13, "256": 4, "512": 1, "64": 355}
maximum selected cutoff: 1024
```

No numerical tail at the cutoff is assumed to be zero.

## Exact identities

All 374 calls have 70/110-digit overlap for the four
rays and complete nonsaddle value.  On every call both rigorous complex
identities contain zero:

```text
Q_exact^- - (P^- - I69^-) = 0,
(Q_exact^- - Q_paper^-) - R_paper = 0.                 (1)
```

Aggregate enclosure bounds are

```text
maximum omitted-origin radius       <= 7.59133194080543378956362172776120227245461839142243E-15
maximum far-tail radius             <= 5.86083000123935813095248830503277364858077013704702E-31
maximum |formula-target gap|        <= 4.15547487074219225496041207890129953739233314990997E-14
maximum |correction-residual gap|    <= 4.15547487603614817529978919807831516664009541273117E-14
minimum |Q_exact-Q_paper|            >= 1.31747793749129559667889767787113397545795123733114E-6
maximum |Q_exact-Q_paper|            <= 1.53020780144673485725912958211984704615105101765948E-3
nonzero exact corrections             374 / 374
```

Thus the residual left after exact equation-(69) integration and published
W2--W5 evaluation is independently reconstructed as the exact full-cubic
nonsaddle correction on all 374 calls.  It is not a fitted term and is not a
source patch.

## Interpretation

This closes the finite block-20 contour-integration obligation: the exact
nonsaddle value, the direct parent-minus-saddle value, and the previously
recorded published residual agree as rigorous complex balls.  The remaining
proof problem is no longer numerical identification of this finite residual.
It is to derive a height-uniform analytic bound for the full-cubic versus
quadratic endpoint model and propagate that bound through the complete
recurrence and outer Hardy representation.

## Pi provenance and boundary

Pi is inherited from the paper's Fourier exponential and is evaluated as the
mathematical Arb constant.  It controls the exact logarithmic kernel, phase
and tail decay; no fitted normalization is inserted.

This is rigorous finite-input work for one low-height block.  It proves no
height-uniform W2--W4 estimate, recursive/global error accumulation, outer
Hardy bound, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
