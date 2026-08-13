# Fourth-Summand Tangency Interval Continuation Certificate

Date: 2026-08-03

Status: rigorous interval continuation certificate for one finite theta homotopy.
It is not a proof of RH and does not represent the complete theta sum.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fourth_summand_tangency_interval_continuation_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_summand_tangency_interval_continuation_certificate.py --progress
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_tangency_interval_continuation_certificate.py
```

## Certified Statement

A unique regular contact branch in the certified local chain exists for 0<=lambda<=0.22, has F_xx<0 and dt/dlambda<0 throughout, and crosses t=0 exactly once in the attached crossing enclosure.

The certificate uses `11` overlapping parameter
charts and `10` strict endpoint connectors.
Every chart proves `H_4>0`, `F_xx<0`, and therefore
`dt/dlambda=H_4/F_xx<0` on its unique contact.

## Crossing

The unique `t=0` crossing lies in lambda ball `[0.2016916799221276872405642198031810373999722401475015637393482937481187 +/- 3.35e-16]`
and x ball `[135.8619892565584600777999051559123149576525618374852434672992302989559 +/- 7.38e-16]`.

## Method

Direct interval substitution erases the high-frequency cancellation. The certificate
instead evaluates point oscillatory moments with ACB, encloses each parameter box by a
time-order `10` and frequency-order `28` Taylor model,
adds positive-moment remainders and the explicit cutoff tail, and applies a uniform
two-dimensional Krawczyk inclusion on every lambda slab.

## Pi Provenance

pi is fixed by the Jacobi theta normalization sum_n exp(-pi*n^2*y), with y=exp(4u); changing pi changes the zeta/Xi kernel.

## Proof Boundary

This certifies one connected contact branch for the four-term finite theta homotopy. It does not count all contacts outside these boxes, control lambda>0.22, add n>=5, exclude a complete-Xi contact, prove Lambda<=0, RH, or the Clay theorem.
