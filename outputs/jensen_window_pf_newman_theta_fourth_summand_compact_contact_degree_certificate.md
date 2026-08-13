# Four-Summand Compact Contact Degree Certificate

Date: 2026-08-03

Status: rigorous finite-homotopy compact contact count. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate.py --progress
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate.py
```

## Compact Theorem

For every `0<=lambda<=0.22`, the tracked four-summand branch contact is
the unique common zero of `(F_lambda,F_lambda,x)` in

```text
-0.06<=t<=0.45,  135.5<=x<=136.
```

The uniform boundary cover has `352` interval cells;
its minimum one-component separation is `[1.516233069921239552187963644321759989481189835386458838485765283173776e-23 +/- 3.78e-93]`.
At lambda zero, `32` convex image segments give
the exact `(t,x)` boundary winding `-1`.

The all-multiplicity heat-contact theorem assigns local index
`-floor(m/2)` in `(t,x)`. The certified regular branch already contributes
`-1`; any hidden contact would make the total degree more negative.

## Why No Triple-Contact Search Is Needed

A separate exclusion of F=F_x=F_xx=0 is unnecessary because the all-multiplicity index theorem already assigns every triple or higher contact a nonzero charge of the same sign.

## Proof Boundary

This counts contacts only for the finite homotopy H_1+H_2+H_3+lambda H_4 inside the declared compact cylinder and only for 0<=lambda<=0.22. It does not control contacts outside that cylinder, lambda>0.22, summands n>=5, the complete theta/Xi transform, an Xi collision, a degree-uniform Jensen remainder, Lambda<=0, RH, or a prize-level conclusion.
