# Four-Summand Positive-Time Degree Certificate

Date: 2026-08-03

Status: rigorous finite-homotopy contact-exclusion theorem. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fourth_summand_positive_time_degree_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_summand_positive_time_degree_certificate.py --progress
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_positive_time_degree_certificate.py
```

## Compact Theorem

For every `0.22<=lambda<=1`, the finite four-term transform has no
common zero of `(F_lambda,F_lambda,x)` in

```text
0<=t<=0.45,  135.5<=x<=136.
```

The lateral-boundary cover has `1248` cells;
its minimum one-component separation is `[4.330285654244989032542320539192868800659758091418775408939610330095095e-23 +/- 4.02e-93]`.
At lambda `0.22`, `32` convex image segments
give exact `(t,x)` winding zero. Uniform boundary nonvanishing preserves
that degree through the full fourth coefficient.

Every possible heat contact has local index `-floor(m/2)` in `(t,x)`.
A zero total degree therefore permits no regular or degenerate contact.

## Proof Boundary

This excludes contacts only for the finite homotopy H_1+H_2+H_3+lambda H_4 inside the declared positive-time rectangle and only for 0.22<=lambda<=1. It does not control contacts outside that rectangle, negative heat time, summands n>=5, the complete theta/Xi transform, an Xi collision, a degree-uniform Jensen remainder, Lambda<=0, RH, or a prize-level conclusion.
