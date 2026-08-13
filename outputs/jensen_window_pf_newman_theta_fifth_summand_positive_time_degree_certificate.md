# Fifth-Summand Positive-Time Degree Certificate

Date: 2026-08-04

Status: rigorous finite-homotopy contact-exclusion theorem. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fifth_summand_positive_time_degree_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fifth_summand_positive_time_degree_certificate.py --progress
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fifth_summand_positive_time_degree_certificate.py
```

## Compact Theorem

For every `0<=mu<=1`, the finite transform
`H_1+H_2+H_3+H_4+mu H_5` has no common zero of `(F,F_x)` in

```text
0<=t<=0.45,  135.5<=x<=136.
```

The lateral-boundary cover has `32` cells;
its minimum separation is `[1.946581008278572725404437041561276623575684854292304814661743878009124e-21 +/- 4.27e-91]`.
At mu zero, `32` convex image segments give
winding zero. Every possible contact has negative local index, so degree
zero excludes regular and degenerate contacts alike.

## Scaling Signal

The initial coefficient width is the whole interval `[0,1]`. The stored
fifth-component boxes and perturbation-to-margin ratios provide the first
calibration row for a uniform later-summand tail theorem.

## Proof Boundary

This excludes contacts only for the finite fifth-summand homotopy inside the declared local positive-time rectangle. It does not control other frequency rectangles, negative heat time, summands n>=6, the complete theta/Xi transform, an Xi collision, a degree-uniform Jensen remainder, Lambda<=0, RH, or a prize-level conclusion.
