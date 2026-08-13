# Four-Summand Positive-Time Boundary Scout

Date: 2026-08-03

Status: finite high-precision diagnostic, not an interval certificate and not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fourth_summand_positive_time_boundary_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_summand_positive_time_boundary_scout.py --progress
```

## Result

The sampled winding values are `[0]` on all
`40` coefficient nodes from `0.22` through `1`.
The weakest sampled boundary norm is `4.63425574118937870e-23`.
The largest polygon angle increment is `0.34542605` radians.

This scout decides whether the declared positive-time rectangle is a viable
target for a parameter-uniform ACB/Taylor boundary certificate. It does not
exclude an event between sampled vertices or coefficient nodes.

## Proof Boundary

Finite 192-bit point diagnostic on a sampled boundary and forty coefficient nodes only. It is not interval coverage between vertices or lambda nodes, not a Brouwer-degree theorem, not control outside the declared rectangle, not a fifth-summand or complete-Xi result, and not a proof of Lambda<=0 or RH.
