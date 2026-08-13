# Complete-Theta Positive-Time Spatial Tile Degree Certificate

Date: 2026-08-04

Status: rigorous local complete-theta contact exclusion. This is not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate.py --progress
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate.py
```

## Theorem

The complete theta-kernel heat transform and its first x derivative have no
common zero in

```text
0<=t<=0.45,  135<=x<=136.5.
```

The two new boundary covers use `64` cells
and have minimum separation `[1.560952109220542532101137559922580376403724570329703322323981650085889e-21 +/- 1.82e-91]`.
Both endpoint polygons have winding zero. The all-multiplicity same-sign
index theorem excludes every interior contact, and the central parent fills
the interval between the two new tiles.

## Proof Boundary

This is a rigorous complete-theta contact exclusion only on the declared positive-time band 0<=t<=0.45 and 135<=x<=136.5. It does not control other frequencies, negative heat time, all real x, the global de Bruijn-Newman constant, a degree-uniform Jensen remainder, Lambda<=0, RH, or a prize-level conclusion.
