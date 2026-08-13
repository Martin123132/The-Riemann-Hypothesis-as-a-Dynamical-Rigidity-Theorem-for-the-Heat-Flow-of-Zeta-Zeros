# Complete-Theta Positive-Time Spatial Tile Boundary Scout

Date: 2026-08-04

Status: finite high-precision boundary diagnostic, not an interval tile certificate and not a proof of RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_complete_spatial_tile_boundary_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_complete_spatial_tile_boundary_scout.py --progress
```

## Result

The run requested `384` boundary
vertices and evaluated `318` unique
points, reusing `66` shared-edge references.

- `tile_0` on x=`135` to `135.5`: winding `0`, minimum sampled complete norm lower `1.93800824274295364e-21`, maximum angle increment `0.015792`.
- `tile_1` on x=`135.5` to `136`: winding `0`, minimum sampled complete norm lower `1.99105782921552717e-21`, maximum angle increment `0.010644`.
- `tile_2` on x=`136` to `136.5`: winding `0`, minimum sampled complete norm lower `1.90521590531241579e-21`, maximum angle increment `0.007401`.

The rigorous n>=6 derivative budget is applied at every sampled point,
but no claim is made between adjacent boundary vertices.

## Proof Boundary

Finite 192-bit point diagnostic on the sampled boundaries of the declared tiles. The n>=6 derivative budget is rigorous at each sampled point, but there is no interval coverage between vertices, no certified segment homotopy, and no Brouwer-degree theorem for the new tiles. It does not control other frequencies, negative heat time, all real x, the global de Bruijn-Newman constant, Lambda<=0, RH, or a prize-level conclusion.
