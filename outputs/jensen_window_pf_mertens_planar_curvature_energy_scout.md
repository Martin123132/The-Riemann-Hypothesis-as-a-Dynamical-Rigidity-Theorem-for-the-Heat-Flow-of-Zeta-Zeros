# Jensen-Window PF Mertens Planar Curvature-Energy Scout

Date: 2026-07-24

Status: finite double-precision scaling diagnostic for the open
planar curvature-energy gate; not a proof artifact.
This diagnostic does not prove an asymptotic theorem, RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json
python work/rh_compute/scripts/jensen_window_pf_mertens_planar_curvature_energy_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py
```

## Purpose

The exact handoff
`outputs/jensen_window_pf_mertens_planar_abel_handoff.md` proves

```text
V_(alpha,K)<3*pi^2*(1+log(2R))/K,
|O_(alpha,K)|^2<=V_(alpha,K)E_(alpha,K),
```

and leaves `E_(alpha,K)=O_epsilon(K^(1+epsilon))` open. This
scout asks only whether the finite actual-Mobius values immediately
falsify that scale.

## Grid

| K | alpha | R | K V_K | E_K/K | M_K/K | O_K | Abel error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 0.125 | 14 | 4.76092841 | 0.697715829 | 1.5625 | -0.632270499 | 1.11e-16 |
| 16 | 0.25 | 12 | 4.79043638 | 0.680267373 | 1.5625 | -0.61826626 | 2.22e-16 |
| 16 | 0.5 | 8 | 4.88583275 | 0.654116549 | 1.5625 | -0.594845505 | 2.22e-16 |
| 16 | 0.75 | 6 | 5.03882852 | 0.634584804 | 1.5625 | -0.572119646 | 1.11e-16 |
| 32 | 0.125 | 26 | 5.22164387 | 0.629035544 | 1.8125 | -0.635407742 | 2.22e-16 |
| 32 | 0.25 | 21 | 5.33302036 | 0.603172744 | 1.8125 | -0.622022631 | 0 |
| 32 | 0.5 | 14 | 5.53154216 | 0.557968463 | 1.8125 | -0.596867655 | 1.11e-16 |
| 32 | 0.75 | 9 | 5.71159561 | 0.521251789 | 1.8125 | -0.578482258 | 1.11e-16 |
| 64 | 0.125 | 50 | 5.48323561 | 0.78955672 | 1.640625 | -0.611874369 | 1.11e-16 |
| 64 | 0.25 | 39 | 5.65209823 | 0.765636508 | 1.640625 | -0.598719211 | 0 |
| 64 | 0.5 | 23 | 5.93534634 | 0.725280746 | 1.640625 | -0.574828911 | 2.22e-16 |
| 64 | 0.75 | 14 | 6.10102256 | 0.69415308 | 1.640625 | -0.552877888 | 1.11e-16 |
| 128 | 0.125 | 95 | 5.65270339 | 0.731044486 | 1.28125 | -0.622326147 | 1.11e-16 |
| 128 | 0.25 | 70 | 5.88624881 | 0.714076224 | 1.28125 | -0.609407038 | 1.11e-16 |
| 128 | 0.5 | 39 | 6.17728323 | 0.682145553 | 1.28125 | -0.585908085 | 2.22e-16 |
| 128 | 0.75 | 21 | 6.33275466 | 0.65675789 | 1.28125 | -0.565009797 | 2.22e-16 |
| 256 | 0.125 | 182 | 5.76062921 | 0.777975139 | 1.7578125 | -0.643531086 | 0 |
| 256 | 0.25 | 128 | 6.0227647 | 0.758867493 | 1.7578125 | -0.630248061 | 2.22e-16 |
| 256 | 0.5 | 64 | 6.33009502 | 0.722807075 | 1.7578125 | -0.605717928 | 1.11e-16 |
| 256 | 0.75 | 32 | 6.46081191 | 0.693121114 | 1.7578125 | -0.583320326 | 1.11e-16 |
| 512 | 0.125 | 347 | 5.8442382 | 0.688396743 | 1.625 | -0.646249222 | 1.11e-16 |
| 512 | 0.25 | 235 | 6.14476202 | 0.670625165 | 1.625 | -0.632310783 | 0 |
| 512 | 0.5 | 108 | 6.4269777 | 0.636524269 | 1.625 | -0.606748187 | 3.33e-16 |
| 512 | 0.75 | 50 | 6.53109861 | 0.608157081 | 1.625 | -0.583442757 | 1.11e-16 |
| 1024 | 0.125 | 664 | 5.90801869 | 0.830405727 | 2.09570312 | -0.651938715 | 6.66e-16 |
| 1024 | 0.25 | 431 | 6.22805398 | 0.803596523 | 2.09570312 | -0.638486659 | 5.55e-16 |
| 1024 | 0.5 | 182 | 6.48468477 | 0.751810512 | 2.09570312 | -0.614032788 | 1.11e-16 |
| 1024 | 0.75 | 77 | 6.57030104 | 0.70857728 | 2.09570312 | -0.592205314 | 3.33e-16 |

## Finite Observation

Across this recorded grid:

- maximum `K*V_K`: `6.57030103865`
- maximum `E_K/K`: `0.830405726975`
- maximum `M_K/K`: `2.095703125`
- maximum direct/double-Abel discrepancy: `6.66133814775e-16`

Thus the tested values through `K=1024` are compatible with
`V_K=O(1/K)` and `E_K=O(K)` on this grid. This is a finite
observation, not evidence for a uniform constant and not an
extrapolation theorem. The proved analytic statement remains
`V_K=O((1+log R)/K)`, and the all-scale Mobius energy estimate
remains open.

## Boundary

The calculations use ordinary float64 arithmetic and finitely
many dyadic blocks. Passing the checker confirms reproduction
on selected small rows, the exact finite Abel identity to
floating tolerance, and honest status language. It does not
prove asymptotic scaling, the curvature-energy target, the full
Burnol bound, RH, PF-infinity, or `Lambda <= 0`.
