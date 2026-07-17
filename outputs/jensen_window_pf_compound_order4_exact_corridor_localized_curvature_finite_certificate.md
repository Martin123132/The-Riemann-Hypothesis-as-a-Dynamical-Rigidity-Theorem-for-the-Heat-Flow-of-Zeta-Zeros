# Jensen-Window PF Order-Four Exact-Corridor Localized-Curvature Finite Certificate

Date: 2026-07-13

Status: rigorous exact-corridor localized-curvature theorem on `2<=u<=20`.
This is not a proof of the remaining asymptotic ray, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.json
work/rh_compute/results/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_chunks.jsonl
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.py
```

## Correlated Cover

The cover uses width `10^-4` on `2<=u<=2.3` and width `10^-3` on
`2.3<=u<=20`. Each central block has a rational outer pad whose endpoint
potential slopes certify the complete `t+-2` collar. The exact corridor
boxes are evaluated with the same interval mode as `q`, `t`, curvature,
and the Hurwitz-zeta derivatives, preserving their correlations.

```text
mode blocks:                  20700
t-collar endpoint gates:      41400
maximum t^2*U upper:          3.49148869585704706421401378217589691306143341465438292500535E+0
minimum scaled margin:        8.511304142952935785986217823999999999999999999999999999999999999999999E-3
minimum localized J lower:    2.80129194625256811822821472708524207910912247123837821877901E-37
```

Every block proves `j_0>E_0` and `t^2 U(t)<7/2`. Therefore

```text
K_1(t)<=7/(2t^2) for every mode 2<=u<=20.
```

Composed with the earlier compact theorem, the ceiling now holds from
`t=319` through `t=V'(20)`.

## Remaining Boundary

The only curvature segment left is `u>=20`. The finite cover is not
extrapolated. The next theorem must use q-leading mode geometry and the
exact corridor constants to prove the localized inequality analytically
on that full ray.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/formal_core.md
```
