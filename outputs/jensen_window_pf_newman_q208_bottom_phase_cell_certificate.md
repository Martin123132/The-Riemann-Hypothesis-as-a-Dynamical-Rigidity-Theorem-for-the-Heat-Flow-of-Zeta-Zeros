# Newman Q208 Bottom Phase-Cell Certificate

Date: 2026-07-25

Status: rigorous complete fixed-time bottom certificate.
This is not a proof of Q208, `Lambda<=0`, RH, or the Clay prize.

```text
work/rh_compute/results/jensen_window_pf_newman_q208_bottom_phase_cell_certificate.json
work/rh_compute/results/jensen_window_pf_newman_q208_bottom_phase_cell_certificate.jsonl
python work/rh_compute/scripts/jensen_window_pf_newman_q208_bottom_phase_cell_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_q208_bottom_phase_cell_certificate.py
```

## Contract

The certified proxy is

```text
F_t(x)=16*(1+x^4)*H_t(x),
F_t'(x)=64*x^3*H_t(x)+16*(1+x^4)*H_t'(x).
```

Its triangular map from `(H,H')` has positive determinant and a
canonical homotopy to the identity. It therefore has the same
contact set and closed-boundary winding as the ordinary first jet,
while remaining nonsingular at `x=0` and asymptotic to `J=16x^4H`.

The retained transform uses six ordinary theta terms at the single
time `t=1/1040`, Taylor order 24 in `x`, and the independently
proved raw `E0(7),E1(7)` omitted-tail bars.

## Progress

```text
panels=492/492
certified cells=492
unresolved=0
branches={'value': 313, 'derivative': 179}
minimum closest-point norm lower=[2.96403625960278050652358262105431969370893790069555246382700767905657488040984531840749127695744391530752182006835937500e-29 +/- 5.62e-45]
maximum subdivision depth=0
stop reason=deferred_baseline_cpu
```

## Exact Open Phase Chain

Every cell encloses the whole path segment and excludes the
origin. Each adjacent pair has an exact dyadic rational
intersection witness. The resulting open polygon has

```text
cells=492
witnesses=493
positive-ray crossing count=-19
```

This crossing count is an open-edge contribution, not a
winding number. The right-corner witness must still be
matched to the transformed right-edge cells, and the top
edge must be traversed in reverse before the cyclic integer
can be evaluated.

## Scope

A complete bottom certificate does not by itself certify Q208.
The top phase chain, right-edge coordinate transfer, corner
intersections, and one exact closed polygon winding remain open.
