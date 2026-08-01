# Newman Convex Phase-Cell Unwrapping Lemma

Date: 2026-07-25

Status: exact reusable phase/winding certificate with a rigorous
Q208 right-edge calibration. This is not a proof of Q208,
`Lambda<=0`, RH, or the Clay prize.

```text
work/rh_compute/results/jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma.json
python work/rh_compute/scripts/jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma.py
```

## Convex Cell Lemma

Let a closed continuous path be split into finitely many pieces,
with the kth piece enclosed by a compact convex set `C_k` that
does not contain the origin. If `p_k` is the closest point of
`C_k` to the origin, then

```text
p_k dot v >= ||p_k||^2 > 0  for every v in C_k.
```

Thus every cell lies in a strict open half-plane through the
origin and has one argument branch of width less than `pi`.

Convexity gives a relative-endpoint homotopy from each path
piece to its endpoint chord. If an exact rational point `q_k`
lies in `C_(k-1) intersect C_k` for every cyclic join, those
chords are in turn homotopic to the rational polygon through
the `q_k`. Therefore the continuous path and witness polygon
have exactly the same winding.

## Exact Winding

The polygon integer is evaluated by signed crossings of the
positive real ray:

```text
+1: p_y <= 0 < q_y and det(p,q) > 0
-1: q_y <= 0 < p_y and det(p,q) < 0.
```

All coordinates and determinants are rational. No floating-point
`atan2`, branch cut, or tolerance enters the winding integer.
A cell containing the origin, a missing adjacent-cell witness,
or an edge through the origin is rejected as unresolved.

The checker exercises exact winding `+1`, `-1`, and `0` examples
and independent rejection cases.

## Q208 Right Edge

The stored six-term Arb/Taylor theorem has been replayed at the
geometric interface and transferred to the globally regular
proxy `F=16(1+x^4)H`:

```text
right strip: [1/1040,1/5]x[245,246]
strip cells: 40
right edge: x=246, 1/1040<=t<=1/5
edge cells: 20
minimum full right-edge derivative lower bound: [2.74842995422409717150882800350618905119382156826663560755859516271491876751094511518300491904084419573354e-29 +/- 4.09e-134]
minimum transformed proxy derivative lower bound: [2.74842994114538475612679723538771225333469368077701874652160281001741346414414028782485097030428402853945e-29 +/- 1.85e-134]
```

Every right-edge cell has both `J_t'(246)>0` and
`F_t'(246)>0`, and the 20 time cells form a contiguous cover
from `1/1040` to `1/5`. Hence `F_t(246)+i F_t'(246)` stays in
the open upper half-plane and admits one principal phase branch.
This calibrates the method but does not determine the complete
closed-boundary winding.

## Live Handoff

The next certificate must cover the entire fixed-time bottom
edge `t=1/1040`, `0<=x<=246` by origin-free first-jet cells,
store exact rational witnesses at every join, compose them with
the axis, top, and right pieces, and compute one exact cyclic
polygon winding. Until then Q208 remains open.
