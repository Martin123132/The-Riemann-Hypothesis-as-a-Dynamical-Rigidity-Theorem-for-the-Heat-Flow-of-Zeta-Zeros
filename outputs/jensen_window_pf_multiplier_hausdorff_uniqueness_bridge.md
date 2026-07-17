# Jensen-Window PF Multiplier Hausdorff-Uniqueness Bridge

Date: 2026-07-11

Status: exact Hausdorff uniqueness and unit-atomic characterization. This is not a proof
of the counting-measure target, PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_multiplier_hausdorff_uniqueness_bridge`.

```text
work/rh_compute/results/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.json
python work/rh_compute/scripts/jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py
```

Current result:

```text
validated Jensen-window PF multiplier Hausdorff-uniqueness bridge: 10 rows, 0 issues, 1 Hausdorff measure theorem, 1 unit-atomic characterization, 1 interpolation guard, 1 open recovery handoff
```

## Exact Measure Bridge

For `z>1`, Frullani's identity gives

```text
-log(1-z^(-2))
 =integral_0^infinity exp(-z*t)*q(t)dt,
q(t)=2*(cosh(t)-1)/t>0.
```

If the target multiset exists, put

```text
S(t)=sum_j exp(-alpha_j*t).
```

The convergence condition `sum_j alpha_j^(-2)<infinity` makes the
following positive measure finite. With `r=exp(-t)`,

```text
y_k=-log(x_k)=integral_[0,1] r^(k-1)dnu(r),
dnu(r)=q(-log r)*S(-log r)dr.
```

The values `y_1,y_2,...` are every moment of `nu`. Polynomials are dense
in `C([0,1])`, so the finite positive measure `nu` is unique. Consequently

```text
(-1)^m*Delta^m y_k
 =integral_[0,1] r^(k-1)*(1-r)^m dnu(r)>=0.
```

More importantly, the integer-only multiplier target is equivalent to the
unique measure having an absolutely continuous density of the rigid form

```text
dnu/dr=q(-log r)*sum_j r^alpha_j
```

with unit integer multiplicities and `sum_j alpha_j^(-2)<infinity`.
If such a representation exists, division by `q(t)>0` followed by Laplace
uniqueness also makes the multiset `{alpha_j}` unique.

## Interpolation Guard

Hausdorff uniqueness does not by itself identify the natural continuous
Mellin interpolation used by the power-sum obstruction. The exact ambiguity
already appears in

```text
f(s)=sin(2*pi*s),
f(n)=0,
f(s+1)-2*f(s)+f(s-1)=0.
```

Thus integer values, and even the same centered continuous log curvature,
permit a nontrivial periodic interpolation factor. A canonical-interpolation
or growth theorem is still needed before the continuous Mellin determinant
can rule out the integer-only product.

## New Handoff

The target is now a unique recovery problem:

```text
{y_k} -> unique Hausdorff measure nu
      -> test absolute continuity
      -> divide by q
      -> recover or rule out a unit counting Laplace measure.
```

The order-55 Arb frontier supports positivity of `nu` finitely but supplies
neither the all-order measure theorem nor the unit-multiplicity recovery.

```text
outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md
outputs/jensen_window_pf_multiplier_counting_measure_target.md
outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md
outputs/signed_hankel_jensen_dependency_graph.md
```
