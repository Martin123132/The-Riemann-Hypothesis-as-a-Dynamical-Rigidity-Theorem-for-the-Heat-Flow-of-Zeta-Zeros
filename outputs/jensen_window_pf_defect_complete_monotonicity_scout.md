# Jensen-Window PF Defect Complete-Monotonicity Scout

Date: 2026-07-10

Status: finite interval diagnostic with exact countermodel gate. This is not a proof of
all-k complete monotonicity, PF-infinity, Jensen hyperbolicity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_defect_complete_monotonicity_scout`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_defect_complete_monotonicity_scout.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_defect_complete_monotonicity_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_defect_complete_monotonicity_scout.py
```

Current result:

```text
validated Jensen-window PF defect complete-monotonicity scout: 3284 defect positives, 3288 log positives, 838 inconclusive, both certified through order 8, 5 lambdas, 1 exact all-shape countermodel, 0 issues
```

## Finite Arb Pattern

For `d_k=1-x_k`, every cached interval through order 8 satisfies

```text
(-1)^m*Delta^m d_k > 0, 0<=m<=8.
```

Across all orders 0 through 12 there are `3284` certified
positive intervals and `421` intervals containing zero; no interval is
strictly negative. The order 9 through 12 midpoints remain positive, but are not promoted.
This is finite interval evidence only; no all-`k` theorem is claimed.

For the multiplier-kernel variable `y_k=-log(x_k)`, every cached interval through
order 8 also satisfies

```text
(-1)^m*Delta^m y_k > 0, 0<=m<=8.
```

Across orders 0 through 12, `3288` log intervals
are certified positive and `417` contain zero; none is
strictly negative. This is the finite necessary pattern for a positive multiplier-kernel
or counting-measure representation, not a proof that such a representation exists.

## Exact Guard

Take the Hausdorff defect measure `(1/2)*delta_0`. Then

```text
d_1=1/2 and d_k=0 for every k>=2
x_1=1/2 and x_k=1 for every k>=2
1/3 <= x_1 <= x_2 <= ... <= 1
```

so the infinite sequence satisfies the full static ratio cone and complete defect
monotonicity. Its degree-3 Jensen polynomial is

```text
1+3*z+(3/2)*z^2+(1/4)*z^3
discriminant = -27/16 < 0.
```

The full static ratio cone plus complete monotonicity of d_k does not imply all-shape Jensen-window PF or degree-3 hyperbolicity.

## Consequence

Every cached zeta defect difference through order 8 is strictly positive after the alternating sign; orders 9 through 12 retain positive midpoints but include 421 enclosure-width inconclusives. The multiplier-kernel sequence y_k=-log(x_k) is likewise certified through order 8, with 3288 positive intervals and 417 high-order inconclusives. An exact Hausdorff-defect sequence inside the full static ratio cone has a degree-3 Jensen discriminant -27/16. Complete defect monotonicity is therefore a live column-structure diagnostic, not an all-shape bridge.
The next viable use is a column-recurrence theorem or a stronger heat-flow-specific
hierarchy; this condition cannot be promoted directly to the all-shape bridge.
