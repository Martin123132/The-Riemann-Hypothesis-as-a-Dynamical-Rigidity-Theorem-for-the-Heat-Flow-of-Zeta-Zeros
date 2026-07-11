# Jensen-Window PF Heat-Flow Ratio-Cone Invariance Lemma

Date: 2026-07-06

Status: exact conditional ratio-cone invariance lemma. This is not a proof of
the monotone-contraction theorem for the zeta coefficients, Jensen-window
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_ratio_cone_invariance_lemma`.

Proof boundary: this artifact proves inward-pointing boundary algebra for the
ratio cone only, conditional on positive smooth coefficients and the required
infinite or collared finite ratio variables. It does not prove that the actual
zeta coefficients enter the cone at an initial or terminal lambda, does not
prove all-shift monotone contractions, does not prove `jwpf_06`, and does not
establish `Lambda <= 0`.

Machine-readable lemma:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py
```

Current result:

```text
validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues
```

## Cone

Use:

```text
r_k=A_{k+1}/A_k
x_k=r_k/r_{k-1}
a_k=(2*k-1)/(2*k+1)
```

The ratio cone is:

```text
a_k <= x_k <= 1
x_{k+1} >= x_k
```

The exact threshold lemma supplies the lower wall:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

## Boundary Algebra

Write

```text
dlog(x_k)/dlambda = 2*r_k*F_k
F_k=(2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1)
```

Lower wall:

```text
at x_k=a_k, F_k=(2*k+3)*x_{k+1}-(2*k+1)>=0
```

provided `x_{k+1}>=(2*k+1)/(2*k+3)`.

Upper wall:

```text
at x_k=1, F_k=(2*k+3)*(x_{k+1}-1)<=0
```

provided `x_{k+1}<=1`.

Monotone wall:

```text
at x_{k+1}=x_k=q,
dlog(x_{k+1}/x_k)/dlambda >= 0
```

provided `x_{k+2}>=q` and `q >= (2*k-1)/(2*k+5)`.

Thus the cone has exact forward-invariance boundary algebra. A finite prefix
needs collar variables: `x_{K+1}` for lower/upper wall tests and `x_{K+2}`
for monotone wall tests.

## Finite Sanity Diagnostics

The existing coefficient enclosure cache gives:

```text
lambdas: 0.0, 0.000001, 0.0001, 0.01, 0.1
max coefficient index: 64
coordinate lower rows: 315
coordinate upper rows: 315
coordinate monotone rows: 310
lower wall rows: 310
upper wall rows: 310
monotone wall rows: 305
```

Minimum coordinate margins:

```text
lower: 8.170843282044743874E-3 at lambda=0.1, k=63
upper: 7.573304754109889330E-3 at lambda=0.0, k=63
monotone: 9.895636183563448458E-5 at lambda=0.0, k=62
```

Minimum wall margins:

```text
lower wall: 1.037697096819682472E+0 at lambda=0.1, k=62
upper wall: 2.951966135520726476E-1 at lambda=0.0, k=1
monotone wall: 2.665474275226736438E-4 at lambda=0.1, k=61
```

## Remaining Gap

```text
The heat-flow ratio cone has an exact conditional forward-invariance algebra: the lower wall, upper wall, and monotone wall all point inward under the stated neighbor and threshold hypotheses. The remaining theorem gap is no longer the boundary algebra; it is to prove that the actual zeta heat-flow coefficient sequence enters the full infinite cone at a suitable lambda and that the infinite/collared flow argument is analytically legitimate.
```

Integration:

```text
outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```
