# Jensen-Window PF Monotone-Contraction Stress

Date: 2026-07-06

Status: finite Arb ratio-curvature stress. This is not a proof of an analytic
monotone-contraction theorem, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_monotone_contraction_stress_summary`.

Proof boundary: this artifact checks adjacent log-concavity and increasing
ratio contractions on the existing finite zeta coefficient enclosure cache. It
does not cover all shifts beyond the cache, all lambda values, all Schur
shapes, or the missing Newman-direction theorem.

Machine-readable summary:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json
```

Row-level JSONL:

```text
work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64.jsonl
```

Builder:

```text
python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_stress.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_stress.py
```

Current result:

```text
validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues
```

## Condition

For positive coefficients `A_k(lambda)`, define:

```text
x_k = (A_{k+1}/A_k) / (A_k/A_{k-1})
```

Adjacent log-concavity is `0 < x_k <= 1`. The new frontier condition is:

```text
x_{k+1} >= x_k
```

Equivalently:

```text
A_{k+2}*A_k**3 >= A_{k+1}**3*A_{k-1}
```

So the stress test is probing a finite third-log-difference sign pattern of
the zeta heat-flow coefficient sequence.

## Grid

The run uses the existing coefficient enclosures:

```text
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl
work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl
```

It checks degrees d=3..12 using the k<=64 cache. For each degree `d`, it
checks every shift `0 <= n <= 64-d`, on the five lambda values in the cache.
This gives 2875 finite Arb-classified rows.

## Result

Every checked row satisfies:

```text
adjacent_log_concavity_ok = true
monotone_contractions_ok = true
contains_zero = false
ok = true
```

Counts by degree:

```text
d=3: 310/310
d=4: 305/305
d=5: 300/300
d=6: 295/295
d=7: 290/290
d=8: 285/285
d=9: 280/280
d=10: 275/275
d=11: 270/270
d=12: 265/265
```

The global min monotone gap sample is:

```text
9.895636183563448458E-5
```

The global min `1-x_k` sample is:

```text
7.573304754109889330E-3
```

These samples are for orientation only; pass/fail classification uses Arb
balls.

## Relation To Frontier Scout

This broadens the finite evidence behind:

```text
outputs/jensen_window_pf_monotone_contraction_frontier_scout.md
```

The frontier scout showed exact positivity of the first two hard
column-frontier polynomials under `x1 <= x2 <= x3`. This stress test shows
that the same increasing-contraction pattern holds across a much wider finite
zeta grid. It still leaves the analytic theorem open:

```text
Prove noncircularly that the zeta heat-flow coefficients satisfy the required
increasing-contraction inequalities in the all-shift/all-lambda regime needed
for the Jensen-window PF bridge.
```

## Boundary

This stress artifact is useful evidence for a theorem target. It is not a
determinant-integral construction, not an all-shape Cauchy-Binet kernel, not a
Laguerre-Polya theorem, and not a proof of `Lambda <= 0`.
