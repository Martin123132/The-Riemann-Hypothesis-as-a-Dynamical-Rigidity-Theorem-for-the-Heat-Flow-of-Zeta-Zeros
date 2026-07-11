# Jensen-Window PF Negative-Lambda Defect-Recurrence Scout

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
defect-tail theorem, cone entry, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_defect_recurrence_scout`.

Proof boundary: this artifact tests recurrence-style sufficient
conditions against the certified finite prefix. It does not prove an
all-`k` recurrence or tail theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_defect_recurrence_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_recurrence_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues
```

## Sufficient Condition

A finite-compatible sufficient tail theorem is:

```text
0 <= d_k <= 2/(3*(2*k+1))
d_(k+1) <= d_k
```

This is stronger than the defect-tail target, but it is compatible with
the certified prefix and would imply the required barriers.

## Rejected Recurrence

The tempting stepwise wall-preserving recurrence

```text
d_(k+1) <= ((2*k+1)/(2*k+3))*d_k
```

is rejected on every checked adjacent pair.

Finite diagnostics:

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_22
checked contractions: x_1..x_21
buffered sufficient rows: 63 / 63
defect monotone rows: 60 / 60
width-preserving recurrence rows: 0 / 60
width-preserving rejected rows: 60 / 60
```

Extrema:

```text
min buffer margin: 4.569139109160357364E-2 at lambda=-25.0, k=21
min defect monotone gap: 7.539431011168075104E-5 at lambda=-100.0, k=20
max width-preserving excess: 1.062756314392564011E-2 at lambda=-25.0, k=1
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
```

Summary:

A buffered sufficient tail theorem is compatible with the finite prefix: prove 0<=d_k<=2/(3*(2*k+1)) for k>=22 and d_(k+1)<=d_k for k>=21. The tempting width-preserving one-step recurrence d_(k+1)<=((2*k+1)/(2*k+3))*d_k is rejected on all checked adjacent pairs, so it should not be used as the direct recurrence route.
