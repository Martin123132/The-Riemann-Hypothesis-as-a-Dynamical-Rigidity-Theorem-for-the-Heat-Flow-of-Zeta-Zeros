# Jensen-Window PF Negative-Lambda Log-Curvature Bridge

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
defect-tail theorem, cone entry, Jensen-window PF-infinity, RH, or
`Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_log_curvature_bridge`.

Proof boundary: this artifact translates the buffered defect tail into
log-curvature inequalities and checks finite compatibility. It does
not prove the required all-`k` log-curvature theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_log_curvature_bridge.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_log_curvature_bridge.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_log_curvature_bridge.py
```

Current result:

```text
validated Jensen-window PF negative-lambda log-curvature bridge: 63 simple log-buffer rows, 63 exact defect-buffer rows, 60 curvature-monotone rows, 5 bridge rows, 0 issues
```

## Log-Curvature Translation

Let:

```text
B_k = -Delta^2 log A_k = -log(x_k)
x_k = (A_{k+1}/A_k)/(A_k/A_{k-1})
d_k = 1 - x_k
```

A sufficient all-tail theorem for the buffered route is:

```text
0 <= B_k <= 2/(3*(2*k+1))
B_(k+1) <= B_k
```

The second line is the existing monotone-contraction sign target in
log-curvature language:

```text
Delta^3 log A_{k-1} >= 0
```

Thus the remaining new analytic burden is the bounded-curvature upper
estimate `B_k <= 2/(3*(2*k+1))` on the actual zeta heat-flow tail.

Finite diagnostics:

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_22
checked contractions: x_1..x_21
positive curvature rows: 63 / 63
exact defect-buffer rows: 63 / 63
simple log-buffer rows: 63 / 63
curvature-monotone rows: 60 / 60
```

Extrema:

```text
min simple log slack: 2.034879929997686858E-3 at lambda=-25.0, k=21
min exact log slack: 2.156321864086280617E-3 at lambda=-25.0, k=21
max log curvature: 3.524954269540071057E-2 at lambda=-25.0, k=1
min curvature-monotone gap: 7.539431011168075104E-5 at lambda=-100.0, k=20
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md
outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```

Summary:

The buffered defect-tail condition reduces to a cleaner log-curvature target: prove 0<=B_k<=2/(3*(2*k+1)) for k>=22, where B_k=-Delta^2 log A_k=-log(x_k), and prove B_(k+1)<=B_k from k>=21. The second clause is exactly the Delta^3 log A monotone-contraction target, while the bounded-curvature upper estimate remains open.
