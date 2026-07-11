# Jensen-Window PF Negative-Lambda Gaussian Curvature Matrix

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_gaussian_curvature_matrix`.

Proof boundary: this artifact compares the actual negative-lambda prefix
with the Gaussian raw-moment baseline. It does not prove any all-`k`
tail theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda Gaussian curvature matrix: 7 matrix rows, 63 positive-deficit rows, 63 bounded-deficit rows, 63 raw-threshold rows, 0 issues
```

## Baseline

For Gaussian raw moments:

```text
M_k^G proportional to Gamma(k+1/2) T^(-k-1/2)
M_(k+1)^G*M_(k-1)^G/(M_k^G)^2 = (2*k+1)/(2*k-1)
x_k = 1
B_k = -log(x_k) = 0
```

Thus the bounded log-curvature target asks for a small positive deficit
from the Gaussian baseline:

```text
0 < B_k <= 2/(3*(2*k+1))
```

A positive Gaussian scale mixture is not the right upper-wall proof
template, because it adds nonnegative curvature to the Gaussian baseline
and pushes `x_k>=1` rather than `x_k<=1`.

Finite diagnostics:

```text
lambdas: -100.0, -50.0, -25.0
coefficient range: A_0..A_22
checked contractions: x_1..x_21
positive Gaussian-deficit rows: 63 / 63
bounded Gaussian-deficit rows: 63 / 63
raw threshold rows: 63 / 63
positive threshold rows: 63 / 63
```

Extrema:

```text
min deficit-bound slack: 2.034879929997686858E-3 at lambda=-25.0, k=21
min raw-threshold margin: 2.034879929997686858E-3 at lambda=-25.0, k=21
max Gaussian deficit: 3.524954269540071057E-2 at lambda=-25.0, k=1
max deficit fraction of bound: 8.687502445151491977E-1 at lambda=-25.0, k=21
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
```

Summary:

The negative-lambda curvature target is a controlled deficit from the Gaussian moment baseline: Gaussian moments give raw curvature log((2*k+1)/(2*k-1)) and x_k=1, while the actual zeta prefix has a positive deficit B_k=-log(x_k) on 63 rows but keeps it below 2/(3*(2*k+1)) on all 63 checked rows. Positive Gaussian scale-mixture arguments point the wrong way for the upper wall because they make the deficit nonpositive.
