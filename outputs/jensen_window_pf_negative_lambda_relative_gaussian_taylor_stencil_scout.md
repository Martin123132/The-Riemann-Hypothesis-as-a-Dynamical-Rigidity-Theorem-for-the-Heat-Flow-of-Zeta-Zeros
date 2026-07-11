# Jensen-Window PF Negative-Lambda Relative-Gaussian Taylor Stencil Scout

Date: 2026-07-07

Status: finite theorem-search diagnostic. This is not a proof of
scaled-curvature monotonicity, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout`.

Proof boundary: this artifact certifies fixed-`k` leading Taylor
signs and finite truncation behavior for the relative-Gaussian
stencil. It does not prove a uniform Taylor-tail remainder theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: 8 rows, 0 issues, 3 certified leading-sign rows, 35 truncation rows, 4 all-positive stencil rows, 0 ready-to-apply rows
```

## Fixed-k Leading Signs

For the normalized local multiplier:

```text
F_k = 1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+...
q = k+1/2
f_k = log(F_k)
```

the formal fixed-`k` stencil gives:

```text
B_k = (a^2-2*b)/T^2 + O_k(T^-3)
B_k-B_(k+1) = 2*(a^3-3*a*b+3*c)/T^3 + O_k(T^-4)
C_(k+1)-C_k = 2*(a^2-2*b)/T^2 + O_k(T^-3)
```

Certified leading coefficients:

```text
a = [-37.45380985900664233683260390434643972954 +/- 2.44e-39]
b = [606.0506446363092395662753197700082678903 +/- 9.11e-39]
c = [-5047.989249594588274331383288623799728021 +/- 3.86e-37]
a^2-2*b = [190.686583682004682389944962471127216573 +/- 6.58e-37]
2*(a^3-3*a*b+3*c) = [825.9976249272062741326759643225349915 +/- 4.60e-35]
2*(a^2-2*b) = [381.373167364009364779889924942254433147 +/- 4.86e-37]
```

## Truncation Stencil Matrix

Rows classify degrees 6, 8, 10, 12, and 14 at `k=22` and
`T=25,50,100,200,500,1000,2000`:

```text
truncation rows: 35
invalid normalizers: 9
positive normalizers: 26
B-positive rows: 24
B-decrease rows: 18
C-increase rows: 12
all-positive stencil rows: 4
upper-wall violation rows: 2
companion-failure rows: 8
weighted-gap failure rows: 12
```

Interpretation:

The leading fixed-`k` signs point in the right direction, but finite
Taylor truncations are not stable proof objects. A promoted proof must
bound the full Taylor tail in the `B_k`, `B_k-B_(k+1)`, and
`C_(k+1)-C_k` stencils over a stated `q/T` range.

Integration:

```text
outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md
outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md
outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
```

Summary:

The fixed-k Taylor stencil has certified positive leading signs for B_k, B_k-B_(k+1), and the relative-Gaussian weighted gap C_(k+1)-C_k. However, finite truncations through degree 14 are not stable proof objects: only 4/35 sampled truncation rows satisfy all three stencil inequalities, so the route still needs a uniform Taylor-tail remainder theorem.
