# Jensen-Window PF Negative-Lambda Taylor Moment Budget

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_taylor_moment_budget`.

Proof boundary: this artifact derives the Gaussian-moment Taylor budget
for the local/mesoscopic route. It does not prove the Taylor-tail
remainder theorem or close the all-`k` negative-lambda tail.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_taylor_moment_budget.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_taylor_moment_budget.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py
```

Current result:

```text
validated Jensen-window PF negative-lambda Taylor moment budget: 9 budget rows, 7 tail-start samples, 4 invalid truncation rows, 2 bounded truncation rows, 0 ready-to-apply rows, 0 issues
```

## Moment Budget

For the normalized local Taylor multiplier:

```text
q = k+1/2
F_k = 1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+R_k^(>=8)
x_k = F_(k+1)*F_(k-1)/F_k^2
B_k = -log(x_k)
```

the expansion parameter is `q/T`, not fixed `k` alone.

Certified local ratios:

```text
a = [-37.45380985900664 +/- 7.39e-15]
b = [606.050644636309 +/- 2.44e-13]
c = [-5047.98924959459 +/- 3.10e-12]
```

Tail-start samples for `k=22`:

```text
T=25: q/T=9.000000000000000000E-1, |a|q/T=3.370842887310597575E+1, F_k=-3.705176510564669320E+3, status=invalid_truncation_normalizer, B/bound=n/a
T=50: q/T=4.500000000000000000E-1, |a|q/T=1.685421443655298788E+1, F_k=-4.108228689777087154E+2, status=invalid_truncation_normalizer, B/bound=n/a
T=100: q/T=2.250000000000000000E-1, |a|q/T=8.427107218276493938E+0, F_k=-4.077572511834853726E+1, status=invalid_truncation_normalizer, B/bound=n/a
T=200: q/T=1.125000000000000000E-1, |a|q/T=4.213553609138246969E+0, F_k=-3.376514867254146729E+0, status=invalid_truncation_normalizer, B/bound=n/a
T=500: q/T=4.500000000000000000E-2, |a|q/T=1.685421443655298788E+0, F_k=7.322730386875991354E-2, status=overbound_positive_truncation, B/bound=1.362086182073576929E+1
T=1000: q/T=2.250000000000000000E-2, |a|q/T=8.427107218276493938E-1, F_k=4.123450107885821700E-1, status=bounded_positive_truncation, B/bound=4.583846541015252075E-2
T=2000: q/T=1.125000000000000000E-2, |a|q/T=4.213553609138246969E-1, F_k=6.505827654571353051E-1, status=bounded_positive_truncation, B/bound=3.929354468000893553E-3
```

Simple first-correction stability budget:

```text
abs(a)*(k+1/2)/T <= 1/2
q/T <= 1.334977674853986504E-2
k=22 enters only after T >= 1.685421443655298788E+3
```

Interpretation: the low-order Taylor truncation is not a finite proof
model for the certified `lambda=-25,-50,-100` prefix. The actual prefix
is certified by ACB enclosures elsewhere; this artifact only sizes the
local Taylor/remainder theorem needed for an analytic route.

Required analytic upgrades:

```text
1. prove F_k^(<=6)+R_k^(>=8)>0 over a stated q/T range
2. bound the second log-difference of R_k^(>=8) against the B_k target margin
3. control the monotone third log-difference or use the separate monotone-contraction theorem
4. hand off to a global far-tail moving-saddle theorem or finite collar plus analytic tail
```

Integration:

```text
outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md
outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
```

Summary:

The exact Gaussian-moment budget makes q/T the relevant local parameter. The low-order Taylor truncation is not a finite proof model for the actual finite prefix: at the k=22 tail start it has invalid positive-moment normalizers for T=25,50,100,200, and the Taylor remainder must prove positivity plus log-curvature control before the signed perturbation route can enter the bounded-curvature target.
