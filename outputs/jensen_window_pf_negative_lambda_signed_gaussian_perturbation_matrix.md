# Jensen-Window PF Negative-Lambda Signed Gaussian Perturbation Matrix

Date: 2026-07-06

Status: finite theorem-search diagnostic. This is not a proof of the
bounded log-curvature theorem, cone entry, Jensen-window PF-infinity,
RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix`.

Proof boundary: this artifact packages the fixed-k signed Gaussian
perturbation route. It does not prove a uniform asymptotic remainder
or an all-`k` tail theorem.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py
```

Current result:

```text
validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: 8 matrix rows, 2 certified Taylor signs, 1 fixed-k activation estimates, 0 ready-to-apply rows, 0 issues
```

## Perturbation Route

For:

```text
Phi(u)=c0*(1+a*u^2+b*u^4+c*u^6+...)
lambda=-T
```

the fixed-k expansion gives:

```text
log x_k = (2*b-a^2)/T^2 + O_k(T^-3)
B_k = -log(x_k) = (a^2-2*b)/T^2 + O_k(T^-3)
log(x_(k+1)/x_k) = 2*(a^3-3*a*b+3*c)/T^3 + O_k(T^-4)
```

Certified signs:

```text
a = [-37.45380985900664 +/- 7.39e-15]
b = [606.050644636309 +/- 2.44e-13]
c = [-5047.98924959459 +/- 3.10e-12]
2*b-a^2 = [-190.68658368200 +/- 5.53e-12]
D2=a^2-2*b = [190.68658368200 +/- 5.53e-12]
2*(a^3-3*a*b+3*c) = [825.997624927 +/- 3.55e-10]
```

Leading fixed-k activation estimate:

```text
K = 21
T >= [110.902139959017 +/- 7.87e-13]
```

This estimate ignores remainders and therefore cannot be used as a proof.
It is a scale diagnostic for the signed perturbation route only.

Integration:

```text
outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md
outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md
outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md
```

Summary:

The surviving Gaussian-baseline route is a signed or tilted perturbation, not a positive Gaussian mixture: the certified Phi Taylor signs give a positive leading deficit D2=a^2-2*b and positive monotone correction 2*(a^3-3*a*b+3*c), but the fixed-k expansion still needs uniform remainder control or a finite-collar tail theorem before it can close the bounded log-curvature target.
