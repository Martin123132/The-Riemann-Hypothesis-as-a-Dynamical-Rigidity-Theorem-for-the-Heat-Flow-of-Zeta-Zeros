# Jensen-Window PF Negative-Lambda Relative-Gaussian Full-Kernel Evenness/Cauchy Lemma

Date: 2026-07-10

Status: exact full-kernel evenness and Cauchy lemma. This is not a proof
of a Gamma-expectation bound, cone entry, RH, or `Lambda <= 0`.

Machine-readable artifact:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.json
```

Generator and checker:

```text
work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.py
work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.py
```

## Exact Derivation

Set

```text
theta(x)=sum_{n in Z} exp(-pi*n^2*x)
Psi(x)=(theta(x)-1)/2=sum_{n>=1}exp(-pi*n^2*x)
L_x = 2*x^2*d_x^2 + 3*x*d_x.
```

The Jacobi transformation gives

```text
Psi(x)=x^(-1/2)*Psi(1/x)+(x^(-1/2)-1)/2
```

while direct differentiation gives

```text
Phi(u)=exp(u)*(L_x Psi)(x), x=exp(4u), L_x=2*x^2*d_x^2+3*x*d_x
L_x[x^(-1/2)f(1/x)]=x^(-1/2)*(L_y f)(1/x)
L_x[1]=L_x[x^(-1/2)]=0.
```

Therefore `Phi(-u)=Phi(u)` exactly. This is an infinite-kernel theorem;
it is not inferred from the tiny odd coefficients of `Phi_30`.

## Order-42 Residual

On every closed disk `|z|<=R<pi/8`, the full defining series and its
derivatives converge uniformly because

```text
Re(exp(4z)) >= exp(-4R)*cos(4R) > 0.
```

At `R=0.38`, the recorded lower bound is `0.011104983408849904885559698512335034910163134523327`.

Exact evenness makes every odd Taylor coefficient vanish. Hence

```text
G_20(z)=Phi(z)/Phi(0)-sum_{j=0}^{20}r_j*z^(2j)
G_20^(m)(0)=0 for 0<=m<=41
```

so `G_20` has a zero of order at least `42`.

## Cauchy Factors

If sup_|z|<=R |Phi(z)|<=M and Phi(0)>=c0, then for 0<=x<R, |G_20(x)|<=M*x^42/(c0*R^42*(1-x/R)).

|x*G_20'(x)/2|<=M*x^42/(2*c0*R^42)*(42/(1-q)+q/(1-q)^2), q=x/R.

## Proof Boundary

Exact full-kernel parity, local analyticity, order-42 residual zero, and factored Cauchy bounds only. It does not bound a Gamma expectation by itself, does not prove cone entry or sign regularity, and does not prove RH or Lambda <= 0.

## Reproduction

```powershell
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.py
```

Current result:

```text
validated Jensen-window PF negative-lambda relative-Gaussian full-kernel evenness/Cauchy lemma: 9 rows, 0 issues, 3 symbolic identities, order>=42 residual zero, 0 ready-to-apply rows
```
