# Jensen-Window PF Newman Theta/Bessel Higher-Shift Regularization Gate

Date: 2026-07-11

Status: exact higher-shift expansion with spectral summation obstruction.
This is not a proof or disproof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.py
```

Builds on `outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md`.

Current result:

```text
validated Jensen-window PF Newman theta/Bessel higher-shift regularization gate: 9 rows, 0 issues, 3 exact expansion identities, 1 coefficient sign theorem, 1 fixed-block Bessel theorem, 3 spectral non-promotion gates, 1 coupled modular handoff
```

## Exact Higher Shifts

Set

```text
phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u))
g_n(u)=(phi_n(u)+phi_n(-u))/2
Phi(u)=sum_(n>=1)g_n(u)
```

Symmetrizing before expanding gives, for every fixed `n>=1`,

```text
c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3*m)/m!
g_n(u)=exp(-2*pi*n^2*cosh(4u))*sum_(m>=0) c_(n,m)*cosh((9-4*m)*u)
```

The series and all kernel derivatives converge locally uniformly. Its
absolute `m`-majorant is integrable on the positive half-line for each
fixed `n`, so the individual block transform is ordinary, not formal.

## Sign Law

```text
sign(c_(n,m))=sign(2*pi^2*n^4-3*m)
c_(n,m)>0 iff m<2*pi^2*n^4/3
For n=1, c_(1,m)>0 for 0<=m<=6 and c_(1,m)<0 for m>=7; the bounds 3<pi<22/7 place 2*pi^2/3 strictly between 6 and 7.
```

Thus higher shifts really do enter, but with infinitely many signed
coefficients. The sign change is intrinsic and cannot be removed by
treating every shifted block as an independent positive residual.

## Fixed-Block Transform

The standard Bessel integral gives

```text
P_(n,m)(x)=1/8*(K_(9/4-m+i*x/4)(2*pi*n^2)+K_(9/4-m-i*x/4)(2*pi*n^2))
integral_0^infinity exp(-2*pi*n^2*cosh(4u))*cosh((9-4m)u)*cos(xu)du=P_(n,m)(x)
I_n(x):=integral_0^infinity g_n(u)cos(xu)du=sum_(m>=0)c_(n,m)*P_(n,m)(x)
```

This theorem is deliberately indexed by one fixed `n`. The arithmetic
sum cannot then be moved through the Fourier integral.

## Spectral Obstruction

The translate identity gives exactly

```text
Fourier[phi_n](x)=n^(-1/2-i*x/2)*Qhat(x)
I_n(0)=Qhat(0)/(2*sqrt(n)), Qhat(0)=-Gamma(1/4)/(32*pi^(1/4))<0
sum_(n>=1)I_n(0)=-infinity
```

But `H_0(0)=integral_0^infinity Phi(u)du` is finite. Therefore the
exact pointwise kernel expansion is not an ordinary termwise Bessel
decomposition of the Xi transform. This is a concrete failure of
Fubini/Tonelli, not a failure of the theta identity.

Analytic continuation does not repair the proof route:

```text
Qhat(x)*sum_(n>=1)n^(-1/2-i*x/2)
Qhat(x)*zeta((1+i*x)/2)=xi((1+i*x)/2)/4
```

It reconstructs Xi itself and supplies no positive block summation.

## Live Handoff

```text
M_(i,j)=F_i'*F_j'-(F_i*F_j''+F_j*F_i'')/2
L[sum_j a_j*F_j]=sum_(i,j)a_i*a_j*M_(i,j)
```

Group the theta summands through the modular identity before the spectral transform, then prove positivity of the resulting coupled mixed-term quadratic form. Alternatively supply a noncircular renormalized summation theorem that preserves the needed sign.

The theta cell-renormalization gate supplies an ordinary normally convergent coupled matrix at t=0, but proves that every individual cell block retains an exp(-5u) tail and cannot be deformed to t>0.
See `outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md`.

Naive ordinary higher-shift summation is now closed. A surviving proof
must retain the modular cancellation inside the coupled mixed terms
before taking the spectral transform.

Reference: https://dlmf.nist.gov/10.32.
