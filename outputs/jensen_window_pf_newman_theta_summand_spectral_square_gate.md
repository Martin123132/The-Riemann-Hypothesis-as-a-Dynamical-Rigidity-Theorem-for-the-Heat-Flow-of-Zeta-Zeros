# Jensen-Window PF Newman Theta-Summand Spectral-Square Gate

Date: 2026-07-11

Status: exact theta differential/Mellin reduction with finite-summand
spectral obstruction. This is not a proof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_theta_summand_spectral_square_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_summand_spectral_square_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_summand_spectral_square_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_summand_spectral_square_gate.py
```

Current result:

```text
validated Jensen-window PF Newman theta-summand spectral-square gate: 12 rows, 0 issues, 7 exact transform identities, 1 xi-reconstruction non-promotion gate, 2 exact finite-truncation theorems, 1 numerical sign diagnostic, 1 infinite-cancellation handoff
```

## Theta Primitive

For `u>=0`,

```text
phi_n(u)=n^(-1/2)*phi_1(u+(log n)/2)
h(u)=exp(u-pi*exp(4u))
phi_1(u)=(h''(u)-h(u))/8
R(u)=sum_(n>=1) exp(u-pi*n^2*exp(4u)), Phi(u)=(R''(u)-R(u))/8 for u>=0
```

The theta functional equation gives `R'(0)=-1/2`. Two integrations
by parts therefore yield

```text
H_0(x)=1/16-(1+x^2)*C(x)/8
```

This formula exposes a cancellation, not a square: the constant cancels
the full algebraic high-frequency expansion of the second term.
The same integration by parts works at every target time:

```text
C_t(x)=integral_0^infinity exp(tu^2)*R(u)*cos(xu)du
D_t=-4*t^2*partial_x^2+4*t*x*partial_x+(2*t-1-x^2)
H_t(x)=1/16+D_t[C_t](x)/8
With A_t=D_t[C_t], L_t(x)=(A_t'(x)^2-(A_t(x)+1/2)*A_t''(x))/64.
```

## Mellin Audit

The bilateral profile transform is

```text
Qhat(z)=-(1+z^2)/32*pi^(-(1+i*z)/4)*Gamma((1+i*z)/4)
sum_(n>=1)n^(-1/2)*exp(-i*z*log(n)/2)=zeta(s)
Qhat(z)*zeta(s)=xi(s)/4
```

Thus the most direct arithmetic spectral assembly reconstructs the
completed zeta function itself. It does not factor xi into a new
zero-free term, and using its zeros here would be circular.

## Finite-Sum Theorem

For the first `N` summands,

```text
phi_n'(0)=pi*n^2*exp(-pi*n^2)*(-8*pi^2*n^4+30*pi*n^2-15)
A_N=f_(N,t)'(0)=sum_(n=1)^N phi_n'(0)=-sum_(n>N)phi_n'(0)>0, independent of t
H_(N,t)(x)=-A_N/x^2+O_(N,t)(x^-4), H_(N,t)'(x)=2A_N/x^3+O(x^-5), H_(N,t)''(x)=-6A_N/x^4+O(x^-6)
L[H_(N,t)](x)=-2*A_N^2/x^6+O_(N,t)(x^-8)<0 for all sufficiently large x
```

So every finite truncation fails the global first-Laguerre criterion.
The obstruction moves to larger frequency as `N` grows because `A_N`
tends to zero, but it never disappears at finite `N`.

## Pairwise Guard

L[sum_n H_n]=sum_n L[H_n]+sum_(m<n)B[H_m,H_n]
The hoped-for entrywise sign is false. Independent 60-digit quadrature
gives:

```text
t=0, x=10.5, cross_1_2=-9.47240435132102504e-11, full L=1.40662152605622681e-05
t=0, x=16.75, cross_1_2=-5.45580704229529059e-09, full L=2.24311139002416139e-06
t=0, x=39.5, self_1=-1.43539553077037718e-11, full L=8.11247906394412734e-12
t=0.5, x=10.5, cross_1_2=-2.00989272997002332e-10, full L=1.42347941351116568e-05
t=0.5, x=16.75, cross_1_2=-5.57165768004164914e-09, full L=2.23254898230694973e-06
t=0.5, x=39.5, self_1=-1.49130918219351315e-11, full L=7.78250223341584428e-12
```

The full sampled Laguerre expression is positive precisely where its
first self or cross component is negative. Any proof must therefore
use global arithmetic cancellation rather than termwise positivity.

## Live Handoff

For A_t=D_t[C_t], prove A_t'(x)^2-(A_t(x)+1/2)*A_t''(x)>0 for every real x and 0<t<=1/2 using the positive theta primitive R and its modular identity. The formula already performs the infinite odd-endpoint cancellation; the missing step is a noncircular curvature estimate. Finite summand squares and pairwise-positive cross terms are closed routes.

The endpoint-subtracted transform is now exact. The remaining task is
the displayed curvature inequality, uniformly for `0<t<=1/2`; finite
theta blocks are a closed route.

References: https://dlmf.nist.gov/20.7 and https://dlmf.nist.gov/25.4.
