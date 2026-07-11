# Jensen-Window PF Newman Correlation Hierarchy Gaussian-Mixture Gate

Date: 2026-07-11

Status: exact correlation hierarchy with Gaussian-mixture and direct-PF
obstructions. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.py
```

Current result:

```text
validated Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate: 11 rows, 0 issues, 3 exact hierarchy identities, 1 universal boundary-contact signature, 1 Gaussian-mixture sufficient theorem, 2 numerical diagnostics, 1 exact super-Gaussian tail theorem, 2 non-promotion gates, 1 tail-compatible handoff
```

## Exact Hierarchy

Set

```text
K_(n,t)(v)=integral_R phi_t(s+v)*phi_t(s-v)*s^(2n) ds
F_(n,t)(xi)=Fourier[K_(n,t)](xi)
partial_t K_(n,t)(v)=2*v^2*K_(n,t)(v)+2*K_(n+1,t)(v)
partial_t F_(n,t)(xi)=-2*partial_xi^2 F_(n,t)(xi)+2*F_(n+1,t)(xi)
L_(n,t)(x)=2^(2n-1)/(2n)!*F_(n,t)(2x)
```

The hierarchy also resolves the exact contact made by a hypothetical
positive Newman boundary:

```text
If H_t has a real root c of multiplicity m, then F_(n,t)(2c)=0 for n<m and F_(m,t)(2c)=(2m)!/2^(2m-1)*(H_t^(m)(c)/m!)^2>0.
F_(2,t)(2c)=3*partial_xi^2 F_(1,t)(2c), partial_t F_(1,t)(2c)=H_t''(c)^2
```

For a double root this is a nonnegative quadratic touch of `F_1` with
a strictly positive `F_2` source. The hierarchy is therefore compatible
with birth; it becomes useful only with an additional coercive estimate.

## Gaussian-Mixture Test

If g_t is nonzero and completely monotone on [0,infinity), then g_t(r)=integral_[0,infinity) exp(-a*r)dmu_t(a).
Then

```text
Fourier[K_(1,t)](xi)=sqrt(pi)*integral_(a>0) a^(-1/2)*exp(-xi^2/(4a))dmu_t(a)>0
```

so this single property would close the strict Laguerre/Wiener target.
It survives a misleadingly strong local scout:

```text
t=0: alternating derivatives through n=8; g0*g2-g1^2=-7.86301502267677636e-06
t=0.10000000000000001: alternating derivatives through n=8; g0*g2-g1^2=-7.87135706726551851e-06
t=0.25: alternating derivatives through n=8; g0*g2-g1^2=-7.88311869073625437e-06
t=0.5: alternating derivatives through n=8; g0*g2-g1^2=-7.90063124955146627e-06
```

The determinant is negative, not positive: complete monotonicity
requires log-convexity. Independent 55-digit mpmath quadrature gives
`-7.8630150226764958e-6` at `t=0` and
`-7.9006312495514071e-6` at `t=1/2`.

## Exact Tail Obstruction

The Phi series gives

```text
0<Phi(u)<=C_Phi*exp(9u-pi*exp(4u)) for u>=0
0<K_(1,t)(v)<=C_Phi^2*M*exp(2*t*v^2+18v-pi*exp(4v)) for v>=0 and 0<=t<=1/2
lim_(r->infinity) -log(g_t(r))/r=+infinity, uniformly for 0<=t<=1/2
```

A nonzero completely monotone function is a positive Laplace mixture
and therefore has an exponential lower bound along the positive axis.
The displayed limit is impossible for such a mixture. Hence the
Gaussian-mixture route is rejected exactly for every `0<=t<=1/2`;
the numerical determinant merely shows that the failure is visible
already at the first moment-minor test.

## Closed Generic Routes

A differentiable even kernel has K'(0)=0. If it were convex on [0,infinity), K' would be nondecreasing and hence nonnegative; a nonconstant positive kernel tending to zero cannot satisfy the classical decreasing-convex Polya criterion.

A nondegenerate even PF-infinity function has a Gaussian lower tail when the Gaussian factor is present, and an exponential lower tail otherwise. It cannot decay faster than every Gaussian on both half-lines.
Therefore `K_(1,t)` itself cannot be PF-infinity. This sharpens the
previous handoff: an eventual total-positive argument must be a
different structured factorization, not direct PF membership.

Primary sources: https://doi.org/10.2307/1968466,
https://doi.org/10.1073/pnas.33.1.11, and
https://arxiv.org/abs/2006.16213.

## Live Handoff

Seek an Xi/theta-summand spectral-square decomposition or a relative correlation-hierarchy coercivity estimate that excludes the universal contact F_1=F_1'=0, F_2=3*F_1''>0. Gaussian mixtures, the smooth Polya convex criterion, and direct PF-infinity membership of K_(1,t) are closed.

The most concrete next test is a theta-summand spectral-square
decomposition. Failing that, the hierarchy needs a relative bound
strong enough to make `F_2=3*F_1''` impossible at a zero of `F_1`.
