# Jensen-Window PF Newman Polymath-15 Critical Wronskian Phase Reduction

Date: 2026-07-17

Status: exact crossing-localized reduction of the corrected `C1`
target. The arithmetic separation remains open; this is not a proof
of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.py
```

## Complex Main

Absorb one half of the explicit endpoint lift into the complex finite
sum:

```text
Choose the sharp endpoint half q_(N,t) so Q_(N,t)=2*Re(q_(N,t)) on the real axis, and set E_(N,t)=exp(i*beta_t)D_(N,t)-q_(N,t); then J_(N,t)=2*Re(E_(N,t)) and J_(N,t)'=2*Re(E_(N,t)')
For E=X+iY and E'=U+iV, T_L[J]=4*(X^2+(U/L)^2)
```

For `E=a*exp(i*theta)`, the exact polar form is

```text
If E=a*exp(i*theta), u=(log(a))', and v=theta', then T_L[J]=4*a^2*(cos(theta)^2+(u*cos(theta)-v*sin(theta))^2/L^2) and W_E=a^2*v
At Re(E)=0 with E!=0, T_L[J]=4*W_E^2/(|E|^2*L^2); J=J'=0 if and only if theta'=0
```

Thus a corrected real-part crossing is multiple precisely when its
complex phase is stationary, apart from the separately visible case
`E=0`.

## Robust Collision Gate

The certified remainder caps imply

```text
If Z=J+r=Z'=0, |r|<=epsilon_0, and |r'|<=L*epsilon_1, then 2*|W_E|<=epsilon_0*|E'|+L*epsilon_1*|E|
```

Consequently it is enough to establish the disjunction

```text
For every critical point, either |Re(E)|>epsilon_0/2 or 2*|W_E|>epsilon_0*|E'|+L*epsilon_1*|E|
With epsilon_0=2500*exp(-3L/4) and epsilon_1=5000*exp(-3L/4), either |Re(E)|>1250*exp(-3L/4) or 2*|W_E|>exp(-3L/4)*(2500*|E'|+5000*L*|E|)
```

The target formulas remain valid on their original broader domain,
while the current proof partition is

```text
The explicit formulas remain valid on the original domain L>=50, 0<tL<=25. After the oscillatory-zeta theorem, the live asymptotic high-frequency layer is 0<tL<=c_*+o(1), with c_*=4911678521/1933561194; uncovered bounded L remains a separate compact-certificate obligation
```

This localizes the arithmetic theorem to an exponentially thin collar
around corrected real-part crossings. It does not require a global sign
for the phase derivative.

## Sign Shortcut Rejected

At moderate frequency and `t=0`, direct corrected evaluations give
both signs of `theta'`. The same occurs at two certified-to-high-precision
corrected crossings:

- `x=517.2201245916915343451981788255318076677`: `theta'=-0.55142810044966476597942521956113449`, residual `2.635758643014824499667915e-74`
- `x=519.7486269586060360224240422895051497992`: `theta'=1.2989222762566859009109017488767268`, residual `1.385241436237148917911971e-74`

These finite diagnostics only shape the theorem search. They show that
one-sided phase monotonicity is not supplied by the algebra alone; they
do not prove or disprove the `L>=50` thin-collar separation.

## Dynamical Bridge

At a corrected crossing the phase Wronskian is also the exact
denominator of the zero-flow generator:

```text
At J=0 with E=iY and Y!=0, W_E=-Y*J'/2; along a moving simple zero x_*(t), x_*'=-J_t/J'=J_t*Y/(2*W_E)
If J(x)=(x-alpha)*(x-beta)*G(x), then |W_E(alpha)|=|E(alpha)|*(beta-alpha)*|G(alpha)|/2; the Wronskian margin is linearly close-pair sensitive
```

The corrected Lehmer-pair diagnostic has gap
`0.074957277806787566257949131543400340` and minimum absolute
phase velocity
`0.156356967747187753117412881978`.
This is the finite-main signature of the same square-root velocity
blow-up identified by the original dynamical-rigidity programme.

## Live Target

```text
Where E!=0, collision means theta is pi/2 modulo pi and theta'=0; the open arithmetic theorem is quantitative avoidance of those phase critical levels
```

The problem is now a quantitative phase-critical-value avoidance theorem
for the corrected Riemann-Siegel complex main. That theorem remains
RH-level and is not established here.
