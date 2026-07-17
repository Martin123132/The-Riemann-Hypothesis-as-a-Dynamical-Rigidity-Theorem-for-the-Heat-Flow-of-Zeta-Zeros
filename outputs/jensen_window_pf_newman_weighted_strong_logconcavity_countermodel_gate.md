# Jensen-Window PF Newman Weighted Strong-Log-Concavity Countermodel Gate

Date: 2026-07-17

Status: explicit root-concave theta-tail admissible weighted-correlation countermodel to generic shape-and-tail promotion. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.py
```

Current result:

```text
validated Jensen-window PF Newman weighted strong-log-concavity countermodel gate: 9 rows, 0 issues, 1 exact strong-curvature bound, 1 exact root-variable concavity theorem, 1 explicit theta-tail admissible kernel, 1 Gaussian endpoint identity, 1 exact endpoint witness, 1 Arb theta-tail witness, 1 weighted-correlation countermodel, 1 Xi-specific handoff
```

## Countermodel Family

```text
phi_(delta,epsilon)(u)=exp(-u^2-delta*u^4-epsilon*(cosh(4*u)-1))*(1+u^4/10), delta,epsilon>=0
```

For `z=u^2/sqrt(10)>=0`,

```text
z*(3-z^2)/(1+z^2)^2<=3*z/(1+z^2)<=3/2
(log phi_(delta,epsilon))''<=-(2-6/sqrt(10))-12*delta*u^2-16*epsilon*cosh(4*u)
2-6/sqrt(10)>0
```

Thus every `delta>=0` member is uniformly strongly log-concave. For every
`delta>0`, quartic damping makes the kernel admissible and
super-Gaussian; evenness and the curvature bound make it strictly decreasing
on the positive half-line.

At the explicit values `delta=1/10`, `epsilon=1/1000`, it also has the
published Xi kernel's stronger root-variable shape and tail properties:

```text
log(phi_(1/10,1/1000)(sqrt(r)))=-r-r^2/10+log(1+r^2/10)-(cosh(4*sqrt(r))-1)/1000
base root curvature=-r**2*(r**2 + 30)/(5*(r**2 + 10)**2)<0
d^2/dr^2 cosh(4*sqrt(r))=-(-4*r**(3/2)*cosh(4*sqrt(r)) + r*sinh(4*sqrt(r)))/r**(5/2)>0
phi_(1/10,1/1000)(u)<=exp(1/1000)*(1+u^4/10)*exp(-u^2-u^4/10-exp(4*abs(u))/2000)
```

## Exact Endpoint Witness

At `delta=0`, Gaussian differentiation gives

```text
sqrt(pi)*(x**4 - 12*x**2 + 172)*exp(-x**2/4)/160
L_1[F_0](x)=pi*(x**8 - 16*x**6 + 440*x**4 - 7680*x**2 + 37840)*exp(-x**2/2)/51200
L_1[F_0](3)=-743*pi*exp(-9/2)/51200<0
```

This endpoint identity is a cross-check. The admissible root-concave member is
certified directly at `delta=1/10`, `epsilon=1/1000`, `x=21/5`. Acb
integration on `[0,6]` and upper incomplete-gamma tail balls give

```text
F=[0.010079077726118413514324903055736354411454698489725216 +/- 3.24E-55]
F'=[-0.02883692696821660760997548984142684809891905250938752 +/- 3.43E-54]
F''=[0.1036847681860763492844176464243140394504271561765499 +/- 3.77E-53]
L_1[F]=[-0.000213478480591774964507646766853475076598648832979809 +/- 5.09E-55]<0
```

## Weighted-Correlation Failure

For the first Csordas correlation,

```text
K_(1,delta,epsilon)(v)=integral_R phi_(delta,epsilon)(s+v)*phi_(delta,epsilon)(s-v)*s^2 ds
L_1[F_(delta,epsilon)](x)=4*Fourier[K_(1,delta,epsilon)](2*x)
Fourier[K_(1,1/10,1/1000)](42/5)=L_1[F_(1/10,1/1000)](21/5)/4<0
```

Hence a smooth positive even admissible kernel with uniform strong
log-concavity, strict root-variable log-concavity, and the Xi theta tail class
can have a first `s^2`-weighted correlation that is not positive definite.
This is stronger than a zeroth-correlation shape guard: even the stronger
Xi-style shape and tail properties plus the generic weighted construction do
not supply the missing sign.

Primary source for the correlation identity:
https://arxiv.org/abs/1309.0055.

## Nonpromotion Gate

Even strict concavity of log(phi(sqrt(r))), uniform strong log-concavity, theta-type double-exponential decay, and admissibility of phi do not imply positive definiteness, let alone strict Fourier positivity, of K_(1). The Xi route must use arithmetic or modular structure beyond these generic shape and tail properties and the s^2-weighted correlation construction.

## Live Handoff

For Xi, seek a relative F_2/F_1 coercivity inequality or an s^2-weighted modular square that is absent from the countermodel.
