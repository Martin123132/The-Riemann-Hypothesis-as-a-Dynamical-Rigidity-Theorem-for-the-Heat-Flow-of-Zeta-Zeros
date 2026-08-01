# Newman Theta Fixed-Block Cofinal Obstruction Gate

Date: 2026-07-24

Status: exact method obstruction with a live cancellation-aware handoff.
This is not a proof or disproof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fixed_block_cofinal_obstruction_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fixed_block_cofinal_obstruction_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fixed_block_cofinal_obstruction_gate.py
```

## The Fixed-Block Setup

For a fixed finite `N`, put

```text
phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u))
f_(N,t)(u)=exp(tu^2)*sum_(1<=n<=N)phi_n(u)
H_(N,t)(x)=integral_0^infinity f_(N,t)(u)*cos(xu)du
```

For fixed N, every u-derivative of f_(N,t) is uniformly integrable for 0<=t<=1/5, and every integration-by-parts boundary term at infinity vanishes.

The double-exponential theta tail makes all derivatives used below uniformly integrable for `0<=t<=1/5`.

## Two Integrations By Parts

Define the finite uniform constants

```text
A_N=sup_(0<=t<=1/5)(|f_(N,t)'(0)|+||f_(N,t)''||_1)<infinity; B_N=sup_(0<=t<=1/5)||(u*f_(N,t))''||_1<infinity
```

Then, for every `x>0`,

```text
H_(N,t)(x)=-f_(N,t)'(0)/x^2-x^(-2)*integral_0^infinity f_(N,t)''(u)*cos(xu)du
|H_(N,t)(x)|<=A_N/x^2
partial_x H_(N,t)(x)=x^(-2)*integral_0^infinity (u*f_(N,t)(u))''*sin(xu)du
|partial_x H_(N,t)(x)|<=B_N/x^2
```

The derivative identity uses `g(u)=u*f_(N,t)(u)`, for which `g(0)=0`.

## Exact Cofinal Obstruction

The existing direct certificate uses static bars

```text
E_0(x)=16*x^4*epsilon_0; E_1(x)=64*x^3*epsilon_0+16*x^4*epsilon_1, where epsilon_0>0 and epsilon_1>=0 are fixed
```

Since `J_(N,t)=16*x^4*H_(N,t)`, normalization cancels the apparent `x^4` advantage:

```text
|J_(N,t)(x)|/E_0(x)=|H_(N,t)(x)|/epsilon_0<=A_N/(epsilon_0*x^2)
|J_(N,t)'(x)|/E_1(x)=|4H_(N,t)(x)+x*partial_x H_(N,t)(x)|/(4epsilon_0+xepsilon_1)<=(4A_N/x^2+B_N/x)/(4epsilon_0+xepsilon_1)
```

Take

```text
X_N=max(1,sqrt(2A_N/epsilon_0),B_N/(2epsilon_0))
```

For `x>X_N`, the value bound is below `1/2`. Also `4A_N/x^2<2epsilon_0` and `B_N/x<2epsilon_0`, while the derivative denominator is at least `4epsilon_0`. Therefore

```text
For every x>X_N and every 0<=t<=1/5, |J_(N,t)(x)|<E_0(x) and |J_(N,t)'(x)|<E_1(x)
```

Every sufficiently large diagonal shell

```text
S_j=[1/(5j),1/5]x[38,38+j]
```

contains such an `x`. Hence a fixed number of retained theta blocks plus fixed positive absolute-moment bars cannot close the cofinal family.

## What The Bars Erase

At the half-line endpoint,

```text
phi_n'(0)=pi*n^2*exp(-pi*n^2)*(-8*pi^2*n^4+30*pi*n^2-15)
a_2=f_(2,t)'(0)=sum_(n=1)^2 phi_n'(0)>0
a_2 approximately 0.000000082652795777273394114371932109122889696862795688896760119283063044133246939086341
For n>=3, 4*pi*n^2-15>0, hence -8*pi^2*n^4+30*pi*n^2-15<0. The differentiated theta series converges locally uniformly. Since Phi is even, sum_(n>=1)phi_n'(0)=0, so a_2=-sum_(n>=3)phi_n'(0)>0.
Phi is even, so (exp(tu^2)Phi(u))'(0)=0 and sum_(n>=3)phi_n'(0)=-a_2<0
```

The omitted theta tail cancels the algebraic x^(-2) endpoint term of the retained two-block transform. Static absolute moments discard this cancellation.

This is the structural point: the full even Xi kernel has no odd endpoint jet, but every fixed arithmetic truncation generally does. The omitted tail cancels that jet. Taking its absolute mass first destroys the cancellation before the oscillatory transform is used.

## Finite Scout

The separate resumable scout records

```text
shells=5..30
certified rows=26
evaluated boxes=20262
subdivisions=1292
smallest stored ratio occurs at j=23: [1.3709314869015039163669568274457936340118136691538 +/- 2.92e-50]
```

Those are rigorous stored finite diagnostics, but they are not promoted here as new `Q_j` theorems. Passing through `j=30` and failing cofinally are compatible: the obstruction is asymptotic.

## Live Handoff

The exact repeated endpoint-jet formula is

```text
integral_0^infinity f(u)cos(xu)du=sum_(r=0)^(m-1)(-1)^(r+1)*f^(2r+1)(0)/x^(2r+2)+(-1)^m*x^(-2m)*integral_0^infinity f^(2m)(u)cos(xu)du
```

A cofinal theta-shell proof must let N grow with x, replace static moment bars by oscillatory remainder bounds that retain endpoint-jet cancellation, or use a modular grouping that preserves the cancellation before transformation.

Open target:

```text
Construct a uniform cancellation-aware first-jet remainder theorem strong enough to close every S_j without assuming zero simplicity or RH.
```

## Proof Boundary

This is an obstruction to the fixed-N static absolute-tail proof method, not evidence of a common zero of the full Xi first jet.
No common zero, new `Q_j`, `Lambda<=0`, RH, or Clay-prize conclusion follows.

Current validation:

```text
validated Newman theta fixed-block cofinal obstruction gate: 9 rows, 0 issues, 3 exact integration-by-parts bounds, 1 cofinal method obstruction, 1 endpoint-cancellation theorem, 26 finite diagnostic shells, 1 cancellation-aware open handoff
```
