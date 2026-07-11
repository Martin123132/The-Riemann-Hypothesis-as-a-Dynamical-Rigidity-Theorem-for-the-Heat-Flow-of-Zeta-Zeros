# Jensen-Window PF Newman Classical-Field Balance Gate

Date: 2026-07-11

Status: exact classical-field balance reduction with fixed-time
compactness gate. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_classical_field_balance_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_classical_field_balance_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_classical_field_balance_gate.py
```

Current result:

```text
validated Jensen-window PF Newman classical-field balance gate: 10 rows, 0 issues, 3 exact field identities, 1 arithmetic equilibrium, 1 continuum -pi/8 benchmark, 1 quantile-drift match, 1 published fixed-time localization theorem, 2 exact sensitivity countermodels, 1 compactness reduction, 1 open lambda-uniform handoff
```

## Exact Balance Identity

At a double zero with the colliding pair removed,

```text
For H(z)=(z-c)^2*V(z), B_c=V'(c)/V(c)=PV sum_(signed roots y outside the pair) 1/(c-y).
If c=c0+delta_c and y_j=r_j+delta_j, then 1/(c-y_j)-1/(c0-r_j)=(delta_j-delta_c)/((c-y_j)(c0-r_j)).
|B(c;y)-B(c0;r)| <= sum_j |delta_j-delta_c|/(|c-y_j|*|c0-r_j|), provided the paired sums converge.
```

The constant-density reference is exact arithmetic equilibrium:

```text
For roots c+(j+1/2)h, j in Z, removing the pair c+/-h/2 gives PV sum 1/(c-y)=0 by termwise reflection.
```

## Classical Density Benchmark

For the Riemann-von Mangoldt density,

```text
rho_0(y)=Psi'(y)=log(y/(4*pi))/(4*pi), y>=a>4*pi
B_rho(x)=2*x*PV integral_a^infinity rho_0(y)/(x^2-y^2)dy
B_rho(x)=-pi/8+O_a(1/x)
```

The two exact principal-value integrals behind the constant are

```text
PV integral_0^infinity du/(1-u^2)=0
PV integral_0^infinity log(u)/(1-u^2)du=-pi^2/4
```

Adding the positive-time density correction leaves the same limit:

```text
time*log(sqrt(-a**2 + x**2)/a)/(8*x)
B_(rho,t)(x)=-pi/8+O_a((1+t*log(x/a))/x)
```

Implicit differentiation of the published quantile law gives

```text
-pi*level*(-log(level) + log(pi) + 2*log(2))/(-4*level*log(level) + 4*level*log(pi) + 8*level*log(2) - pi*time)
dx_n/dt=-pi/4+o(1), so B_n=(1/2)dx_n/dt -> -pi/8
```

The classical reference field is therefore -pi/8 at leading order.

This matches the high-positive-time zero drift described in the
Polymath asymptotics.

## Fixed-Time Compactness

For each fixed 0<t<=1/2, sufficiently high zeros are real, simple, and uniquely localized near the quantiles g(x_n,t)=n once x_n>=exp(C/t).

Consequently, if `Lambda>0`, every multiple zero at the boundary
is confined to a region of scale `exp(C/Lambda)`. The threshold
The threshold diverges as `Lambda->0+`, so this is a compactness reduction rather
than the desired exclusion theorem.

Primary sources: https://arxiv.org/abs/1801.05914 and
https://arxiv.org/abs/1904.12438.

## Macroscopic Guard

Signed integer lattice with double roots at +/-n; after removing the +n pair its field is 1/(2*n).

```text
+/-(n+1) to +/-(n+eps): B=(4*eps**2*n**2 + 2*eps**2*n + eps**2 + 8*eps*n**3 + 4*eps*n**2 + 2*eps*n - 8*n**3 - 4*n**2)/(2*eps*n*(eps + 2*n)*(2*n + 1)) -> -infinity
+/-(n-1) to +/-(n-eps): B=-(-4*eps**2*n**2 + 2*eps**2*n - eps**2 + 8*eps*n**3 - 4*eps*n**2 + 2*eps*n - 8*n**3 + 4*n**2)/(2*eps*n*(-eps + 2*n)*(2*n - 1)) -> +infinity
```

Both deformations preserve even symmetry and move every zero by less than one.
Therefore bounded classical-location error, and a fortiori
macroscopic counting control, cannot determine the collision field.

## Live Handoff

Prove a lambda-uniform paired reciprocal-gap discrepancy bound strong enough for the weighted perturbation sum, together with an explicit far-tail bound, or prove a different analytic exclusion that remains effective as lambda->0+.

This is the next strict gate: the estimate must remain effective as
`lambda->0+`; fixed-lambda asymptotics do not prove `Lambda<=0`.
