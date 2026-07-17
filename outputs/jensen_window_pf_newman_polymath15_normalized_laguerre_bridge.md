# Jensen-Window PF Newman Polymath-15 Normalized Laguerre Bridge

Date: 2026-07-16

Status: exact normalized-curvature transfer from a published
positive-time Riemann-Siegel theorem. This is not a proof of
`Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.py
```

Published input: D. H. J. Polymath, [Effective approximation of heat
flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438), Theorem 1.3.

## Real Normalization

For `s=(1-i*x)/2`, put

```text
B_t(x+i*y)=M_t((1+y-i*x)/2)
M_t(s)=exp((t/4)*alpha(s)^2)*M_0(s)
alpha(s)=1/(2s)+1/(s-1)+(1/2)Log(s/(2*pi))
N=floor(sqrt(x/(4*pi)+t/16))
D_(N,t)(x)=sum_(n<=N) exp((t/4)*log(n)^2)*n^(-s_*)
```

On the real axis the published two-sum approximation has `kappa=0`
and its second sum is the phase-adjusted conjugate of the first:

```text
f_t(x)=D_(N,t)(x)+(conj(B_t(x))/B_t(x))*conj(D_(N,t)(x))
P_(N,t)(x)=2*Re(exp(i*beta_t(x))*D_(N,t)(x)), B_t(x)=A_t(x)*exp(i*beta_t(x))
```

Writing `Z_t=H_t/A_t`, exact product differentiation gives

```text
L[H_t]=A_t^2*C_t[Z_t], Z_t=H_t/A_t, C_t[Z]=Z'^2-Z*Z''+V_t*Z^2, V_t=-(log A_t)''
```

Thus the complex approximation becomes a real finite curvature
problem with an explicit normalizer potential.

## Jet Transfer

If one fixed-`N` approximation is holomorphic on a disk of radius
`rho` and its remainder is bounded there by `E_rho`, Cauchy's
estimate and the phase product give

```text
eps_0=E_rho
eps_1=(1/rho+|beta'|)*E_rho
eps_2=(2/rho^2+2*|beta'|/rho+|beta''|+beta'^2)*E_rho
Err=2*|P'|*eps_1+eps_1^2+|P|*eps_2+eps_0*|P''|+eps_0*eps_2+|V|*(2*|P|*eps_0+eps_0^2); C[P]>Err implies L[H_t]>0
```

This is the missing derivative bridge: a scalar Riemann-Siegel
remainder can certify the Laguerre expression once the finite main
sum clears the explicit error.

## Phase Floor

The `n=1` term is exactly `P_0=2*cos(beta_t)`, with

```text
P_0=2*cos(beta); C[P_0]=4*beta'^2+2*beta''*sin(2*beta)+4*V*cos(beta)^2
```

| t | x | beta' | beta'' | V | phase-only lower |
|---:|---:|---:|---:|---:|---:|
| 0 | 200 | -0.691813905 | -0.0012500937 | 4.3747031e-5 | 1.91192573 |
| 0 | 1000 | -1.09418238 | -0.00025000075 | 1.7499953e-6 | 4.78844035 |
| 0 | 10000 | -1.66982903 | -2.5000001e-5 | 1.75e-8 | 11.1532659 |
| 0 | 1000000 | -2.82112158 | -2.5e-7 | 1.75e-12 | 31.8349073 |
| 0.01 | 200 | -0.69182411 | -0.0012500419 | 4.3798952e-5 | 1.91198231 |
| 0.01 | 1000 | -1.09418438 | -0.00024999874 | 1.7541889e-6 | 4.7884578 |
| 0.01 | 10000 | -1.66982922 | -2.4999981e-5 | 1.7570964e-8 | 11.1532686 |
| 0.01 | 1000000 | -2.82112158 | -2.5e-7 | 1.7628556e-12 | 31.8349074 |
| 0.1 | 200 | -0.691915953 | -0.001249575 | 4.4266238e-5 | 1.91249159 |
| 0.1 | 1000 | -1.09420231 | -0.00024998061 | 1.7919317e-6 | 4.78861485 |
| 0.1 | 10000 | -1.669831 | -2.4999803e-5 | 1.8209639e-8 | 11.1532922 |
| 0.1 | 1000000 | -2.8211216 | -2.4999998e-7 | 1.8785558e-12 | 31.8349078 |
| 0.5 | 200 | -0.692324143 | -0.0012475003 | 4.6343065e-5 | 1.91475587 |
| 0.5 | 1000 | -1.09428204 | -0.00024990006 | 1.9596777e-6 | 4.78931289 |
| 0.5 | 10000 | -1.66983887 | -2.4999014e-5 | 2.1048199e-8 | 11.1533974 |
| 0.5 | 1000000 | -2.82112168 | -2.499999e-7 | 2.392779e-12 | 31.8349095 |

These 16 rows are diagnostics, not a ray proof. They show that the
normalizer itself supplies a large positive phase-speed floor; the
remaining issue is controlling the arithmetic terms and the theorem
remainder in `C^2`.

## Remaining Geometry

The published cutoff exactly matches the square-root scale found by
the modular-blend saddle audit. Cauchy disks cannot cross its discrete
transition points without an overlapping adjacent-`N` estimate:

```text
For fixed N, 4*pi*(N^2-t/16)<=x<4*pi*((N+1)^2-t/16); a Cauchy disk must remain in one cell.
```

For the nonconstant arithmetic terms, the published lower bound on
`Re(s_*)` yields

```text
|b_n^t*n^(-s_*)|<=n^(-1/2)*exp(-(t/4)*log(n)*log((x/(4*pi))/n)+(t/(2*x^2))*log(n)); hence n>=2 dominance is controlled by t*log(x), and the nonuniform regime is t*log(x)=O(1).
```

So the global tail now has two sharply separated tasks: prove
single-saddle dominance when `t*log(x)` is large, and handle the
shrinking but RH-critical boundary layer `t*log(x)=O(1)` without
discarding the coupled Dirichlet cancellation.
