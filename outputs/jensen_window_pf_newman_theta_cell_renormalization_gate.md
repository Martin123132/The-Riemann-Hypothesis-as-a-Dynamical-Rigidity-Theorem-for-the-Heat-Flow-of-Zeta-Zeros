# Jensen-Window PF Newman Theta Cell-Renormalization Gate

Date: 2026-07-11

Status: exact endpoint cell renormalization with positive-time obstruction.
This is not a proof or disproof of RH or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_newman_theta_cell_renormalization_gate`.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_cell_renormalization_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_cell_renormalization_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_cell_renormalization_gate.py
```

Current result:

```text
validated Jensen-window PF Newman theta cell-renormalization gate: 10 rows, 0 issues, 4 exact kernel/transform identities, 3 exact convergence/sign theorems, 1 coupled Laguerre identity, 1 positive-time obstruction, 1 modular handoff
```

## Vanishing Continuum

Extend the arithmetic summand to every real `a>0`:

```text
phi_a(u)=(2*pi^2*a^4*exp(9u)-3*pi*a^2*exp(5u))*exp(-pi*a^2*exp(4u)), a>0
phi_a(u)=a^(-1/2)*phi_1(u+(log a)/2)
integral_0^infinity phi_a(u)da=0 and integral_0^infinity g_a(u)da=0
```

The last identity is an exact cancellation of the second and fourth
Gaussian moments. It suggests subtracting the continuous index cell
that ends at each integer:

```text
r_n(u)=g_n(u)-integral_(n-1)^n g_a(u)da
Phi(u)=sum_(n>=1)r_n(u)
```

## Normal Convergence

The continuous translate law differentiates to

```text
partial_a phi_a(u)=a^(-3/2)/2*(phi_1'(v)-phi_1(v)), v=u+(log a)/2
integral_R |u|^k*|partial_a g_a(u)|du<=C_k*a^(-3/2)*(1+|log a|^k), a>=1
sum_(n>=1)integral_R |u|^k*|r_n(u)|du<infinity for every integer k>=0
```

Thus the renormalized kernel series converges in every polynomially
weighted `L1` norm. Its cosine transforms may be differentiated term
by term to every fixed order:

```text
J_n(x)=integral_0^infinity r_n(u)cos(xu)du=I_n(x)-integral_(n-1)^n I_a(x)da
H_0^(k)(x)=sum_(n>=1)J_n^(k)(x), k=0,1,2,...
```

## Euler-Zeta Assembly

For `s=(1+i*x)/2`, set

```text
e_n(s)=n^(-s)-(n^(1-s)-(n-1)^(1-s))/(1-s)
J_n(x)=1/4*(Qhat(x)*e_n(s)+Qhat(-x)*e_n(conjugate(s)))
sum_(n>=1)e_n(s)=lim_(N->infinity)(sum_(n=1)^N n^(-s)-N^(1-s)/(1-s))=zeta(s), 0<Re(s)<1
sum_(n>=1)J_n(x)=1/2*Re(Qhat(x)*zeta((1+i*x)/2))=xi((1+i*x)/2)/8=H_0(x)
```

This is now an ordinary convergent transform decomposition derived
from the kernel, not the divergent termwise Bessel sum. At the origin
each renormalized block has the same strict sign:

```text
J_n(0)=Qhat(0)/2*(n^(-1/2)-2*(sqrt(n)-sqrt(n-1)))
2*(sqrt(n)-sqrt(n-1))=integral_(n-1)^n a^(-1/2)da>n^(-1/2)
J_n(0)>0 for every n>=1
```

## Coupled Matrix

Normal convergence through the second jet gives

```text
M_(m,n)(x)=J_m'(x)*J_n'(x)-(J_m(x)*J_n''(x)+J_n(x)*J_m''(x))/2
L[H_0](x)=sum_(m,n>=1)M_(m,n)(x)
```

The identity retains all mixed signs. No positive-semidefinite or
diagonal-dominance theorem is asserted.

## Positive-Time Boundary

The endpoint renormalization is not Newman-time compatible:

```text
n^2-integral_(n-1)^n a^2 da=n-1/3
r_n(u)=-(pi/2)*(3*n-1)*exp(-5u)+O_n(exp(-9u)) as u->infinity
For every fixed n>=1 and every t>0, exp(t*u^2)*r_n(u) is not integrable and has no ordinary half-line Fourier transform.
```

So the endpoint matrix cannot be deformed block by block to `t>0`.

## Live Handoff

Construct a t-compatible modular grouping whose individual blocks cancel the full negative-side theta tail before multiplication by exp(t*u^2), or work directly with the endpoint-subtracted theta primitive A_t and prove positivity of its coupled curvature expression.
