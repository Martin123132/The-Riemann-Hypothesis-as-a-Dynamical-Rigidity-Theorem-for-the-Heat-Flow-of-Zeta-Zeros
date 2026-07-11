# Jensen-Window PF Newman Boundary-Energy Direction Gate

Date: 2026-07-11

Status: exact boundary-energy singularity and published-directionality
gate. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_boundary_energy_direction_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_boundary_energy_direction_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_boundary_energy_direction_gate.py
```

Current result:

```text
validated Jensen-window PF Newman boundary-energy direction gate: 10 rows, 0 issues, 1 universal gap law, 1 exact higher-jet birth model, 1 nonintegrable collision-energy asymptotic, 1 conditional exclusion criterion, 1 published-scope audit, 1 open Xi boundary-energy handoff
```

## Higher Gap Jet

The real-zero ODE gives

```text
x_k'=2*PV sum_(j!=k) 1/(x_k-x_j); for a pair gap g and q=g^2, q'=8-4*q*S, S=sum_out 1/((x_+-y)(x_--y)).
At a finite double-zero birth t_*, K=sum_out 1/(c-y)^2, so q(t)=8*(t-t_*)-16*K*(t-t_*)^2+O((t-t_*)^3).
```

For the exact classical-field countermodel from the previous gate,

```text
a^2=1 + 16/(pi + 8)
P(z)=(z^2-1)^2*(z^2-a^2)
F_tau=exp(-tau*d_z^2)P
B(1)=-pi/8
K=(pi**2 + 24*pi + 160)/64=3.83230981386319369
z_s=1+s*sqrt(2)*e-pi*e^2/4-s*sqrt(2)*K*e^3+P4*e^4+s*sqrt(2)*L5*e^5+O(e^6), e=sqrt(tau), s=+/-1
q(tau)=8*tau-16*K*tau^2+C3*tau^3+O(tau^(7/2))
C3=(pi**4 + 48*pi**3 + 800*pi**2 + 14848 + 5632*pi)/32=1313.21386617429988
K>0, q''(0)=-32*K<0, q'''(0)=6*C3>0
```

Thus the classical drift and the first two nontrivial gap-jet signs
still do not distinguish a forbidden positive boundary from ordinary
square-root birth.

## Boundary Energy

For a fixed positive reference gap `Delta`, Rodgers-Tao use

```text
V(r)=r^(-2)-1+2*(r-1)=r^(-2)+2*r-3 for r>0
V(1)=V'(1)=0 and V''(r)=6/r^4>0, hence V(r)>=0
E_pair=Delta^(-2)*V(g/Delta)=1/g^2+2*g/Delta^3-3/Delta^2
```

Substitution of the collision gap gives

```text
E_pair(tau)=1/(8*tau)+K/4-3/Delta^2+4*sqrt(2)*sqrt(tau)/Delta^3+O(tau)
dE_pair/dtau=-1/(8*tau^2)+2*sqrt(2)/(Delta^3*sqrt(tau))+O(1)<0 for small tau>0
integral_0^epsilon E_pair(tau) dtau=+infinity
```

This is the decisive local fact. The energy is finite at every time
after birth and decreases immediately, but it has a logarithmically
nonintegrable trace at the birth time. Nonnegative forward relaxation
therefore does not by itself forbid a boundary collision.

A valid conditional criterion remains:

```text
If a nonnegative renormalized energy contains the colliding interaction with a reference gap Delta bounded away from zero and a weight bounded below, then local time-integrability down to t_* excludes the collision, because the interaction is asymptotic to 1/(8*(t-t_*)).
```

## Published Direction

The primary source is Rodgers-Tao, arXiv:1801.05914v5.

```text
The paper assumes Lambda<0 and studies real zeros for Lambda<t<=0.
Theorem 11 states the zero ODE for Lambda<t<=0.
Theorem 17 proves integral_(Lambda/4)^0 E_tilde_[0.5*T*log(T),3*T*log(T)](t) dt=o(T*log_+(T)^3) as T tends to infinity.
For Lambda<0, Lambda/4-Lambda=3*abs(Lambda)/4>0; the theorem does not integrate down to the collision boundary t=Lambda.
The estimate controls relaxation inside an already-real interval and uses zeta-zero information at t=0. It is unavailable in a hypothetical Lambda>0 regime, where 0 lies before the real-zero boundary.
```

The theorem is exactly suited to proving relaxation between an
already-real negative time and `t=0`. It does not provide finite energy
on `(Lambda,Lambda+epsilon)` and cannot be imported into the opposite
hypothetical regime `Lambda>0`.

## Live Handoff

Prove an Xi-specific boundary-uniform estimate giving finite local integrated renormalized energy on (Lambda,Lambda+epsilon), or find another global invariant that forbids the universal 1/(t-Lambda) collision singularity, without assuming Lambda<0 or importing t=0 real-zero information.

This narrows the energy route to a boundary-trace theorem. Until that
theorem is proved from the Xi kernel, energy language is diagnostic rather
than a completed Newman obstruction.
