# Jensen-Window PF Newman Modular-Blend Adaptive-Saddle Gate

Date: 2026-07-16

Status: exact saddle geometry with finite frequency-adaptive diagnostics.
This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.py
```

## Saddle Law

The two exponential components of each theta summand have `a=9` and
`a=5`. After the Newman weight and spectral phase,

```text
F_(a,n,t,x)(u)=t*u^2+(a+i*x)*u-pi*n^2*exp(4u), a in {5,9}
2*t*u+a+i*x-4*pi*n^2*exp(4u)=0
```

At `t=0`,

```text
u_*=(1/4)*Log((a+i*x)/(4*pi*n^2))
Re(u_*)=(1/4)*log(sqrt(a^2+x^2)/(4*pi*n^2))
n_*(a,0,x)=(a^2+x^2)^(1/4)/(2*sqrt(pi))=sqrt(x/(4*pi))*(1+O_a(x^-2))
```

Thus the saddle crosses the modular blend's `u=0` transition when
`n` is on the Riemann-Siegel scale `sqrt(x/(4*pi))`.

At positive Newman time the crossing saddle is `u=i*y` and obeys

```text
At u=i*y: 4*pi*n^2*cos(4y)=a and 4*pi*n^2*sin(4y)=x+2*t*y.
4y=atan((x+2*t*y)/a), 0<y<pi/8; n_*^2=sqrt(a^2+(x+2*t*y)^2)/(4*pi)
For 0<=t<=1/2 and a in {5,9}, n_*(a,t,x)=sqrt(x/(4*pi))*(1+O(x^-1)) as x->infinity.
```

## Adaptive Scout

For each sampled point, the table compares the exact transition index
with the first blended block counts whose Laguerre expression and
strict monotonicity margin `-L_t'` are within 50% of their full
high-precision values:

| t | x | n* (a=5) | n* (a=9) | first N for L | first N for -L' | ceil(max n*)+2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 80 | 2.525593 | 2.531078 | 3 | 3 | 5 |
| 0 | 100 | 2.822709 | 2.826643 | 4 | 4 | 5 |
| 0 | 120 | 3.091534 | 3.09453 | 4 | 4 | 6 |
| 0 | 150 | 3.455901 | 3.458047 | 5 | 5 | 6 |
| 0 | 200 | 3.990046 | 3.991441 | 6 | 6 | 6 |
| 0.5 | 80 | 2.531516 | 2.536771 | 3 | 3 | 5 |
| 0.5 | 100 | 2.828058 | 2.83183 | 4 | 4 | 5 |
| 0.5 | 120 | 3.096446 | 3.099322 | 4 | 4 | 6 |
| 0.5 | 150 | 3.460321 | 3.462383 | 5 | 5 | 6 |
| 0.5 | 200 | 3.993897 | 3.995238 | 6 | 6 | 6 |

Both observed counts agree row by row and grow from three to six, while
a two-index collar beyond the larger `a=5,9` transition covers all ten
rows. This is a stable truncation diagnostic only. The exact Xi Lehmer
counterexample shows that the sampled `-L_t'` signs are not globally
representative.

## Decision

The fixed finite truncation is closed, but the blend identifies its
correct replacement: choose `N=N(x,t)` on the transition scale and
prove a sign-aware remainder after retaining both theta components.
That task now concerns `L_t` only. Arb certifies `-L_0'<0` at the
Lehmer stress point, and time continuity rejects the proposed global
monotonicity condition. Near close pairs, corrected `C1` double-zero
transversality remains the lower-derivative Newman target.
