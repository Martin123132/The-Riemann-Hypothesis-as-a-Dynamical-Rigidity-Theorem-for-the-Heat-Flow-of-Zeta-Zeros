# Newman Positive-Boundary Delta Localization Gate

Date: 2026-07-24

Status: exact positive-boundary localization and route guard. The
delta-dependent compact-strip family remains open; this is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_positive_boundary_delta_localization_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_positive_boundary_delta_localization_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_delta_localization_gate.py
```

## Exact Localization

```text
If Lambda>0, then 0<Lambda<=1/5 and H_Lambda has a finite real multiple zero
Lambda<=0 iff H_t has only simple zeros for every 0<t<=1/5
To contradict Lambda>0 it is enough, for every rational delta in (0,1/5], to exclude common zeros of H_t and H_t' on delta<=t<=1/5; choose any rational delta<Lambda
```

For a fixed positive lower time, the complete dominant ray gives

```text
For fixed delta>0, t>=delta and L=log(|x|/(4*pi))>=25/delta imply tL>=25 and L>=50, so the dominant-saddle theorem gives H_t'(x)^2-H_t(x)H_t''(x)>0
Every multiple zero with delta<=t<=1/5 lies in |x|<R_delta, R_delta=4*pi*exp(25/delta)
```

Therefore one explicit cofinal family is

```text
For delta_j=1/(5j), j>=1, Lambda<=0 iff (H_t,H_t')!=(0,0) on [delta_j,1/5]x[-4*pi*exp(125j),4*pi*exp(125j)] for every j
```

| j | delta_j | L cutoff | compact radius |
|---:|---:|---:|---:|
| 1 | 1/5 | 125 | 4*pi*exp(125) |
| 2 | 1/10 | 250 | 4*pi*exp(250) |
| 3 | 1/15 | 375 | 4*pi*exp(375) |
| 4 | 1/20 | 500 | 4*pi*exp(500) |

These rectangles grow rapidly and are not currently certified.
The reduction removes the need for one margin uniform all the way
to `t=0`; it does not make the compact arithmetic problem easy.

## Endpoint-Uniformity Guard

```text
G_t(x)=x^2-2t solves partial_t G=-partial_x^2 G; for t>0 it has the simple real zeros +/-sqrt(2t), at t=0 it has a double zero, and for t<0 its zeros are nonreal
L[G_t]=G_t'^2-G_t G_t''=2x^2+4t>0 for t>0, while G_t(0)^2+(G_t'(0)/A)^2=4t^2 and at either root the same first-jet quantity is 8t/A^2; both margins tend to zero as t decreases to zero
```

This model has the correct backward-heat collision geometry:
positive-time zeros are real and simple, but their separation and
first-jet margins vanish at the endpoint. Hence

```text
A t-independent positive first-jet floor down to t=0 is a valid stronger sufficient target, but it is not logically necessary for Lambda<=0 or positive-time simplicity
```

The current uniform corrected small-ball target remains a valid
stronger route. It must not be presented as logically necessary
unless endpoint simplicity is explicitly added to the burden.

## Live Handoff

```text
A sufficient nonuniform route is: for each fixed delta>0, prove the corrected C1 finite-main inequality only on delta<=t<=1/5 and 38<|x|<4*pi*exp(25/delta), with constants allowed to deteriorate as delta decreases
Use delta-dependent compact or arithmetic transversality estimates, or keep the uniform scaled small-ball theorem as an optional stronger route; do not require endpoint simplicity without stating that extra burden
```
