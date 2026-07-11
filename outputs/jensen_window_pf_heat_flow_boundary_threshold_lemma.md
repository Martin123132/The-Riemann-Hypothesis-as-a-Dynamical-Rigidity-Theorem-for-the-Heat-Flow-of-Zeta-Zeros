# Jensen-Window PF Heat-Flow Boundary-Threshold Lemma

Date: 2026-07-06

Status: exact boundary-threshold lemma. This is not a proof of the
monotone-contraction theorem, Jensen-window PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_boundary_threshold_lemma`.

Proof boundary: this artifact proves only the lower-bound threshold needed by
the heat-flow boundary bracket. It does not prove adjacent log-concavity
`x_k<=1`, does not prove `x_{k+1}>=x_k` for all zeta windows, does not prove
a closed differential inequality, does not prove `jwpf_06`, and does not prove
`Lambda <= 0`.

Machine-readable lemma:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_heat_flow_boundary_threshold_lemma.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_boundary_threshold_lemma.py
```

Current result:

```text
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

## Exact Lemma

For `u>=0`,

```text
Phi(u)=sum_{n>=1} pi*n^2*exp(5*u)*(2*pi*n^2*exp(4*u)-3)*exp(-pi*n^2*exp(4*u)) > 0
```

because `2*pi*n^2*exp(4*u)-3 >= 2*pi-3 > 0`.

Set:

```text
w_lambda(u)=2*exp(lambda*u^2)*Phi(u)
mu_{2k}(lambda)=integral_0^infty u^(2*k) w_lambda(u) du
A_k(lambda)=mu_{2k}(lambda)*k!/(2*k)!
x_k=(A_{k+1}/A_k)/(A_k/A_{k-1})
```

Cauchy-Schwarz for the positive weight gives raw moment log-convexity:

```text
mu_{2k}(lambda)^2 <= mu_{2k-2}(lambda)*mu_{2k+2}(lambda)
```

The factorial normalization gives:

```text
x_k = ((2*k-1)/(2*k+1))*(mu_{2k+2}*mu_{2k-2}/mu_{2k}^2)
```

Therefore:

```text
x_k >= (2*k-1)/(2*k+1) > (2*k-1)/(2*k+5)
```

So the heat-flow boundary threshold from
`outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md` is discharged.

## Boundary Bracket

At a monotone-gap boundary `x_{k+1}=x_k=q`, with the next monotone-neighbor
condition `x_{k+2}>=q`, the heat-flow bracket has lower bound:

```text
((q-1)^2*((2*k+5)*q-(2*k-1)))/q
```

The exact threshold makes this lower bound nonnegative. This is only a
boundary pointing condition; it is not yet a global flow-invariance theorem.

## Finite Sanity Check

The existing coefficient enclosure cache gives:

```text
lambdas: 0.0, 0.000001, 0.0001, 0.01, 0.1
max coefficient index: 64
strong-threshold rows: 315
strong-threshold positive rows: 315
heat-threshold rows: 315
heat-threshold positive rows: 315
```

Weakest strong-threshold margin:

```text
8.170843282044743874E-3 at lambda=0.1, k=63
```

Weakest heat-threshold margin:

```text
3.822433850353900366E-2 at lambda=0.1, k=63
```

## Remaining Gap

```text
The heat-flow boundary threshold is not an open subtarget: Phi positivity and Cauchy-Schwarz raw-moment log-convexity prove the stronger bound x_k >= (2*k-1)/(2*k+1), which implies the needed x_k >= (2*k-1)/(2*k+5). The remaining heat-flow route still needs adjacent log-concavity, a global flow-invariance argument, and initial or terminal conditions.
```

Integration:

```text
outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
```
