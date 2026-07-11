# Jensen-Window PF Heat-Flow Monotone-Closure Scout

Date: 2026-07-06

Status: finite heat-flow closure scout. This is not a proof of a closed
differential inequality, monotone-contraction theorem, Jensen-window
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_heat_flow_monotone_closure_scout`.

Proof boundary: exact lambda-flow algebra plus finite Arb diagnostics only.
This artifact does not prove the monotone-contraction theorem for all zeta
windows, does not construct a Cauchy-Binet kernel, does not prove `jwpf_06`,
and does not prove `Lambda <= 0`.

Machine-readable scout:

```text
work/rh_compute/results/jensen_window_pf_heat_flow_monotone_closure_scout.json
```

Generator:

```text
python work/rh_compute/scripts/jensen_window_pf_heat_flow_monotone_closure_scout.py
```

Checker:

```text
python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_monotone_closure_scout.py
```

Current result:

```text
validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues
```

## Exact Flow Algebra

With

```text
A_k(lambda) = mu_{2k}(lambda)*k!/(2*k)!
r_k = A_{k+1}/A_k
x_k = r_k/r_{k-1}
```

the heat-flow derivative is:

```text
dA_k/dlambda = 2*(2*k+1)*A_{k+1}
```

So the adjacent contraction evolves by:

```text
dlog(x_k)/dlambda
  = 2*r_k*((2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1))
```

and the monotone-gap log ratio evolves by:

```text
dlog(x_{k+1}/x_k)/dlambda
  = 2*r_k*((2*k+5)*x_{k+1}*x_{k+2}
            - 3*(2*k+3)*x_{k+1}
            + 3*(2*k+1)
            - (2*k-1)/x_k)
```

## Boundary Threshold

At a boundary point `x_{k+1}=x_k=q`, assuming only the next monotone
inequality `x_{k+2}>=q`, the bracket has lower bound:

```text
((q-1)^2*((2*k+5)*q-(2*k-1)))/q
```

Thus a heat-flow invariance proof from this direct bracket route needs the
additional lower-bound subtarget:

```text
q >= (2*k-1)/(2*k+5)
```

This threshold is supplied by ordinary raw-moment log-convexity after keeping
the factorial normalization explicit:

```text
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues
```

## Finite Diagnostics

The existing coefficient enclosure cache gives:

```text
lambdas: 0.0, 0.000001, 0.0001, 0.01, 0.1
max coefficient index: 64
threshold rows: 315
threshold positive rows: 315
flow-bracket rows: 305
flow-bracket positive rows: 305
```

Weakest threshold margin:

```text
3.822433850353900366E-2 at lambda=0.1, k=63
```

Weakest actual flow bracket:

```text
2.665474275226736438E-4 at lambda=0.1, k=61
```

The finite result supports the heat-flow route but does not close it. The
threshold subtarget is now exact; the remaining analytic upgrade is a global
flow-invariance argument, adjacent-log-concavity control, and initial or
terminal conditions for the monotone-contraction cone.

## Integration

This scout sharpens:

```text
outputs/jensen_window_pf_monotone_contraction_theorem_target.md
outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md
outputs/jensen_window_pf_monotone_contraction_stress.md
outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md
```

Summary:

```text
The heat-flow route reduces to exact identities for the contraction flow and a concrete boundary threshold q >= (2*k-1)/(2*k+5). The existing finite zeta cache clears that threshold on 315 Arb-classified rows and has positive actual flow brackets on 305 rows. The boundary threshold is discharged exactly by raw moment log-convexity in the separate boundary-threshold lemma, but this scout is still not a closed analytic differential inequality.
```
