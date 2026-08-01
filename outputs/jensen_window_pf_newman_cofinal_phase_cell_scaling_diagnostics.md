# Newman Cofinal Phase-Cell Scaling Diagnostics

Date: 2026-07-25

Status: reproducible finite diagnostics and exact shell-scaling
identities. This is not a cofinal theorem and not a proof of
`Lambda<=0` or RH.

## Exact Shell Geometry

```text
t_j=1/(5*j), R_j=j+38
t_j*R_j=1/5+38/(5*j) -> 1/5
t_j-t_(j+1)=1/(5*j*(j+1))
t_j*log(R_j/(4*pi)) -> 0
```

The cofinal bottom edge remains in the nonuniform t*L->0 layer. The proved dominant-saddle ray t*L>=25 cannot certify it.

## Branch Stability

| comparison | agreements | panels | disagreements |
|---|---:|---:|---:|
| Q207 lowest-time 2D cells vs Q208 bottom | 408 | 414 | 6 |
| Q208 bottom vs Q208 top | 470 | 492 | 22 |

A branch records which coordinate supplied the stronger interval
separation. Agreement is useful evidence about certificate
conditioning, but it is not a homotopy or no-contact theorem.

## Q208 Cells

| edge | cells | value | derivative | transitions | max run |
|---|---:|---:|---:|---:|---:|
| bottom | 492 | 313 | 179 | 81 | 46 |
| top | 492 | 319 | 173 | 81 | 46 |

Both exact witness chains have positive-ray crossing count
`-19`. The complete closed
Q208 artifact, not this diagnostic, proves winding zero.

## Scaling Tests

The tested positive derivative scales are
`1`, `max(1,L/2)`, and `max(1,L)`, where
`L=log(x/(4*pi))`. Positive scaling preserves the exact
zero set and winding, but the metrics below are floating
conditioning diagnostics only.

| edge | scale | min relative clearance | max relative uncertainty | min tail domination | max panel turn |
|---|---|---:|---:|---:|---:|
| bottom | unit | 0.248707 | 1.14503 | 1.28936e+25 | 1.34753 |
| bottom | half_log | 0.248707 | 1.14503 | 1.02083e+25 | 1.07341 |
| bottom | full_log | 0.21996 | 1.14503 | 5.85046e+24 | 1.32919 |
| top | unit | 0.248838 | 1.1427 | 1.34421e+25 | 1.15721 |
| top | half_log | 0.248838 | 1.1427 | 1.06426e+25 | 1.07375 |
| top | full_log | 0.248838 | 1.1427 | 6.09936e+24 | 1.35255 |

## Route Decision

The half-unit phase cells remain well inside the exact
convex-cell interface at Q208, and the six-term arithmetic
tail is far smaller than the certified first-jet clearance.
The hard quantity is retained first-jet separation in the
small-time layer, not omitted-tail control.

The next theorem should therefore be an adiabatic successor
criterion: transfer the old bottom path across the
`1/(5*j*(j+1))` time collar using heat-equation derivative
bounds, then close the one-unit right strip by a uniform
phase or derivative cone. Computing Q209 alone would not
supply that theorem.

## Proof Boundary

The exact identities describe the chosen exhaustion and the source-derived counts reproduce finite Q207/Q208 certificates. The logarithmic scaling, angular, and conditioning summaries are diagnostics. No uniform time-collar displacement bound, no all-j right-strip cone, no Q209 theorem, no cofinal boundary theorem, no Lambda<=0, and no RH proof is supplied.
