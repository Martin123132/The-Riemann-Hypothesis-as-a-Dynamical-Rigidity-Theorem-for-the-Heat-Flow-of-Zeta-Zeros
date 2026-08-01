# Newman Multiplicity-Compatible Adaptive-Jet Benchmark

Date: 2026-07-25

Status: exact hypothetical heat-flow benchmark, theorem-domain correction, and open Xi handoff; not a proof of an all-stage successor, `Lambda<=0`, or RH.

## Why The Uniform Floor Is Too Strong

Suppose a continuous first jet `V(t,x)` at one fixed `x` obeys

```text
||V(t,x)||>=b(x)>0 for every 0<t<=t_0.
```

Taking `t` down to zero gives `||V(0,x)||>=b(x)`. Thus a positive
fixed-`x` floor uniform through `0<tL<=c_*` excludes a multiple endpoint
zero. RH itself does not assert simplicity, and the cofinal successor does
not need such a floor. Its bottom times are `t_j=1/(5j)`, and an old cell
may have clearance tending to zero while remaining nonzero at every
positive stage.

The existing positive-boundary attainment and delta-localization gates
already establish this quantifier warning and the arbitrary-multiplicity
Hermite split. The new contribution below is to put that split into the
adiabatic successor coordinates and derive its relative transport cost.

## Multiple-Zero Thought Experiment

Let an endpoint profile have an `m`-fold local zero and evolve by the same
backward heat equation:

```text
P_m(t,y)=exp(-t D_y^2)y^m
        =m! sum_(k=0)^floor(m/2)
          (-t)^k y^(m-2k)/(k!(m-2k)!).             (1)
```

For `s=sqrt(2t)` and `z=y/s`, (1) is exactly

```text
P_m(t,y)=s^m He_m(z),                               (2)
```

where `He_m` is the probabilists' Hermite polynomial. Its zeros are real
and simple. Since `He_m'=m He_(m-1)`, consecutive Hermite polynomials have
no common zero. Therefore an `m`-fold zero at `t=0` becomes `m` distinct
real zeros for every `t>0`; it does not create a positive-time contact.

This is precisely the kind of mathematically legitimate limiting
thought experiment that a simplicity-assuming proof would miss.

## Adaptive Jet

Use the time-dependent positive scaling

```text
W_m(t,y)=(P_m(t,y),s partial_y P_m(t,y)),
s=sqrt(2t).
```

Then both components have the same parabolic order:

```text
W_m=s^m(He_m(z),m He_(m-1)(z)).                    (3)
```

At fixed `y`,

```text
partial_t W_m=s^(m-2)(A_m(z),B_m(z)),              (4)
A_m=-m(m-1)He_(m-2),
B_m=m He_(m-1)-m(m-1)(m-2)He_(m-3),
```

with negative-index Hermite terms interpreted as zero. Define

```text
C_m=sup_(z in R)
 sqrt(A_m(z)^2+B_m(z)^2)
 /sqrt(He_m(z)^2+m^2 He_(m-1)(z)^2).              (5)
```

The denominator in (5) never vanishes, and the quotient tends to zero as
`|z|` tends to infinity, so `C_m<infinity`. Equations (3)-(5) give

```text
||partial_t W_m||/||W_m||<=C_m/(2t).               (6)
```

The singularity in (6) is the correct one: it permits the endpoint
clearance to vanish as a power of `t` without permitting a contact at
positive `t`.

## Cofinal Cost

Put `K_m=C_m/2`. Gronwall between consecutive bottom times gives

```text
||W_m(t_(j+1),y)||
 >=(t_(j+1)/t_j)^K_m ||W_m(t_j,y)||
 =(j/(j+1))^K_m ||W_m(t_j,y)||.                    (7)
```

The logarithmic cost is

```text
K_m log(1+1/j)->0.                                 (8)
```

Every factor in (7) is positive, although their infinite product may tend
to zero. That is enough for the successor induction, whose stages only
need positive-time origin exclusion.

For a general heat solution `H_t=-H_xx`, the same adaptive jet
`W=(H,sH_x)` obeys

```text
partial_t W=(-H_xx,H_x/s-sH_xxx),  s=sqrt(2t).     (9)
```

For every `t>0`, `diag(1,s)` lies in `GL+(2,R)`, so this scaling preserves
the contact set, orientation, and first-jet degree.

## Corrected Xi Strategy

The exact benchmark suggests a two-regime theorem.

```text
entry regime:
  certify each newly added high-frequency right strip by a corrected
  Riemann-Siegel half-plane, slope-gap, or equivalent arithmetic theorem;

descendant regime:
  transport previously certified cells down their later collars using
  the adaptive jet and a relative K/t estimate, or a rigorously comparable
  local factorization estimate.
```

This would use the history carried by the successor induction instead of
reproving one absolute lower bound on the entire old edge at every stage.
It is compatible with finite-multiplicity zeros at `t=0`.

## Proof Boundary

The finite heat-polynomial formula, Hermite rescaling, simple real
splitting, adaptive-jet identities, `O(1/t)` condition number, cofinal
relative transport cost, and fixed-`x` uniform-floor guard are exact. They
are a benchmark, not a local factorization theorem for Xi. No uniform
adaptive relative bound, entry-strip cone, all-`j` composition,
`Lambda<=0`, RH, PF-infinity, or Clay-prize conclusion is proved.

Machine-audited files:

```text
work/rh_compute/results/jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark.json
work/rh_compute/scripts/jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark.py
work/rh_compute/scripts/check_jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark.py
```
