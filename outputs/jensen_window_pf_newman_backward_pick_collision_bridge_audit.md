# Newman Backward-Pick Collision-Bridge Audit

Date: 2026-07-24

Status: exact collision countermodel and literature audit. This
is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_backward_pick_collision_bridge_audit.json
python work/rh_compute/scripts/jensen_window_pf_newman_backward_pick_collision_bridge_audit.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_backward_pick_collision_bridge_audit.py
```

Current result:

```text
validated Newman backward-Pick collision-bridge audit: 12 rows, 0 issues, 3 exact flow identities, 2 Pick-sign regions, 1 square-root speed blowup, 1 cutoff-hiding theorem, 1 noncommuting-limit obstruction, 1 preprint gap, 1 open Xi repair
```

## Exact Collision Model

Take

```text
E_t(z)=z^2+a-2t, a>0
partial_t E_t=-partial_z^2 E_t
t_*=a/2
t>t_*: zeros are +/-sqrt(2t-a), both real
t=t_*: z=0 is a double zero
t<t_*: zeros are +/-i*sqrt(a-2t), so one lies in C+
```

This is an even real entire Cartwright-class solution of exactly
the same backward heat equation. Its logarithmic derivative obeys

```text
g_t(z)=-E_t'(z)/E_t(z), h_t=Im(g_t)
partial_t g=-partial_z^2 g+2g*partial_z g
h_t(x,y)=2y*(x^2+y^2-b)/((x^2-y^2+b)^2+4x^2y^2), b=a-2t
```

Above the collision its imaginary part is positive throughout
the upper half-plane. Below the collision it is negative on the
newborn half-disc:

```text
For t<t_* one has h_t<0 on x^2+y^2<a-2t inside C+.
At x=0, y=sqrt(a-2t)/2, h_t=-4/(3*sqrt(a-2t))<0.
```

## Uniformity Failure

The upper zero has

```text
rho_+(t)=i*sqrt(a-2t), t<t_*
|rho_+'(t)|=1/sqrt(a-2t), which tends to infinity as t increases to t_*.
```

Thus a speed bound proved on each closed collision-free interval
does not produce one finite constant on an interval ending at the
collision.

There is a second, independent cutoff obstruction:

```text
For a lower cutoff y>=rho, the weighted negative energy is exactly zero whenever 0<t_*-t<rho^2/2, although h_t^- is nonzero below y=rho.
delta_rho=rho^2/2
The zero-energy bridge width tends to zero with rho. One cannot first prove a rho-dependent bridge and then send rho to zero while retaining a fixed backward time interval.
```

For every fixed `rho`, the weighted energy can therefore vanish
on a two-sided collision neighbourhood while the unweighted
negative part is already nonzero below `y=rho`. Sending `rho` to
zero also sends the guaranteed backward width to zero.

## Literature Audit

Kevin Schatz, Riemann Hypothesis: Backward Parabolic Positivity Barriers for the Xi Flow, preprint, 2025

DOI: https://doi.org/10.5281/zenodo.17636625

Manuscript: https://kschatz.github.io/rh-xi-backward-parabolic-barrier/schatz_riemann_hypothesis_backward_parabolic_positivity_barriers_xi_flow.pdf

Lemma 4.7 assumes a closed collision-free window; Lemma 7.3 applies its speed-dependent energy bound on an open interval ending at a collision and concludes global upper-half-plane positivity from a cutoff-supported energy.

The quadratic flow satisfies the generic heat, Burgers, Cartwright, top-time Pick, and isolated-collision inputs but violates the claimed collision bridge. The cited preprint therefore does not presently supply an admissible proof of backward Pick positivity.

## Live Handoff

A viable Xi-specific backward argument must add an estimate uniform in zero separation and lower cutoff through collisions, or an arithmetic invariant that excludes the collision before the cutoff limit. Closed-window speed bounds, weighted-energy continuity, and top-time Pick positivity are insufficient.
