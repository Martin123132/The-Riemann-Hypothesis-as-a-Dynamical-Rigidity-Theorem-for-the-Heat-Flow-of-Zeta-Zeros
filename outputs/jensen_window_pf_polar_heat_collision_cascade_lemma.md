# Jensen-Window PF Polar Heat-Collision Cascade Lemma

Date: 2026-07-11

Status: exact polar heat-collision cascade and unbounded-degree escape
theorem. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json
python work/rh_compute/scripts/jensen_window_pf_polar_heat_collision_cascade_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_polar_heat_collision_cascade_lemma.py
```

Current result:

```text
validated Jensen-window PF polar heat-collision cascade lemma: 10 rows, 0 issues, 3 exact local identities, 1 double-root criterion, 1 higher-multiplicity gate, 1 infinite polar cascade, 1 exponential-polynomial classification, 1 unbounded-degree escape theorem, 1 open scaled-tail handoff
```

## Local Polar Heat Jet

Write `Q=J_(d+1,n,lambda)` and let `xi!=0` be a root of
`J_(d,n,lambda)` of multiplicity `m`. Polar descent gives

```text
J_d=J_(d+1)-w*J_(d+1)'/(d+1)
If mult_xi(J_d)>=m and Q=J_(d+1), then Q^(ell)(xi)=(d+1)_ell*Q(xi)/xi^ell for 0<=ell<=m.
```

Combining this recurrence with the exact heat hierarchy gives

```text
For 0<=r<=m-2, (partial_lambda J_d)^(r)(xi)=(4*n+4*d+2)*(d)_r*Q(xi)/xi^(r+1).
```

The checker verifies 210 Taylor recurrences and
165 heat-jet instances over degrees 2 through 10.
At a double root the complete first-order real-splitting condition is

```text
At a double root xi, forward real splitting requires Q(xi)*J_d''(xi)/xi<=0.
```

For multiplicity at least three, differentiable one-sided hyperbolicity
forces `Q(xi)=0`; otherwise the constant perturbation produces a
nonreal root of the leading local polynomial.

## Infinite Cascade

If `Q` is itself negative-root hyperbolic, the polar reciprocal sum
shows that a multiple lower root must be a root of `Q`, with one extra
unit of multiplicity. Repeating upward gives

```text
If every J_D, D>=d, is hyperbolic and mult_xi(J_d)=m>=2, then mult_xi(J_D)=D-(d-m).
```

Put `L=d-m`. Then `J_D(w)=(1-w/xi)^(D-L) R_D(w)` with
`deg(R_D)<=L`. Scaling `w=z/D` and taking the locally uniform Jensen
limit proves

```text
The cascade forces F(z)=exp(-z/xi)*R(z) with deg(R)<=d-m; conversely this form gives multiplicity at least D-(d-m), with equality when deg(R)=d-m.
```

The converse follows from `J_D(w)=D! [t^D] exp(t)F(wt)`. The
executable sample checks the resulting factor through degree 10.

## Boundary Consequence

A Laguerre-Polya entire function with positive coefficients that is not
exponential-polynomial therefore has strictly hyperbolic Jensen
polynomials in every fixed degree. If a nearby parameter is outside
Laguerre-Polya, the least failing degree must tend to infinity as the
boundary is approached.

This removes fixed-degree first contact as a possible closing argument.
The live problem is a quantitative estimate in the rescaled `D->infinity`
collision layer.
