# Unequal-truncation saddle map for the joined equation-(4) packet

Date: 2026-08-27

Status: exact parameter-map and route-split certificate; quantitative joined-packet enclosure remains open.

## Exact map

Put `p=t/(2*pi)`.  The real phase of one Fresnel label is stationary when

```text
ell=2*y+t/(pi*y)=2*(y+p/y).
```

Therefore `alpha=max(y,p/y)` and `beta=min(y,p/y)` satisfy

```text
t=2*pi*alpha*beta,
lambda=sqrt(alpha/beta),
ell=2*(alpha+beta).
```

This is exactly the unequal-truncation hyperbola used by O'Sullivan, not a fitted analogy.

## Actual-height split

For Theorem 1.5, the ratio of successive displayed remainder scales at `sigma=1/2` is

```text
rho_5=(lambda+lambda^(-1))^3/sqrt(t).
```

At `t=10^10`, `rho_5=1` occurs at

```text
lambda=46.3943339799335062266842987645385937090384426,
y=859.894401273188520021692227054209930487641234.
```

The first half-integer boundary not below this crossover is `860.5`; the previous one, `859.5`, still has ratio above one.  Hence a completed-series implementation has the exact cell split

```text
endpoint quarantine: 622..860   (239 cells),
contracting core:     861..39936 (39076 cells).
```

This split is a route audit, not an error bound.  The paper's constants are implicit.

## Route decision

The better donor is Theorem 2.1 / equation (2.42), before the completed symmetric re-expansion.  Its displayed `sigma=1/2` error scale is `lambda^(1/2)t^(-1/4-N/6)` and the stated implied constant is independent of `lambda` for `lambda>=1`.  We must still derive an identity for the fully reassembled heat-flow packet and make every constant explicit.  O'Sullivan's `R(s;alpha,beta)` is not identified here with `J_Z`.

Primary source: https://arxiv.org/abs/1811.01130

Pi provenance: every `pi` above comes from the fixed quadratic Riemann-Siegel/Fresnel phase and the standard relation `t=2*pi*alpha*beta`; none is fitted or inserted geometrically.

## Proof boundary

Exact stationary-phase/unequal-truncation parameter map, actual-height displayed-scale crossover, and a cancellation-preserving route split only. No identity between O'Sullivan's classical zeta remainder and J_Z, no explicit remainder constant, no numerical enclosure of the joined packet, J_Z, or D_K, and no non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion.
