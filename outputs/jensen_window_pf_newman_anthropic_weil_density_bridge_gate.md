# Anthropic Weil Rank-Trace Density-Bridge Audit

Date: 2026-08-10

Status: external theorem statement and sources audited; central constant and finite-dimensional rank-trace assembly independently checked; density bridge open; not a proof of RH.

## Exact External Claim

The newly released paper claims unconditionally that

```text
liminf N_0^s(T,2T)/N(T,2T) >= c_*
c_*=c_*=3/2-cot(1/sqrt(2))/sqrt(2)=2-1/c_1*
c_*=0.67250070367941164573437979080329518859340302862626407892958886987266196710366666
```

Here `N_0^s` counts simple zeros on the critical line and `N` counts all
nontrivial zeros with multiplicity.  The result is a newly released external
claim.  The source paper, informal note, and public Lean companion make it a
serious input, but this local gate does not replace broad independent scrutiny
or independently reprove every analytic estimate.

## Rank-Trace Core

For P>=0, rank(P)<=r and n_+(Q)<=b, ||P+Q||_F^2 >= c tr(P)-c^2 r/4+2c tr(Q)-c^2 b.

At `c=2`,

```text
r>=2tr(P)+4tr(Q)-4b-||P+Q||_F^2.
```

The equality model is `P=(c/2)Pi_1 and Q=c Pi_2 for orthogonal projections of ranks r and b.`  Combining the paper's
zero-side block structure with its first and second moments gives

```text
tr(A)~N, ||A||_F^2~4N/3, and tr(P)+2b<=N give rank(P)>=(2/3-o(1))N.
```

This gate checks that algebra and the optimized constant.  It treats the
explicit-formula localization and prime-side moment evaluation as audited
external inputs rather than claiming a second complete proof.

## Exceptional Cap

Put

```text
E(T)=N(T,2T)-N_0^s(T,2T).
```

The claimed theorem gives

```text
limsup E(T)/N(T,2T) <= delta_*=1-c_*
delta_*=0.32749929632058835426562020919670481140659697137373592107041113012733803289633334
delta_*/2=0.16374964816029417713281010459835240570329848568686796053520556506366901644816667
```

Every off-line zero is exceptional, but exceptional also includes a repeated
on-line zero.  Therefore an off-line family can contradict the theorem only
by exceeding the same upper cap.

## Required Density Amplification

The exact candidate bridge is:

```text
If one off-line zero exists, prove that for some epsilon>0 and an unbounded sequence T_j, N_off(T_j,2T_j)/N(T_j,2T_j)>=delta_*+epsilon.
```

If symmetric pairs are counted instead, the target is:

```text
If C_off counts symmetric off-line pairs with multiplicity, force C_off(T_j,2T_j)/N(T_j,2T_j)>delta_*/2.
```

This is a very large requirement.  Since

```text
N(T,2T) ~ T log(T)/(2 pi),
```

N(T,2T)~T log(T)/(2 pi). O(1), O(log T), O(T^alpha) for alpha<1, and even O(T) descendants all have zero relative density.  In particular, one collision per height unit is still
zero density.  An infinite descendant chain is nowhere near sufficient.

## Comparison With Our Corpus

- Fixed-root shift tangency: The fixed-root shift chain propagates a common Jensen root and a finite Hankel defect; it does not manufacture distinct zeta zeros in height.
- Polar collision cascade: The polar cascade increases Jensen degree and multiplicity at one scaled root; degree growth is not zero-count growth.
- Cofinal contact route: The cofinal boundary programme seeks direct absence of contacts. Its local index is sign-definite but presently has no lower count proportional to T log T.
- Active Airy route: The current Airy branch controls a grouped source remainder. It supplies neither a dyadic off-line-zero count nor an evolution law for a Weil compression.

Thus none of the current exact cascades closes the new density target.  Calling
Jensen-degree growth, shift growth, or local contact index a density of zeta
zeros would be a category error.

## Heat-Flow Inertia Option

The other possible interface would be a heat-dependent family of Weil
compressions whose negative inertia is quantitatively observable from the
source side.  At `t=0`, the classical explicit formula supplies the prime-side
representation used in the paper.  Away from `t=0`, our present heat-flow
source formulas do not supply the corresponding trace and Hilbert-Schmidt
moment theorem.  A naive congruence transport would preserve inertia and is
incompatible with the collision mechanism, so no such transport is promoted.

A useful future theorem would need either:

1. a differential inequality for negative spectral mass of a heat-dependent
   Weil compression with source-side errors under control; or
2. a uniform observability theorem forcing any off-line pair to contribute a
   non-negligible negative eigenvalue across enough independent windows.

Neither theorem is currently available, and sparse zeros arbitrarily close to
the line are precisely what first and second aggregate moments fail to see.

## Route Decision

Do not pivot the main proof programme.  Preserve this result as a strong new
external constraint and a precisely stated optional bridge.  Continue the
grouped Airy amplitude, contour-tail, interior-join, and `T_upper` assembly.
Reopen the density branch only if a mechanism produces order `T log T`
distinct off-line zeros in dyadic windows, or if a genuine heat-dependent
negative-inertia estimate appears.

## Pi Provenance

The proportion constant `c_*` contains no `pi`; it comes from the optimized
Montgomery-Taylor cosine window with parameter `theta=1/sqrt(2)`.  The `pi` in
the density scale comes from the Riemann-von Mangoldt main term
`N(T,2T)~T log(T)/(2 pi)`.  In the paper's Gabor compression, `2 pi` also comes
from the Fourier sampling spacing `h=2 pi/L`.  These are separate appearances
and are not inserted as free normalization constants.

## Proof Boundary

This audit proves the constant identities, rank-trace equality model,
exceptional-cap logic, and zero-density scale guard.  It does not independently
reprove every analytic estimate in the external paper, establish community
acceptance, prove a defect-density amplification theorem, transport Weil
inertia through the heat flow, complete the grouped Airy remainder, prove
`Lambda<=0`, prove RH, or establish a prize-level conclusion.
