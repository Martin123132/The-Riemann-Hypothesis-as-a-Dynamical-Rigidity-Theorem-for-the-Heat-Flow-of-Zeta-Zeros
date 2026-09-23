# Growing arithmetic groups: a collective centering theorem

22 September 2026. TWO HANDS NETWORK LTD. Private research. No publication.

## What advanced

The predecessor centered one fixed pair. This stage derives centering for
genuinely growing groups of actual Dirichlet terms, with their nonlinear
cross terms retained. For N=k^q, fixed q and finitely many fixed power layers,
the group may contain indices c k^h with

    sqrt(C_N) log(C_N+1) = o(log N),   1 <= c <= C_N.

The original S-roughness, both quarter windows, whole profile I=[1/8,2^24],
norm atom, fixed-height coefficients, gate and GLOBAL prime subtraction remain.
The estimate is uniform over height-independent masks, including the existing
off-band deletion. It does not assume independent phases for composite c's.

One explicit choice admits order (log N)^2/(log log N)^4 candidate indices
per layer. After natural group-amplitude scaling, its centered response is
O(1/log log N)+N^(-2/3+o(1)). Constants and onset depend on the fixed prime set.
This is an analytic asymptotic assertion, not a numerical onset certificate.

## Why the whole-group argument works

Instead of adding isolated pair responses, one polynomial is applied to the
ENTIRE field. It preserves the angular mode and every cross term and stays
well-defined at source zeros.

The first construction has paid error O(B^2/m) on |U|<=B. A sharper
Chebyshev construction improves that to O(B^2/L^2). The third angular mode
requires removing a linear term EXACTLY; otherwise it would not be a genuine
polynomial at zero. That correction is included and proved.

For odd L, the precise radius-one errors are

    E1 = 2/[pi L(L+2)],
    E3 = (6L+4)/[pi L^2(L+2)].

The polynomials have total degrees L and2L-1. Unique factorization and the
global prime subtraction then give a quantitative oscillatory estimate for
the whole growing group, not merely for each term in isolation.

Full proofs and all scope conditions are in DERIVATION.md, especially
sections4,7,8. The proofs are human-readable analytic arguments, not a
machine-formalized theorem library.

## Finite checks and their limit

The protocol fixed the complete retained N=512 field, p=257, x=1, both
original windows, radius4 and degrees8/32/128 before evaluation. Four
whole-window Arb builds use two precisions and two grids. The exact response
agrees with the authenticated predecessor in both windows.

The degree-128 binomial polynomial tracks the response closely, but its
guaranteed uniform payment is about0.00031908. That is too large to certify
either sign after restoration. This negative finite outcome is retained;
the radius and protocol were not fitted to it.

The sharper Chebyshev construction was checked separately at L=3,9,33,129.
At L=129 its mode1/mode3 errors are about0.0000376721/0.000113601 at unit
radius. These are smaller error constants, NOT a newly evaluated integrated
source response. No finite sign is claimed for the supplemental construction.

The finite controls use ONE prime and ONE profile, not the original S or a
numerical integral over I. The asymptotic theorem's full-profile statement
comes from its uniform analytic estimates.

## What is still open

The group is growing but remains sparse relative to all N source terms.
For the complete spectrum, the improved envelope payment is O(N/(d L^2)),
while the elementary frequency-separation argument deteriorates with degree.
This stage does not close that conflict or prove the full signed bound.

The next target is a collective estimate for the complete polynomial twisted
moments, or an amplitude-adapted approximation with a rigorously paid remainder.
Neither a sum of nonlinear group responses nor an extrapolation from the
one-prime finite signs is valid.

The core target .00092027989 and total budget .0011250999 are unchanged.
Inherited exterior payments remain attached exactly once. RH remains open.

## Evidence and reproducibility

CHECKS.json, INDEPENDENT_CHECKS.json and SPECTRAL_CHECKS.json retain exact
coefficient identities and numerical controls. VERIFICATION.json checks
complete cache coverage, two-grid overlaps, predecessor agreement and nine
tamper rejections. ACCEPTANCE.json binds both constructions without confusing
their distinct finite evidence. The original scalar/source tests and the
supplemental tests total406 named checks; integration acceptance adds559.

The source integrations share Arb. The supplemental exact-rational and
high-precision Cartesian checks are separately coded, but are not independent
special-function implementations. Finite tests do not prove the all-N claims.

One intake filename collision was repaired additively. One supplemental
checker stopped on a syntax error before calculation; its source, log and
receipt remain. ITERATIONS.md records both. No mathematical input was repaired
using observed outcomes.

The sealed archive is verified by MANIFEST.json. External CLEAN_REPLAY.json
records fresh-extraction rebuilds and byte comparisons; FINAL_CHECKPOINT.json
records the historical integrity audit. No provider calls, Forge experiments,
GitHub updates, external publication or unrelated process changes were made.
