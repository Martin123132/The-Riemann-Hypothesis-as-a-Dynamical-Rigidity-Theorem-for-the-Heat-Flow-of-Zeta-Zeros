# Diagonal audit and short radial pairing

19 September 2026. TWO HANDS NETWORK LTD. Private local research.

## Result

Both supplied batches pass this audit. The diagonal replacement and the
low-amplitude cancellation are retained with their original common inverse
and weighting. The high-amplitude diagonal has the stated vanishing
payment. None establishes the full signed inequality.

The continuation adds a further deterministic reduction: only
O(log log N) radial inverse scales need the complete middle-amplitude
numerator. The original radial roster has O(log N) scales. Above the new
cutoff, the entire diagonal has an explicit o(log N) bound.

This reduces the size of the cancellation problem, not the length of the
Dirichlet source or the size of the remaining signed complementary term.

## New theorem

Set d=ceil(log(N+1)), ell=ceil(log_2(d+1)), S=d ell. Keep the original
radial nodes s_j=2^(j/2), and split their kernel at s_j>=S into L_<+L_>.
For C0=log(2)/[pi(1-1/sqrt(2))],

    L_>(r) <= C0/S

at every actual source amplitude. Thus any [0,1]-weighted high-node
diagonal is at most

    alpha gamma_128 C0 (d+1)^2 mu_I/(bar W_N S)
      = O_I(d/ell)=o_I(d).

The exact short-node bound is ceil(log_2 d)+2ceil(log_2 ell).
Examples of the integer allocation bound: 11 nodes at d=64, 16 at d=1024,
19 at d=4096. These are radial bookkeeping examples, not new source runs.

The remaining signed expression is exactly reduced, with finite errors, to

    V_pair = bar W^-1 int g_mid H L_< dnu,
    V_comp = bar W^-1 int (H-A)(g_hi L_< + a0 L_>) dnu.

Only V_pair requires the complete numerator H. Both its amplitude window
and its radial roster are controlled explicitly. V_comp remains signed
and coupled to that SAME physical source and inverse. It is not discarded.

The new high-node subtraction has an amplitude-Hessian coefficient
O(1/ell), so the reduction does not introduce a growing derivative cost.
On the already audited profile compact [1/2,3/2], its sharper normalized
coefficient is bounded by 0.0000304/ell asymptotically. This is NOT a
full-profile or full-source sign bound, and no numerical onset is supplied.

## Verification

- Incoming manifests: 21/21 and 11/11 payloads; four additional manifest
  and verification files authenticate.
- Four original successful checkers rerun unchanged. Four JSON results and
  two 2,048-cell interval ledgers reproduce byte-for-byte.
- The original failed checker remains preserved. Its successor changes
  only an incorrect overlap requirement for rounded tail upper bounds.
- New primary controls: 4,828 exact allocations, 1,600 exact derivative
  cases, 96 complex-rational ensembles of 24 points, 32 cross-precision
  interval overlaps, and six rejected wrong-accounting variants.
- New independent controls: 30 exact Q(sqrt(2)) lattice/derivative cases
  and 32 overlaps with rational-series enclosures, without Arb or mpmath.
- The archived final acceptance and clean-replay receipts record the
  extraction/replay and final historical-pin audit.

These are implementation and scalar proof controls. The asymptotic
deductions are written analytic arguments, not computer-formalized proofs.
No defining-source integrations, source scans, zero counts, provider calls,
new Forge experiments, repository push, or external publication occurred.

## Exact remaining obligation

The joint bound for V_pair+V_comp must still fall below 0.0011250999 per
log N, with the original window weight, fixed Euler polynomials, shared
core, common inverse, and inherited endpoint/count obligations retained.

The finite error is E_low+E_amp+E_node, not zero. Each is o_I(log N)
before endpoint restoration, so the asymptotic target receives no new
fixed charge. Separate absolute O(log N) bounds do not close that target.

The next useful work is the arithmetic sign of this short-block pairing
JOINTLY with its high complementary channels. More isolated diagonal
estimates cannot substitute for that sign. RH remains open.
