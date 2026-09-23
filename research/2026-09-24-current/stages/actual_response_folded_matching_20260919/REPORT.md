# Actual-source test and response-specific pairing

19 September 2026. TWO HANDS NETWORK LTD. Private local RH research.

## Outcome

The same complete exact-Euler all-node observable was evaluated on the
prescribed N=64 rough core, both original windows and I=[1/2,3/2]. All
1,048,576 height/profile cells in each final window are covered by
interval enclosures. There is no sampled-height selection.

Normalized by log(64), conservative outward upper bounds are:

| Reduced observable | First window | Second window |
| --- | ---: | ---: |
| Exact all-node response G | 0.00000474 | 0.00000164 |
| Finite-potential response G_(1/d) | 0.0000525 | 0.0000408 |

Both finite reduced observables are below the numerical comparison target
0.0011250999. This is NOT a full-source finite certificate: the compact
and cutoff are fixed, upstream asymptotic deletions/payments are not
given finite bounds here, and no unbounded or RH assertion follows.
The first finite-potential response is in fact strictly positive; the
required upper bound must not be confused with a claim of negativity.

## What changed mathematically

The response is odd in X and even in Y. Therefore the old full-opposite
matching cost can be replaced by

    integral [(X+X')^2 + (|Y|-|Y'|)^2] dpi,

with BOTH marginals still exactly the actual weighted source law. This
cost is no greater than the old cost. The signed Taylor and sign-free
matching bounds survive with their original constants. It is a weaker
sufficient condition, not an established property of the source.

An exact four-point countercontrol has zero complex mean, balanced
quadratic moments and zero signed response at every node, but its best
old matching cost is 8/3 while the new cost is zero. Thus full opposite
symmetry was an unnecessary proof requirement. A separate two-point
control still obstructs deducing the new condition from moments alone.

The local gradient also has a cubic-scale bound, giving a smaller paid
interval error below the first radial node. No third-moment hypothesis
is used or introduced.

## What the arithmetic test says about matching

For the explicit product coupling, the folded G upper bounds are
approximately 0.00521459 and 0.00239981 before division by log(64).
The first does NOT pass the numerical comparison, although direct G
does. The second passes, as does the old product bound in that window.
Neither global product-coupling G_h bound passes; both direct G_h bounds
do. These are failures of sufficient estimates, not source failures or
claims that an optimal folded coupling fails.

The implication is modest but useful: use the response's exact symmetry
when matching is helpful, but do not make matching an experiment-wide
or proof-wide requirement. Direct signed arithmetic remains the target.

## Verification and limits

- 2,400 exact gradient/symmetry controls and 1,000 exact product-law controls.
- Three exact primal/dual optimal-cost controls; 12 exact zero-response checks.
- Complete coarse/fine Arb cell calculations retained at 96, 128 and160 bits.
- Independent direct-phase evaluation at 85 decimal digits, and exact
  rational reaggregation of all 8,192 final height rows.
- Fresh extracted-package replay and historical-pin audit are recorded in
  the delivery receipts. Raw rows, source and configuration hashes are sealed.

Two initial numerical-wrapper failures are retained. The first exported
no rows because a radius square root included zero with interval padding;
the second saved all512 rows but failed an aggregate diagnostic square
root. Successor code uses explicit nonnegative bounds. A newly created
analysis script was initially placed one directory too high, then moved
within the RH workspace before execution. No historical input changed.

## Next mathematical obligation

The inherited small first moment already makes the contribution from
|z|<=A sqrt(d) negligible on the logarithmic scale for every fixed A.
The surviving problem is the signed contribution of rare large-amplitude
source values, with the exact gate, Euler factor and all middle nodes.
The L2 bound cannot discard those tails. This finite low-cutoff test does
not resolve that obstruction, and more similar finite passes would not
replace an arithmetic theorem for it.

No provider calls, Forge experiments, Git publication or external sharing.
All owned workers are waited before the final checkpoint.
