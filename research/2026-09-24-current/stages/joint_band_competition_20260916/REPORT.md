# The joint selector passes; the real competition is retained

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
Private local research, 16 September 2026. No publication performed.

## Main conclusion

All three supplied packages authenticate (75/75 members), and all 91 incoming
named tests pass. I reviewed the collar, moving-tail and ladder derivations,
including their inherited mean-square dependency. No missing middle region or
unpaid replacement of the complete gate was introduced.

The new step handles the competing band values at the SAME height using a
paid soft minimum. With beta=d^2 and B=O(log d) bands, the comparison error is

    O(log log d / d^(3/2)) = o(d^-1),

under the inherited central L2 bound. This is smaller than the existing collar
payment. It neither assumes a winning band nor discards ambiguous switches.
It does NOT bound the actual arithmetic competition by that small error.

## Finite results

Both complete original N=64 windows, all 2048 height cells, and the unchanged
four-band calibration are retained. Beta=25 was frozen before the calculation.

| Lower certificate | First window | Second window |
| --- | ---: | ---: |
| Joint soft minimum, before tail restoration | >0.8400762 | >0.8887137 |
| After one-sided adverse-tail restoration | >0.8043304 | >0.8544546 |
| After absolute-tail restoration | >0.3426597 | >0.5653413 |

The source is unchanged. These are slightly weaker than the original hard
minimum certificates, as expected for a rigorously paid lower comparison.
No new heights, zeros or physical margins are claimed.

The observed whole-cell gaps bound the artificial smoothing cost by less than
0.003610 and 0.002589, respectively. The generic log(4)/25 bound would pay up
to 0.045261 and 0.042069 after the gate average. The smaller bounds keep the
actual simultaneous band geometry rather than assuming four independent ties.

## The bands really do compete

On the first window, each of the four hybrid bands is the strict winner on
some whole cells with a strictly active gate:

| Band | Certified whole cells |
| --- | ---: |
| Full inner region | 87 |
| K=32 band | 5 |
| K=16 band | 2 |
| K=8 band | 228 |

These are counts of certified cells, not counts of switches. Ambiguous cells
are retained. In the second window K=32 has no certified strict-winning cell;
that is not proof it never wins.

Using the outermost shortened-band minimum b_e as reference, define

    C(r)=max(0,max_(j != e)(b_e-b_j)).

Then min_j b_j=b_e-C exactly. The normalized genuine competition cost is

    first window:  [0.01951971, 0.03313684],
    second window: [0.01571961, 0.03217518].

Its strictly positive lower bounds certify that the penalty cannot simply be
deleted. This is a different quantity from the much smaller entropy allowance.
Even paying this competition independently against the outer reference leaves
a positive finite hybrid score. All remaining finite tail restoration stays
separate and is included in the table above.

## What remains open

The next arithmetic theorem must bound the ACTUAL joint competition C against
its short outer reference on an unbounded family, or directly lower-bound the
Gibbs-weighted arithmetic expression with its entropy term retained. Positive
averages for individual bands are insufficient. The exact middle region is
still present inside C, and intra-band minimization has not been removed.

The new uniform result concerns only the accuracy of a comparison. It does
not prove the unbounded sign. Winding, multiplicities and absolute-count
closure remain additional obligations; no RH proof is claimed.

## Verification

An independently implemented de Casteljau restriction and directed-integer
reaggregation completed all 2048 cells, filling the previous ladder replay
gap without restarting its timed-out implementation. Its hard-minimum lower
floors agree with the original acceptance at 20 decimal places, with only
conservative outward rounding differences.

The new suite passes 36/36 tests, including strict dominance, all possible
bands, source restoration, exact coverage, tied/negative controls, the Gibbs
derivative, and a counterexample to inference from positive marginal means.
Arb soft evaluations run at 160 and 224 bits; exact-Fraction integration is
separate. The accompanying independent-soft record describes an additional
mpmath.iv replay of every vector. All share the authenticated input source
enclosures; the 62-term source sums were not recomputed.

The inherited Watt attribution was checked in the primary paper:
[Bourgain and Watt, page 1](https://arxiv.org/pdf/1505.04161).
The new analytic reasoning is written mathematics, not formal proof-assistant
verification. Derivations, logs, failed wrapper/test iterations and source
hashes are retained. See the separate clean-replay receipt for the extraction
and byte-reproduction checks actually completed.
