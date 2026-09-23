# A favourable signed component, not yet the full bound

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
Private local RH research, 17 September 2026.

## Main result

We obtained an unbounded signed estimate for the **retained arithmetic
source**. For every fixed alpha>1/3, uniformly in alpha<=tau<=1/2 and both
original cutoff-window families,

    E_W[X_N F_N(4*tau^2)] = M_delta_N(tau)+o(1),

    M_delta(tau)=exp(delta*(1/2-tau))*cos(pi*tau/4)
                   * [zeta(1+tau)/2+1/12],
    delta_N=ceil(log(N+1))-log N.

The explicit uniform error order is

    O_alpha(N^(1/3-alpha)*log N+N^(-2/27)*log N).

M_delta is strictly decreasing and -M_delta'>1/10. Thus

    E_W[X_N*(Y_N-F_N(4*tau^2))]
       <= -(1/2-tau)/10+o(1).

This is eventually strictly negative at every fixed positive distance from
the endpoint. It covers all sufficiently large integer cutoffs, not merely
a finite scan or selected degree-edge subsequence. The direct proof does
not assume the sharper central mean-square bound. It groups actual paired
integer terms, retaining the near-stationary product hyperbola and corner.
The exact retained source and its stated reserve bound remain prerequisites.

The proof is in DIRECT_PAIR_PROOF.md. A wider extension down to tau=1/5 is
also derived in DERIVATION.md, but explicitly depends on the inherited
sharp central mean-square hypothesis. No effective numerical onset is
claimed for either asymptotic estimate. This is a written proof, not a
proof-assistant formalization or an independent external review.

For example, the worst degree-offset limiting gaps are below -0.36082 at
tau=3/8 and below -0.15704 at tau=7/16. These are asymptotic comparison
constants, not finite-N error allowances.

## Complete saved-window check

All 65,536 authenticated N2980 children were retained. At four fixed
transverse probes, the linear, unsigned-central and positive-gated gap
moments are all negative. The latter upper bounds, rounded conservatively,
are:

| tau | Upper bound for E_W[X_+*(Y-F_N(4*tau^2))] |
| --- | ---: |
| 1/5 | -2.79818 |
| 1/4 | -1.95678 |
| 3/8 | -0.66834 |
| 7/16 | -0.27688 |

The endpoint tau=1/2 contributes exactly zero. The first two rows are
finite checks; the direct unbounded theorem does not cover those two points.

Crucially, choosing the largest gap at each height from these four points
and the endpoint instead gives

    0.0017220 < E_W[X_+ max_j(Y-F_N(4*tau_j^2))] < 0.0030020.

So **negative fixed-point means do not remove adaptive selection cost**, even
on this actual source window. This is a finite-roster cost, not an upper bound
for the whole transverse minimum. It is fully compatible with the existing
negative complete-window score, which also survives this replay.

## What is still missing

For each fixed gap V_tau=Y-F_N(4*tau^2), exactly

    E[X_+ V_tau]=(E[X V_tau]+E[|X| V_tau])/2.

We now control the first term favourably on the outer band. The unsigned
correlation remains, followed by the adaptive competition and the exact
middle-region loss. Those cannot be dropped, inferred from the linear mean,
or treated as independent noise.

The full requested inequality

    J_N < 2*mu_edge-2*E_remaining,N

is **not proved**. No whole-strip zero exclusion, absolute zero-count result
or proof of RH follows from this stage. The next mathematical obligation is
to control that same-height unsigned correlation and the adaptive excess,
not to rerationalize the already controlled approximation errors.

## Verification and preservation

* 42 named tests passed, including exact resonance enumeration and guarded
  counterexamples to invalid gating, selection and credit-discarding steps.
* 1,053 exact scalar states and 3,000 interval-containment trials passed.
* Complete source reaggregation at 30- and 42-digit fixed-point scales gave
  180,224 overlapping parent-level interval fields.
* Limiting constants at 160 and 224 Arb bits gave 42 overlapping fields;
  an elementary positive-series/integral bracket independently checked each
  real zeta constant. Both numerical routes use Arb for elementary powers.
* No new Dirichlet-source sum, height scan, zero count, provider execution,
  Forge experiment or external publication was performed.
* All 2,164 previously pinned files were authenticated before the work.
  The final audit and clean-extraction receipt are retained with the delivery.

The numerical checks support implementation correctness. The universal
conclusion comes from the written paired-sum argument, not extrapolation
from the finite data. Every historical sealed package remains unedited.
