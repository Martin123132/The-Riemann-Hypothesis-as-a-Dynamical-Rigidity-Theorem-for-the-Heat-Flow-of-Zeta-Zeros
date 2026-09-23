# Euler phase peeling: an explicit allowance, an unresolved correlation

19 September 2026. Private local RH research. Historical evidence unchanged.

## Result

We derived a fixed-prime factorization of the actual retained arithmetic
sum with an O_S(1) L2 error, hence negligible on its sqrt(log N) scale.
The same replacement is stable for every retained nonlinear angular mode,
including at source zeros. This uses the authenticated ratio-Gram argument;
the asymptotic deductions are written proofs, not machine formalization.

The energy-weighted phase supplied by each removed Euler factor has an
explicit convergent series. Its k-th coefficient has product decay
O_k((log(P+1))^(-k^2/4)) when the fixed prime set is enlarged AFTER the
height limit. The first mode decays slowly. The higher odd modes decay
much faster, so the useful version keeps the actual first mode exactly.

**This does not prove phase independence for the remaining sum.** We retain
its signed nonlinear correlation as a separate term. An explicit abstract
counter-control shows that quadratic decorrelation cannot remove that term.

## A small fixed-set target

Let M_N be the existing first nonlinear moment and D_high,S,N the signed
coupling of the higher odd modes to the removed Euler factors, defined in
DERIVATION.md. The original gate, original two windows, all71 physical
indices, and real-source endpoint payments remain attached.

For S={2,3,5}, the explicit higher-mode allowance is below0.000517086888.
The following bound would suffice for the remaining middle budget:

    limsup max_window [(4/(3*pi))*M_N + D_high,{2,3,5},N]/log N
        < 0.002417.

For the larger, still fixed set of primes through997, the explicit allowance
falls below0.000031665379 and the sufficient right side becomes0.002902.
Neither joint bound is proved. A claim that M_N alone suffices would omit
the unresolved coupling. No effective N threshold is available, and the
large fixed-set constants need not be computationally practical.

All seven preselected cutoffs are reported, including the negative one:

| Prime cutoff | Known higher-mode allowance | Remaining joint allowance |
| --- | ---: | ---: |
| 2 | 0.004581917821 | negative; insufficient |
| 3 | 0.001307065929 | >0.001627285031 |
| 5 | 0.000517086888 | >0.002417264072 |
| 11 | 0.000272365483 | >0.002661985478 |
| 31 | 0.000133287449 | >0.002801063511 |
| 97 | 0.000077375879 | >0.002856975081 |
| 997 | 0.000031665379 | >0.002902685581 |

The first column of allowances bounds the factorized, explicitly known
part, NOT all higher nonlinear modes. The unknown signed coupling can
reinforce or offset it. No independent worst-case replacement of that
coupling has been silently made.

## Why the obstacle has not disappeared

The abstract counter-control has flat conditional first and quadratic
moments and the same fixed-prime energy law, but its cofactor cancels the
Euler phase. It has first nonlinear moment zero and signed half-square1/12.
Thus even first-mode cancellation plus our new quadratic facts do not
guarantee the target. This is a test of logical sufficiency, not an RH
source or a counterexample to RH.

The next mathematical target is the actual signed coupling, initially with
S={2,3,5}. The third angular mode is the first unresolved higher mode;
any estimate must keep its same-height relation to the cofactor and original
physical phase. A generic independence assumption or a variance fit cannot
replace that arithmetic input.

## Verification

-7 exact test families,962 exact checks.
-112- and192-bit series builds:5845 corresponding interval fields overlap.
-8 whole-angle phase integrals and4 Poisson-mass integrals cross-check the
 series through two formulas. Both use Arb, not independent libraries.
-14 exact rational comparisons recheck all seven first-mode budget rows,
 using a separate Machin-series enclosure of pi.
-13 guards reject unsupported strengthened conclusions.

No new Dirichlet-source evaluation, height scan, zero count, provider call,
Forge experiment, external publication or repository push was performed.
The numerical tests validate algebra and scalar certificates, not an RH
proof or a numerical proof of the asymptotic analytic argument.

Clean extraction, replay and integrity receipts are recorded beside the
sealed package after it is created. They are not asserted solely on the
basis of this document. See AUDIT.md for the retained failed iteration.
