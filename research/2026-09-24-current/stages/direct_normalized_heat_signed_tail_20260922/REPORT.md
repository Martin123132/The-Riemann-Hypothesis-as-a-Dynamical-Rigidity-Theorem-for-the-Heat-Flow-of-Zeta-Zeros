# Fixed normalized heat and the signed-tail boundary

22 September2026. TWO HANDS NETWORK LTD. Private additive RH research.

## What this stage establishes

The complete source can now be handled on a FIXED normalized heat interval,
with a uniform, paid cutoff error. The artificial radial shift, inverse,
and omitted lower heat endpoint are unnecessary. The cutoff error does
not require the inherited energy bound, a higher moment, a maximum source
amplitude, or uniform integrability of rare energy.

This does not prove the remaining favorable signed bound or RH. A separate
counterexample shows exactly why convergence at each positive heat time
cannot simply be integrated through time zero.

## 1. Uniform fixed-interval theorem

For Q_l(z)=|z|^2 exp(i l arg z), l=1,3,

    |Q_l(z)-Q_l(z) erf(sqrt(Lambda)|z|)| <1/(6 Lambda).

The sharp constant is0.16571661477885140466..., certified below1/6 by an
analytic uniqueness argument and two outward-rounded root enclosures.
This is a bound on the whole complex plane, including zero and arbitrarily
large values. Both angular modes receive the same radial factor.

The Gaussian integral gives a purely cubic representation with no inverse:

    Q_1(z) erf(sqrt(Lambda)|z|)
       =2/sqrt(pi) integral_0^sqrt(Lambda) |z|^2 z exp(-v^2|z|^2) dv,

with z^3 in place of |z|^2 z for the third mode.

Now normalize the actual source by

    V_N=U_N/(a(x)sqrt(d)),    dmu_N=(a(x)/4) dx dP_N,
    d=ceil(log(N+1)),         mu_N(total)=A_I<9.

This preserves the original observable exactly. Its truncated normalized
heat integral H_N(L) satisfies

                  |D_N-H_N(L)| <= C A_I/(6L),

uniformly in N and both original windows. Here C is the explicit fixed-S
multiplier bound. A fixed L gives a fixed uniform error; any L_N tending
to infinity gives an error tending to zero, however slowly it grows.
In physical heat coordinates the upper limit is L_N/[a(x)^2d].

For example L_N=log d SHRINKS that physical heat interval while giving an
O_I,S(1/log d) error. The stronger choice L_N=d^4 retains O_I,S(d^-4)
accuracy, now without needing an energy assumption for that payment.
At this latter cutoff, the rational scalar payment is one third of the
predecessor's sigma/2 regularization payment; all earlier results stay valid.

Original S, I, the norm atom, gate, global Haar subtraction and all remaining
source/core/exterior payments are unchanged. Fixed-S constants and practical
onset can be very large. This is not a numerical estimate for the original S.

## 2. The small-heat boundary is a signed energy tail

Let Z_N be the original combined real signed integrand expressed using V_N,
and retain its SIGNED amplitude tail

    B_N(R)=integral_(|V_N|>R) Z_N dmu_N.

Then exactly

    H_N(epsilon)=2sqrt(epsilon)/sqrt(pi)
       integral_0^infinity exp(-epsilon R^2) B_N(R) dR.

Thus small heat is a Gaussian average of the signed energy in large source
values. Negative tail contributions are retained; a bound on the probability
of large values alone cannot replace this quantity.

An exact abstract rare-event construction has globally centered multipliers,
zero first and nonconjugated second source moments, and unit quadratic
energy. Its damped cubic correlation tends to zero at every positive heat
time, yet its integrated signed response remains alpha+beta. Changing the
source sign gives the opposite response with identical energy laws.

This rules out that generic limit-interchange shortcut. It is NOT the
literal logarithmic RH source and says nothing adverse about the RH target.
Arithmetic cancellation in the actual B_N remains to be proved. Bounded
positive-heat correlations also still need an actual-source argument.

## 3. Whole-window finite ledger

The same N512,p257,x1 source and both original windows were frozen before
evaluation. All511 literal terms and all mixed terms remain. Each window
was integrated on65536 nodes at128 bits and131072 nodes at160 bits, with
explicit whole-window midpoint errors and complete append-only caches.

The normalized heat bands are

    [0,1/64], [1/64,1/4], [1/4,1], [1,4],
    [4,16], [16,64], [64,256], [256,infinity).

Eleven of the16 window/band contributions are certified positive; five
remain sign-unresolved. No negative sign is certified. In particular:

| Quantity | Window0 | Window1 |
| --- | --- | --- |
| Original response | [0.00041395553,0.00041699554] | [0.00032871593,0.00033173594] |
| Lowest normalized heat band | [0.00002735956,0.00003375957] | [0.00002182773,0.00002820774] |
| Direct strong-cutoff payment | <0.00000073297 | <0.00000073297 |
| Direct strong cutoff, error restored | [0.00041154559,0.00041940560] | [0.00032630592,0.00033414593] |

These are one-prime, one-profile controls, NOT the original fixed S and
full profile integral. They cannot be compared with the original core
threshold as if they evaluated it. Positive band values at N512 do not
prove or disprove a large-N limit. No favorable sign was an acceptance gate.

The finite point-profile normalization has mass a(1)^2, not A_I or a(1)/4;
this is tested explicitly. Both physical and normalized partitions sum
back to the original response without removing any band.

## 4. Verification

-177 named scalar/normalization/derivative checks passed.
-1343 integration and audit checks passed, including10 tamper guards.
-258 two-precision interval comparisons overlap.
-All384 source-cache blocks cover the complete windows.
-Both96-step root brackets replay, certifying the global rational1/6 cap.
-The original responses and energies agree with authenticated predecessors.

Arb supplies the source enclosures and sharp-constant proof enclosures.
Separate mpmath controls test Gaussian quadratures, derivatives and scalar
identities; they are not a second rigorous implementation of the full
source integration. The analytic arguments are documented, not formalized.

INPUTS authenticates8383 historical pins and all91 predecessor payloads.
The external seal, clean-replay and final-checkpoint receipts report the
actual packaging and fresh replay outcomes. This report does not anticipate
that later replay's success. No source-calculation failure occurred.

## Remaining target

The original threshold0.00092027989 and complete0.0011250999 budget remain.
The next source-specific work is to control the globally centered,
directional SIGNED energy tail B_N, together with the bounded positive-heat
correlation, or bound their combined H_N(L) directly. Neither independent
prime phases nor an unproved exchange of limits is available for free.

Status: DIRECT_NORMALIZED_HEAT_PROVED_SIGNED_ENERGY_TAIL_OPEN.
All work is local. No provider calls, Forge experiments, new zero counts,
GitHub update or external publication. RH remains open.
