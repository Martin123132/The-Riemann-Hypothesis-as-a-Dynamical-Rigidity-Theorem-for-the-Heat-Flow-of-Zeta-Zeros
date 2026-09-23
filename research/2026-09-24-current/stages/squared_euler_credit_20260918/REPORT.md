# Logarithmic central credit; the middle competition remains open

18 September 2026. Private research, TWO HANDS NETWORK LTD.

## Result

The available negative-central energy is now bounded below at logarithmic
order on both original shared-height windows:

    liminf_(N->infinity) min_(nu=0,1) E_W[(X_(N,nu))_-^2]/log N >= 1/288.

Together with the authenticated upper energy estimate, this puts that energy
at order log N. The previous sealed lower rate was order log log N. This is
an analytic asymptotic result, not a finite onset or a numerical height scan.

The improvement comes from an explicitly nonnegative, FINITE squared witness.
A short half-divisor polynomial Q_K, with K=floor(N^(1/48)), obeys
|Q_K^2|=|Q_K|^2. Thus

    h_K=(|Q_K|^2-Re(e^(i chi)Q_K^2))/3 >=0

has a finite arithmetic expansion, as does its square. The exact dual
inequality x_-^2>=-2xh-h^2 converts its paid mixed moments into the energy
bound. All phases remain at the same actual physical height. A finite Euler
product is used only to majorize a positive coefficient sum; no growing
infinite product or independent-prime model is substituted for the source.

The leading credit is at least log K/6. The unwanted resonance is
O((log K)^(3/4)), and every nonconstant physical-height transfer error tends
to zero. The proof is in DERIVATION.md, Sections 2-6.

## Effect on the remaining target

The stronger credit pays the following two regions at a negligible relative
cost, while retaining the original source, gate and exact minimum:

* Central collar: 0<=tau<=1/(2d log d), where d=ceil(log(N+1)).
* Outer band: b_N<=tau<=1/2, where eventually
  b_N=1/6+2log d/log N.

The collar uses a new completed-square bound for its CLIPPED competing loss.
It does not reuse the old signed-collapse estimate outside its allowed width.
The outer comparison now needs only a bounded L2 error: its gated cost is
O(sqrt(log N)), negligible against the new logarithmic credit.

The exact middle loss B_M remains. A sufficient next theorem is

    limsup_(N->infinity) max_(nu=0,1) B_M/log N < 1/288.

In particular B_M=o(log N) would suffice for the selected linear component.
Neither bound is proved in this package. The narrower unresolved interval is

    1/(2d log d) <= tau <= 1/6+2log d/log N

for sufficiently large N. Overlapping losses are retained as a maximum, not
added as an exact identity. The full positive-gated score, count closure and
RH remain open. A favourable selected linear component alone would not prove
those further statements.

## Verification and limits

* 17 exact-arithmetic test families passed, including 30 raw-character moment
  comparisons, the half-divisor convolution, signed-resonance controls,
  20,000 collar-square controls and 1,125 exact three-region ledgers.
* Separately coded 50- and 80-digit interval constructions overlap in all
  114 compared fields. They agree with 48 exact rational moment values.
* Both interval builds use mpmath.iv. Their coefficient recursion is different
  from the factorization-based exact implementation; they are not independent
  interval libraries.
* Small finite controls do not assert positive credit at every cutoff. The
  test suite retains the negative witness payoff -1/12 at K=1.
* The analytic inequalities are written proofs, not proof-assistant formalized.
  Tests check finite identities and controls, not an infinite claim by sampling.
* There are no new auxiliary-source evaluations, zero counts, provider calls,
  Forge experiments, repository publications or external transmissions.

The intake authenticated 2,572 historical artifact pins, 1,153 prior checkpoint
records, all 79 members of the immediate predecessor package and seven copied
analytic inputs. Historical packages remain unchanged. The delivery receipt
records the new manifest and ZIP hashes. A separate clean-replay receipt records
whether fresh extraction reproduces all four mathematical result files exactly.

## Next research action

Attack the collected middle loss with the common-height dependence intact.
In particular, look for a signed or gate-sensitive estimate rather than paying
the full middle supremum independently of X_+. Its natural sufficient scale
is now logarithmic. The new 1/288 budget is a conservative asymptotic lower
bound, not an observed finite margin and not an assigned value of A_N.
