# Cubic correction retained; fourth-moment remainder controlled

20 September 2026. TWO HANDS NETWORK LTD. Private RH research.

## Incoming audit

The new factor-allocation identities, pointwise linear cancellation,
variance interface and opposite collected coefficient signs are accepted
with the supplied clarification. All three original checkers replay
unchanged:528985,2997 and64133 exact assertions. All27 incoming files,
62 predecessor payloads and5208 historical pins authenticate.

## New coefficient theorem

Let Q(P) be the allocation variance, theta the EXACT admitted fraction,
V_adm its second moment, and M3_adm its SIGNED third moment, always
normalized by the full ordered factorization count tau_3(P).

The coefficient now has the exact form

    a_P(x)/(tau_3(P)x^3)
      = theta h0^3+(3/2)V_adm K2+M3_adm K3+R4_P(x),
    |R4_P(x)| <=(32/9)W4_adm <=64Q(P)^2.

This holds for0<x<=1 and ybar in[0,1]. K2,K3 and the actual admitted
fourth moment are defined in DERIVATION.md. The cubic contribution is
retained, not charged a guessed sign or silently removed.

The elementary reason is that every prime-exponent allocation has a
negative fourth cumulant. Exact finite-count convolution gives

    Av_all (delta_1^2+delta_2^2+delta_3^2)^2 <=18Q(P)^2.

No physical phase-independence model, Shiu extension or numerical
source-height assumption is used in these new arguments.

Under the SAME incoming variance and cutoff conditions, the sufficient
negative-coefficient threshold improves from65536 to1024 distinct primes:

    r>=1024, Q<=1/r, delta_P>=1/4, ybar in[0,7/20]
      ==> a_P(x)<-tau_3(P)x^3/1024 for0<x<=1/100.

This is a finite coefficient criterion, NOT a source-height threshold,
mass estimate or claim of favourable full-band response.

## A sign trap resolved

For the positive middle-prime family supplied in the incoming package,
the full and admitted third moments have opposite limits:

    full factorization mean of S3 -> -8/45;
    actual admitted M3_adm        -> +1/180.

The latter follows from the complete triangle, the already paid roster
difference and exact polynomial integrals, not extrapolation of samples.
The admitted variance tends to1/72 and fourth moment to1/90.
Four separately counted canonical prime-power rosters also exhibit the
sign reversal. They are diagnostic LOW-prime cases, not middle-class
source samples and not evidence for the limiting theorem by themselves.

Thus the unmasked negative skewness cannot be assigned to the source
coefficient after the actual cutoffs. This matters because K3 is
nonnegative in the small-profile region.

## Remaining target

The fourth-order response error is bounded by

    64 beta integral_I x^3 |Xi_N(x)| dx/(4a(x))
                *sum_P tau_3(P)Q(P)^2/sqrt(P),

or the sharper exact-admitted-fourth-moment version. This total is NOT
proved o_I(log N). The count, variance and admitted-skew terms remain
coupled to the same complex quadratic phase, actual Xi_N and J_out.
No new band is deleted, and the original limit order is unchanged.
The favourable joint signed bound0.0011250999 and RH remain OPEN.

## Verification and preservation

- New primary:13800 exact assertions across35 named families.
- New independent:14548 assertions across25 named families.
- Four canonical rosters, including all1848 admitted ordered triples for
  the L=64 control;20 interval fields overlap at112 and176bits.
- Formal polynomial identities in the prime exponent, independently
  convolved marginal moments, exact triangle integrals and roster census.
- One failed numerical wrapper conversion retained and corrected without
  changing the mathematical input. An output-name collision on a later
  rerun was also retained; no old output was overwritten. The final checker refinement
  replaces a trivial sanity check with an explicit positive-variance,
  zero-error control; both executed versions remain available.

The analytic proofs are not machine-formalized. Two interval precisions
share FLINT/Arb. Acceptance and clean-extraction receipts distinguish
new checks from replayed incoming checks and record exact output hashes.

No provider calls, Forge study, GitHub update, source-height scan, zero
count or external publication. Historical packages remain unchanged.
