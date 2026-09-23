# A universal paired-loss bound without spectral equidistribution

18 September 2026. TWO HANDS NETWORK LTD. Private local RH research.

## Result

The incoming weighted-harmonic reduction survives audit with explicit
endpoint and smoothing corrections. It establishes the actual weighted
central limit

    <X_N^2>_W / log N -> 1/2

in both original window families. The exact two derivative-moment limits
are still open. This stage shows that they are not needed for the following
upper bound on the normalized paired middle excursion:

    limsup_N max_window <P_pair(g_N)>_W / log N
        <= C_* < 0.076825.

This is an asymptotic upper bound, not equality for the actual loss and not
an effective finite-height certificate. The proof uses the inherited source
bridge and tail theorem, the written spectral-measure argument here, and the
new interval certificate. It is not proof-assistant formalization.

## Why the unresolved spectral distribution is no longer needed

Near-diagonal localization makes the limiting Gram form multiplicative.
Its positivity and central mass force every subsequential spectral limit
to have the form

    nu = Lebesgue measure on (1/3,1] + mu,
    mu >= 0, support(mu) subset [0,1/3], mass(mu)=1/3.

The same measure must be used throughout the transverse integral. It cannot
be chosen afresh at each transverse position. Among ALL such measures the
largest paired-variation majorant occurs at

    nu_0 = Lebesgue measure on (1/3,1] + (1/3)*delta_0.

That is a worst-case comparison measure, not a claim about the source's
actual distribution. A concave-tangent inequality, certified on all of
[0,1/3], proves the maximum. No equidistribution premise is inserted.

The exact interval result is enclosed outward by

    0.076824077322 < C_* < 0.076824127600.

The closed-interval tangent bound is J(y)<-0.69. In particular, any later
positive lower bound on the unknown measure's second moment would improve
the bound by at least (69/800)*integral y^2 dmu. That improvement is not
assumed here.

## The next obstacle

The unchanged central credit exceeds 0.0386971260539325. Thus a sufficient
remaining target for the SIGNED odd excursion is

    limsup_N max_window <A_odd(g_N)>_W / log N < 0.000569.

The more precise scalar allowance is greater than 0.000570124508, but the
displayed target conservatively uses the simpler paired cap 0.076825.
It is an upper bound on a signed average, not an absolute-value target.

The odd inequality has NOT been proved. The inherited tail theorem already
reduces it to fixed compact transverse intervals with paid endpoint and
anchor errors. Second moments and formal sign parity alone do not supply
the required arithmetic cancellation. Full gated positivity, zero-count
closure, and RH remain open. No effective onset is asserted.

The next work should retain the actual gate and same-height arithmetic in
that compact odd excursion. It need not first solve the stronger exact
derivative-moment equidistribution problem, nor repeat the tail reduction.

## Audit corrections

- Restoring a sharp interval gives an o(log N) product remainder, not the
  fixed-smooth-cutoff O(N^(-1/3)) rate.
- The hard-core formula retains its endpoint parameter; N tends to infinity
  before that parameter tends to zero.
- Small indices and far product pairs receive separate estimates. Their
  weights cannot all be bounded by the endpoint weight O(1/N).
- The unknown spectral measure is nonnegative, may be singular, and has
  fixed total mass. It is not silently replaced by uniform density.

## Verification and provenance

The 60- and 90-digit interval builds each covered 32 complete parameter
strips using 22,166 complete integration cells. Both include removable
endpoint values, positive-series remainders, differentiated remainders,
midpoint quadrature error and analytic infinite tails. Parameter intervals
are not point samples.

The independent rational ledger replay checked 22,203 overlaps and complete
coverage. All 17 algebra/structural test families and eight tamper-rejection
controls passed. Both numerical builds use mpmath.iv; they are not independent
interval-library implementations. See VERIFICATION.json for exact rationals.

ITERATIONS.md preserves the failed intake, invalid exploratory quadrature,
and one implementation failure. None contributes to the accepted certificate.
Fresh-extraction replay and final historical-pin receipts are stored beside
the sealed package, so completing them does not mutate the package.

No new auxiliary or zeta evaluation, provider call, Forge experiment, GitHub
update or external publication was performed. Incoming and historical files
are preserved byte-identically. The session bookmark is updated additively.
