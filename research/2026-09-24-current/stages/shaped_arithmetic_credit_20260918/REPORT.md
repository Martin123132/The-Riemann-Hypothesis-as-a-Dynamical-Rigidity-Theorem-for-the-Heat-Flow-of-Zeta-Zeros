# Stronger arithmetic credit, with the signed estimate still open

18 September 2026. Private local RH research, TWO HANDS NETWORK LTD.

**A new fixed witness raises the certified negative-central-energy coefficient
from greater than 0.038697 to greater than 0.04023.** The source, original
weight and both physical windows are unchanged. This is an asymptotic
actual-height theorem based on the retained transfer estimates, not an
extrapolation from a numerical height scan.

The corresponding sufficient signed-response budget increases:

| Quantity | Previous stage | This stage |
|---|---:|---:|
| Negative-energy lower coefficient | 0.0386971260... | 0.04023 |
| Universal unsigned upper cap | 0.076825 | unchanged |
| Sufficient signed-response upper target | 0.000569 | 0.00363 |

All figures use the same per-log-N normalization. The response allowance is
about 6.38 times the earlier conservative target because the earlier gap was
small. The energy improvement itself is about four percent, not sixfold.
**The response has not been shown to satisfy the new allowance.**

## What supplies the gain

The old divisor witness is multiplied by the fixed logarithmic weight

    w(u)=1+3u/4,  u=log n/log K,

and uses the same nonnegative construction, with fixed lambda=33/8. Its
positivity is exact at every physical height. Since 1<=w<=7/4, the previously
paid transfer estimates remain valid with fixed constants. No height-specific
fitting, favourable subset or new acquisition is involved.

The main square cost changes. It is evaluated using the genuinely correlated
four-variable arithmetic measure, not independent divisor splits. A finite
Chebyshev expansion and a positive spectral-tail bound give whole-integral
enclosures. The complete finite payoff is independently checked by exact
Laurent algebra; the two weighted lower-order terms are retained separately.

For this explicit fixed witness the limiting coefficient is enclosed by

    [0.0402391343422934, 0.0402551404547850].

The theorem uses only the conservative lower bound 0.04023. The interval is
for this witness coefficient, NOT an upper bound on actual negative energy.
No effective onset for finite N is asserted.

## The remaining target

With the exact normalized source g_N and X_N=g_N(0), it now suffices for the
same middle comparison to prove

    limsup_N max_window (1/(2 log N))
       <integral [X_N |g_N'| + |X_N| g_N'] dx>_W < 0.00363.

The moving-minimum reduction, fixed universal unsigned profile, origin and
tail estimates remain unchanged. The actual shared-height dependence and W
must still be used to control this signed nonlinear correlation. All other
source and outer-region obligations remain attached to the larger proof.

The stronger credit does not rescue a covariance-only argument. The prior
abstract comparison was adjusted to meet the NEW negative-energy floor; it
still permits a signed response above the sufficient budget. This is a
logical control, not an RH counterexample or an actual-source measurement.

## Verification and preservation

- 18 exact test families passed, including complete finite weighted payoff
  identities and the recovery of the unweighted identity.
- Nine altered records were rejected.
- The verifier recomputes all 4,096 final quadrature cells, including each
  whole-cell second-derivative enclosure. All 4,114 interval overlaps pass.
- Angular calculations use 90 and 110 digits; a rational half-angle route
  uses 90 digits. The routes share mpmath.iv and derivative arithmetic.
- The analytic argument is written, not machine-formalized. Final replay and
  historical-integrity receipts are delivered beside the sealed package.
- No source-height scan, provider call, Forge experiment, repository update
  or publication occurred. Historical packages remain unchanged.

The step improves a proved arithmetic bound. It does not close the signed
estimate, establish a full positive score, or prove RH.
