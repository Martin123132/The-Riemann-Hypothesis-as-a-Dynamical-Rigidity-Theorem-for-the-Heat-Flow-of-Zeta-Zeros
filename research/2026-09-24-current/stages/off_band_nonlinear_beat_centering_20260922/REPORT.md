# Off-band nonlinear beat: a paid centering gain and a finite interaction test

22 September 2026. Private RH research, TWO HANDS NETWORK LTD.

## Result

One genuine arithmetic interaction is now controlled at its first nonzero
scale. It does not settle the complete-source bound or RH.

At the actual cutoffs N=k^5, the two literal rough indices k^2 and k^3
both survive the previously paid resonance-band removal. Nevertheless,
their nonlinear third mode satisfies the exact gate resonance

    4 phi_(k^2)-phi_(k^3)=2 chi+psi_N.

The uncentred leading response has a certified favourable sign in BOTH
original windows. On the complete original profile I=[1/8,2^24],

    0.0363141390 < integral rho |c(3/5)| c(2/5) dx < 0.0773141392,

so its scaled limit -beta/12 times that integral is enclosed in

    [-0.0005238609, -0.0002458608].

The component can pass through zero: 68 coefficient/tile incidences do so
in the interval calculation. The zero-safe Taylor proof retains them.

## The important extra cancellation

That raw response is NOT the objective. Applying the ORIGINAL full
prime multiplier and GLOBAL mean subtraction cancels the leading beat.
A quantitative Fourier argument gives, separately for modes 1 and 3,

    |integral (K_l-lambda_l) Q_l(U_pair) dnu|
       = O_I,S(1/[d rho_S k^3]).

Before centring the first nonzero scale is 1/[d rho_S k^(5/2)]. Thus the
centred pair gains a factor k^(-1/2). This is an unbounded fixed-pair
statement, with fixed S and both original windows, not a finite fit.
Its constants and onset can be very large. The proof uses a growing
finite Fourier truncation with its factorial tail paid before averaging;
it does not assume prime-phase independence.

The result is for the full multiplier. It cannot silently replace a
sharp-core restriction; the inherited exterior payment remains attached.

## Complete finite-source controls

We also kept the ENTIRE finite field, not just the pair, in the existing
one-prime, x=1 control. N=512 was frozen in the initial protocol. After
the user doubled temporary compute ceilings, N=1024 was frozen as a
separate stress test before evaluation. Both original windows were used.

The table gives outward enclosures for the retained-field centred
response and for the exact nonadditive remainder C(U)-C(V)-C(U-V), where
V is the old removed-band field. The prime kernel is complete.

| N | Window | Retained response | Mixed response |
|---:|---:|---:|---:|
| 512 | 0 | [0.0000153031, 0.0000173832] | [0.0000963867, 0.0001015868] |
| 512 | 1 | [-0.0000264650, -0.0000243849] | [0.0000715264, 0.0000767065] |
| 1024 | 0 | [0.0001407345, 0.0001462346] | [0.0001962205, 0.0002097006] |
| 1024 | 1 | [0.0001296236, 0.0001351437] | [0.0001356856, 0.0001492057] |

All four mixed responses are strictly positive. The retained response
changes sign between the N=512 windows but is positive in both N=1024
windows. There is no basis for extrapolating a favourable sign or
monotone improvement from the smaller control.

These controls use p=257 and one profile x=1. They do NOT use the
theoretical fixed S or integrate the full profile. Comparing these
numbers directly with the original RH budget would be invalid. The
previous removal theorem remains asymptotic, with no claimed useful
finite onset; these finite differences do not contradict it.

## Verification

* 1,368 explicit acceptance assertions passed, including contiguous
  coverage of every source cache block; these are not 1,368 theorems.
* 106 paired interval comparisons passed.
* 400 high-precision derivative controls and 1,200 zero-safe Taylor
  remainder controls passed, plus exact phase and Hessian checks.
* Full-profile quadrature: 14,336 covered tiles, paid analytic tail,
  Arb at 160/224 bits and separate mpmath.iv at 60/90 digits.
* Eight whole-window source builds: four window/cutoff combinations,
  each at 128/160 bits and two grid sizes. Their midpoint errors are
  paid over the entire windows. All use Arb, not independent libraries.
* A fresh-extraction replay is required by the sealing procedure and
  reported in the external CLEAN_REPLAY.json receipt. Do not infer its
  success from this report alone.

The analytic proofs are documented, not machine-formalized. One vacuous
descriptive marker was removed from the first exact-check script; the
old script/output and corrected acceptance are retained. No mathematical
input was altered to obtain a desired sign.

## Remaining obstacle and next action

The exact original target is still the growing retained-spectrum
centred response, with its .00092027989 core budget and inherited
payments. We now know that (i) off-band arithmetic beats exist,
(ii) one explicit beat has a favourable raw sign, (iii) the original
prime subtraction cancels it more strongly, and (iv) complete-source
mixed terms are not additive and can be adverse.

The next useful theorem must control those interactions collectively.
A pairwise Taylor expansion is not enough unless its remainder is
bounded uniformly after summation, including source zeros and regions
where no one term dominates. A fixed-pair Fourier approximation cannot
be promoted to an N-growing spectrum without that estimate.

Historical packages remain unchanged. The final integrity receipt and
bookmark record the checks. No provider call, Forge experiment, GitHub
update or external publication was made.
