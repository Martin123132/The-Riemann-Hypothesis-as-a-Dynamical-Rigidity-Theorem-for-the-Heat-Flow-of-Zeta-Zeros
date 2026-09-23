# A concrete route through the remaining alignment problem

The divisor transform controls a numerator. The unresolved response also
contains its shared-height, source-dependent inverse. Treating them as
independent is invalid; generic positive-measure constraints do not fix
their signed interaction.

This stage supplies an explicit reference and an affordable finite test,
not a proof of the final asymptotic inequality.

## New analytic results

The source-independent Euler/phase reference Gamma_N has the positive limit

    gamma = 0.0352969104924041...
    Gamma_N = gamma + O(d^2 N^(-1/40)).

This is proved uniformly in both original windows using convergent
binomial expansions, an integer lower bound for nonzero character
frequencies, and integration by parts. The same sharp windows and Euler
primes remain. Multiplication by the retained band majorant still leaves
an o_I(1) reference-replacement payment. A qualitative limit alone would
not have justified that step.

The original nonlinear kernel admits the exact split

    Xi = L0 Gamma_N - Delta,
    |Delta| <= L0 M R_N(Ubar).

Here M is the original positive weight's mass, Ubar its mean of |C|^2,
and R_N is the explicitly retained concave radial attenuation function.
The bound is Jensen's inequality under that SAME measure. No independent
average, fourth moment or source nonvanishing condition is used.

Therefore a sufficient upper bound is

    Re(B Xi) <= L0 [Re(B Gamma_N)+|B| M R_N(Ubar)].

We can now test a signed arithmetic reference against a quantified
second-moment payment. There is also a valid first-moment alternative
for rare amplitudes; no asymptotic payment for it is silently assumed.

## Complete-window finite certificates

At x=1/2, full interval quadrature proves the following upper bounds for
Re(B_N Xi_N), including the inverse payment and all quadrature error:

| N | window centre | nonempty products | certified upper bound |
|---|---|---:|---:|
| 1024 | 1/4 | 3 | < -1.36e-7 |
| 1024 | 3/4 | 4 | < -7.20e-7 |
| 2048 | 1/4 | 26 | < -2.61e-8 |
| 2048 | 3/4 | 19 | < -1.09e-6 |

The numerator is the existing FROZEN middle-band numerator. Its original
nonlinear height kernel is bounded, not replaced. These do not certify
the sign of the pre-freezing exact-height band without its extra finite
phase payment. They are fixed x slices, not profile integrals or RH
zero-count certificates. J_out is not evaluated.

The exploratory roster was fixed at N=512,1024,2048, both windows,
x=1/200,1/2,2. N=512 is empty in both windows: its six slices are NOT
credited as successful sign tests. All twelve nonempty numerical slices
were negative, but only the four interval-certified slices above carry
rigorous sign claims.

## The next proof target

Keep J_out and show that the integrated expression

    J_out + beta L0 integral [gamma Re B_N
                    + |B_N| M_N R_N(Ubar_N)] dx/(4a)

fits the existing 0.0011250999-per-log-N budget, with every previous
payment and the original endpoint limit order retained.

This is a sufficient route, not an established estimate. It could be
too costly in some regions; there the signed defect Delta must be kept.
The coarse absolute divisor bound still gives O(d^(7/4)), too large for
an O(d) budget. We must obtain actual signed arithmetic or a sharper
modulus/attenuation estimate, while retaining cancellation with J_out.

The practical next step is to identify which profile regions can afford
the second- or first-moment payment, and reserve a signed treatment only
for those that cannot. Reformatting the seven divisor quantities alone
will not supply this missing sign.

## Verification and limits

- 2819 exact rational checks, including the Euler constant and coupled
  reference/defect controls; 2020 additional scalar first-moment checks.
- Independent complete-circle Euler integration at 128 and 192 bits;
  both enclosures contain the exact rational constant bounds.
- Four source certificates at both precisions, 147456 complete cells
  per precision; 36 actual-source field overlaps and six Euler overlaps.
- The two coarse unresolved range attempts and two engineering failures
  are retained. None is counted as a passing mathematical certificate.
- Both numerical implementations use FLINT/Arb; this is not independence
  of special-function libraries. Analytic proofs are not machine formalized.

Historical packages remain byte-identical. No provider calls, Forge
experiments, GitHub changes or external publication. Consult the delivery
receipts for the final fresh replay and integrity audit.
