# Exact admitted pair transform

The two incoming batches reconcile. The analytic tail is paid at fixed
profile compactness, and the all-order fluctuation is bounded on the
stated small-profile domain. Neither supplies the missing joint sign.

## New result

After fixing P and its smallest factor n, all remaining pair weights can
be reconstructed exactly from a single positive divisor-ratio transform

    T(lambda) = sum_admitted_pairs cosh(lambda log(w/v)/(2d)).

The complete profile requires only seven quantities:
T(0), T''(0), T(x), T'(x), T(2x), T'(2x), T''(2x).
Equation (2) of DERIVATION.md gives the formula. It holds for every x>0
and introduces no approximation or transport error.

This is a fixed number of quantities PER minimum-factor bin, not a fixed
number for the entire problem. It preserves every admitted pair; it
does not make their arithmetic distribution free or multiplicative.

The transform has a positive shared measure. If Q is its count and U its
maximum allowed log-ratio, its two samples obey

    2 T(x)^2/Q - Q <= T(2x)
      <= 2(cosh(xU)+1)T(x) - (2cosh(xU)+1)Q.

There are also explicit second-moment and derivative bounds. Thus the
seven quantities cannot be assigned unrelated worst cases. On the
incoming small-profile domain the full pair weight is positive, so the
minimum factor alone determines each bin's coefficient sign.

## Still open

The signed coefficient then multiplies the unchanged oscillatory kernel.
Its favourable bound is NOT proved here. The next target is the joint
weighted divisor response, using the shared-transform constraints rather
than separately bounding every slot. The original Xi, J_out, source
payments and target 0.0011250999 remain unchanged.

This gives an exact alternative to the growing Taylor-moment list. It
does not remove the arithmetic difficulty of summing over products and
minimum factors. No gain in runtime or RH coverage is claimed.

## Checks

- All four incoming checkers: 77310 assertions, saved outputs byte-identical.
- New primary: 2010 exact rational checks, including 516 hyperbolic
  pair-stencil controls, full permutation regrouping and unchanged-kernel
  accounting on 12 formal prime-power rosters.
- Independent: 814 exact checks, including a symbolic Laurent identity,
  positive-measure inequalities and restricted-divisor reconstruction.
- Actual norm/log-coordinate controls: 1451 checks at each of 128 and
  192 bits; all 83 exported interval fields overlap.
- Four genuine integer diagnostic rosters have 12, 78, 408 and 1848
  ordered triples. They have one distinct prime and are NOT retained
  middle-class examples. Their role is implementation checking only.

See acceptance and delivery receipts for fresh-extraction replay and
final integrity results. Two pre-test replay-wrapper failures are retained.
Analytic proofs are documented, not machine formalized. Historical evidence
is immutable. No provider, Forge, GitHub or publication action occurred.
