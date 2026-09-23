# Signed transport and the limit of phase-drift optimization

23 September 2026. Private TWO HANDS NETWORK LTD research.

## Result

The remaining arithmetic current has an exact transport representation with
an arbitrary fixed nonresonant phase drift. Its optimal least-squares drift
is determined by the original spectral moments, rather than fitted to target
outcomes:

    q*(x)=G1(x)/G0(x),
    residual energy=G2(x)-G1(x)^2/G0(x).

The remainder is quadratically orthogonal to the original source at this
choice, but still has positive energy. Orthogonality is not angular
independence and supplies no favourable sign on its own.

The next shortcut has now been tested and ruled out: choose any scalar drift,
then pay the signed remainder by the defined modewise Cauchy-Schwarz bound.
It fails even if all five gate channels are combined before estimating.

    norm certificate > 0.0023
    remaining signed target = 0.00092027989
    complete original ledger = 0.0011250999

The floor uses only x in [1,2] of the unchanged profile interval. It holds
for every admissible drift, not just the original 1/2 or the least-squares
choice. Its proof retains the constant gate coefficient and the irreducible
spectral spread. Nothing here says the actual current exceeds 0.0023.

## What is useful

The exact signed current keeps the common source, its logarithmic companion,
the prime multiplier's derivative, the sharp endpoints and the r-to-t
Jacobian. The endpoint and lower-order payments vanish asymptotically at
the inherited arithmetic interfaces. Any two fixed nonresonant choices give
the same signed current to o(1), though their size bounds can differ.

The next estimate therefore needs the SIGN of the source/residual angular
correlation. Further tuning of the scalar phase rotation cannot rescue the
tested size-only bound. This narrows the method choice; it does not close
the signed theorem or prove RH.

## Verification

- 9,009 exact arithmetic assertions across 250 moment and 150 signed-transport
  fixtures, including the Jacobian sign, resonances and zero-safe derivatives.
- 144 further exact assertions for 48 grouped-gate fixtures.
- Complete profile integration at 160 and 224 bits, with 2,048 x cells and
  4,096 y cells in each, retaining the original norm atom.
- Independent rational reaggregation certifies the floor; 16,384 paired
  mathematical interval fields overlap. Nine malformed/tampered cases reject.
- Separate 70- and 90-digit mpmath diagnostics at 17 fixed points confirm the
  original norm formula and moment enclosures. These diagnostics are not
  interval proofs.

All routine failed iterations are retained in AUDIT.md and their receipts.
The written analytic proof uses inherited arithmetic interfaces; the tests
do not prove those interfaces or formalize the RH dependency chain.

MANIFEST.json and the external SEAL_RECEIPT.json bind this package. The
external CLEAN_REPLAY.json records actual fresh-extraction replay results.
The additive bookmark is updated only after those checks finish.

No original RH source-height computation, provider call, Forge experiment,
publication or GitHub update occurred. Historical packages remain immutable.
