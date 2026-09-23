# The smooth fixed-gap moment test passes both complete windows

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.

The continuum minima inside the four bands have been replaced by explicit
smooth LOWER comparisons, constructed from the same shortened arithmetic.
This avoids an unpaid sampling error and avoids transferring squared moments
using only an average absolute-error bound.

All 2048 original height cells, the exact inner region, the original gate and
the original height weight remain. The construction is frozen at eight panels
per band, 48 ratio controls per band and beta=25.

| Restored moment certificate | First window | Second window |
|---|---:|---:|
| Smooth-reference minus moment cap | >0.37804239 | >0.45075924 |
| One-sided tail payment | <0.10161704 | <0.15004025 |
| Final retained-polynomial lower score | **>0.27642536** | **>0.30071899** |

All six signed gap means remain strictly negative. The more expensive absolute
tail payment does not certify the first window; it does certify the second.
No parameter was retuned after seeing these outcomes.

## Why the step is useful

For each band let Z_j be the sum of exponentials of its certified Bernstein
ratio controls and ell_j=-log(Z_j)/beta. Positive denominator controls prove
ell_j is below the TRUE continuum band minimum. Therefore the first and second
moments of the fixed differences

    ell_outer-ell_j = log(Z_j/Z_outer)/beta

directly give a sufficient lower score. There are no hard intra-band minima
left in this sufficient arithmetic expression. The inherited ladder payment
is applied to the final score, not incorrectly inserted into squared moments.

## Audit and verification

- 1218 pre-existing historical artifact pins matched at intake.
- Both incoming packages: all 19 manifest-listed members verified; all 17
  supplied test families passed; full-ledger mathematical results reproduced.
- The record archive includes one unlisted bytecode file. It is retained in
  the original ZIP, recorded, and excluded from execution. Its analyzer changes
  metadata/field names relative to the saved result; all saved mathematical
  fields match under the documented mapping. No byte-identical replay claim.
- The supplied moment analyzer outputs the centered cap only. The optimized
  caps have additionally been reconstructed independently from its raw moments.
- 21 new test families pass, including a full exact-Fraction moment replay.
- 160/224-bit Arb builds agree exactly after outward rational export.
- Independent 70-digit mpmath.iv replay agrees on all 8192 band enclosures.

See the external clean-replay receipt for extraction and reproduction of the
sealed package. All original ZIPs, source cells and historical results remain
unchanged. No new Dirichlet evaluation, height scan, provider call, Forge study,
zero count, repository push or external publication was performed.

## Next mathematical obligation

Prove the signed log-partition-ratio moments and the smooth outer reference
meet the paid inequality on an unbounded cutoff family. The finite negative
means do not prove that. Control-hull conservatism is distinct from entropy
and cannot be claimed to vanish merely because beta grows. The current result
is a fully paid finite reduction, not RH or an all-height positivity theorem.
