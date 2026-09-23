# Signed pair audit: the proposed integration by parts is circular

23 September 2026. Private TWO HANDS NETWORK LTD research.

The user-supplied paired-sum formulas for H and K are correct. K is
independent of the scalar drift q and equals -(d/dt)|V|^2/(2d). The supplied
real-part directional identity is also correct for the expression it states.

For the actual production formula, restore the five-channel gate coefficients
and take the imaginary part. With B and C defined in PAIR_TRANSPORT.md,

    T=H C-2K B, P=|V|^2 B.

The proposed integration of the K term has now been carried out completely,
retaining the exact Jacobian, endpoints, gate derivative, multiplier derivative
and source-angle derivative. It gives

    I_K=B_boundary-Dhat-I_H-R_lower_order.

Adding I_H recovers precisely the original signed-transport identity. The
original unknown Dhat returns with coefficient ONE, not a small multiplier.
The source-angle derivative dH/|V|^2 cancels the H term; the remaining gate
rotation restores the original response. There is no new saving from this
complete integration-by-parts step.

This is a useful exclusion of a circular route, NOT a favourable bound and
not evidence that the desired signed inequality is false. The open target is
still the actual same-height weighted pair correlation. Its angular weights
depend on the same source; the quadratic energy laws alone cannot estimate
that dependence. No new optimized residual norm or independent-phase model
is substituted.

## Verification performed here

- All three last sealed packages authenticated: 121 + 122 + 61 = 304 payloads.
- All 9,758 pre-existing bookmark pins authenticated before work.
- 4,201 exact rational assertions over 200 paired-source fixtures and 200
  complete transport fixtures; 398 nontrivial omitted-term/convention guards.
- The same 60 archived finite source points replayed at 192 and 256 bits.
  These are not new source evaluations, windows, zero counts or quadratures.
- Full verbatim Joint Mode derivation and audit, plus Signed Transport
  derivation and audit, assembled in HANDOFF_FULL.md with individual hashes.

The finite identities do not prove an asymptotic sign. Existing arithmetic
interfaces are inherited, and the full RH proof chain is not audited anew.
The .00092027989 target and .0011250999 full ledger are unchanged. RH remains
open. Final replay and preservation receipts are external to the sealed
package, and the bookmark is updated only after they pass.

No provider calls, subagents, Forge experiments, external publication, GitHub
update or new RH source scan occurred. Only additive local evidence and the
research bookmark were written. All owned workers are awaited and closed.
