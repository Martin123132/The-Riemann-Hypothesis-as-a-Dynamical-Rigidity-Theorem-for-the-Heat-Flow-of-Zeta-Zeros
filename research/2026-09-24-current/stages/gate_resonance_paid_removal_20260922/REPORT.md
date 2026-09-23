# Gate transport audited; the single-term resonant bands are now paid

22 September 2026. Private TWO HANDS NETWORK LTD research.

The incoming companion and gate-cell arguments reconcile, with one external
wording correction: their y_n frequencies are frozen in height. Their displayed
derivative equations already use that convention. All historical files remain
unchanged, and all four accepted checker outputs reproduced byte-for-byte.

## New complete-source result

Let A=N+(1+2nu)/4. Remove the terms n=1,2 and the three bands

    |n-A^(j/3)| <= sqrt(A^(j/3)), j=1,2,3,

intersected with the original cutoff and S-rough indices. These deterministic
bands include the single-term gate resonance locations. Across both original
windows their combined reciprocal square-root mass is less than 15.

The removed field has energy at most 225*C_I/(rho_S*d). Its effect on the
COMPLETE nonlinear centered response, including interaction with everything
left behind, is bounded by

    (exp(162)+1)(2alpha+3beta)
      * [15 sqrt(C_I*E_U/(rho_S*d)) + 225*C_I/(2rho_S*d)].

Since E_U=O_I,S(1), this is O_I,S(d^(-1/2))=o(1). The original gate,
both windows, global centering, original profile, and all fixed payments remain.
The sufficient target stays exactly 0.00092027989. There is no new numerical
reserve and no claimed useful finite onset for the potentially large constants.

## Why the remaining sign still needs arithmetic

The profile/transport identities alone admit both signs of a third-harmonic
response, even with zero quadratic phase average. This is verified across the
full original profile interval, not merely at one point. A common positivity
shortcut also fails: the actual combined response Hessian has a negative
tangential direction on the positive real axis.

Actual arithmetic gives a complementary positive result. Along N=k^3,
fixed finite rough clusters around n=k and n=k^2 have an explicitly locked
third harmonic, but their globally centered selected-prime response vanishes
asymptotically. This handles the entire finite cluster nonlinearly. It does
not justify treating a growing source as a fixed cluster.

## What remains open

After the paid deletion, many remaining terms can still interact to create
nonlinear beat frequencies. A favorable bound for that complete centered
response is still missing. Nothing here proves RH or lowers the numerical
threshold. The next target is this growing-spectrum signed interaction, not
another attempt to delete fast gate derivatives or assume Hessian positivity.

## Evidence

See VERIFICATION_FINAL_V2.json for the exact test totals. The incoming four
implementations reproduce 3,408 assertions. The full-profile control has
14,336 closed tiles per build and a paid analytic tail to 2^24, using both
FLINT/Arb and a separate mpmath.iv implementation. Arithmetic carrier checks
cover both windows at six cubic cutoffs; complete-source band guards include
sixteen outward finite rosters. The isolated integer-conversion failure is
retained and repaired in a successor script, not concealed. A verifier-receipt
handling error is also retained; its correction changes no mathematical result.

No original auxiliary integration, source-zero count, provider call, Forge
experiment, GitHub update or external publication was performed.
