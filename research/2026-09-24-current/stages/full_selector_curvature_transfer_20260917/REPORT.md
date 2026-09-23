# Full-selector angular tail controlled on the actual height curve

Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
17 September 2026. Private local RH research; no publication authorized.

## Main result

The full selection gap, including the central candidate and exact middle
region, becomes a convex support function after adding one explicit scalar
from the retained reserve. That scalar is kept exactly, not discarded.

Positive angular curvature then gives the coefficient bound

    |g_k| <= (g0+eta_N)/(k^2-1), |k|>=2.

Consequently, keeping angular harmonics through Q leaves the whole-angle
error at most

    (g0+eta_N)*(1/Q+1/(Q+1)).

This covers every angle, hence the actual shared-height path, without
counting switches or differentiating a minimizer. The first harmonic is
retained exactly. Corners are allowed; exponential decay is not claimed.

## The full-source payment is proved

Using the complete ratio-pair estimate already in the record, plus a bound
on the actual transverse derivative of the cosh representation, this stage
proves

    E_W[|X|*(g0+eta_N)] = O(N^(1/3)*d^3), d=ceil(log(N+1)).

It uses only the weak central norm O(N^(1/6)*d), not the sharper conditional
O(sqrt d) estimate. The polynomial representation error is paid separately
as a value error; its derivative is never inferred from that bound.

Choose Q_N=d^5*ceil(N^(1/3)). The discarded tail in the selection TAX is
then O(d^(-2)), uniformly in both original windows. This is an analytic
all-sufficiently-large-N reduction, not an extrapolation from a finite run.
No effective numerical onset or finite numerical constant is supplied.

The cutoff grows polynomially, not polylogarithmically. Computing all its
coefficients has not been made cheap, and no huge coefficient table was run.

## What still has to be signed

The exact angular mean survives unchanged, but the actual-height correction
must retain harmonics 1 through Q_N+3 after multiplication by W*X. Their
amplitudes and phases remain linked to the same physical height. The new
curvature moment matrices also impose joint positivity constraints; they
are not independent coefficient intervals.

The earlier moving-envelope counterexamples are recovered exactly at Q=4.
At Q=3 the nonzero missing response is paid, not averaged away. Thus this
reduction respects the obstruction found in the previous stage.

For a single fixed outer candidate, the prior direct-pair proof already
makes the correction o(1). The unresolved problem is its nonlinear adaptive
counterpart. The retained finite sum has not been proved favourable.

This route works directly with the exact full minimum. It does not apply
support curvature to the rational lower approximation without justification,
and it does not pay twice for two alternative approximations of that minimum.
The final unsigned-central and positive-gate obligations remain afterward.

## Verification

All 25 named exact tests pass. The suite covers the reserve shift, every member of a 65-candidate
transverse control, Fourier-tail telescoping, the previous moving examples,
256 low-mode preservation controls, deterministic spectral cutoff bounds,
128 exact corner-tail identities, and invalid-input guards.

Separate 40- and 70-digit interval calculations check 63 coefficients of
max(0,cos(theta)-1/2), a genuine switching affine envelope, plus six exact
corner-tail controls. These are analytic-example implementation checks,
not new RH source evaluations. The infinite-tail and full-source moment
proofs are written analytic arguments, not machine formalizations.

All 266 shared interval quantities overlap across precisions; 126
precision-specific conservative margins pass separately. The first
comparison verifier mistakenly required those different rounded margins
to overlap. Its failed source, log and receipt are retained, and two
regression controls cover the correction. No numerical source record changed.

See TEST_RESULT.json and VERIFICATION.json for exact counts, and the
external clean-extraction receipt for reproducibility. All historical
packages remain byte-identical. No source integration, height scan,
provider call, Forge experiment, repository push or external publication
was performed.

## Next target

Use the actual arithmetic and the shared positive curvature moments to
bound the RETAINED signed physical harmonics together with the angular
center, below the established outer drift reserve. The O(d^(-2)) discarded
tail is now controlled; the sign of the retained contribution is not.

This is progress on the reduction, not a proof of the favourable full
selection inequality or RH.
