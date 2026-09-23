# Fixed-shift intake accepted; a finite signed current target

The incoming Signed Transport and Fixed Shift Fourier batches pass audit
within their stated analytic dependencies. Their four unchanged checker
outputs reproduce byte-for-byte. The concrete q=8192*pi remains fixed.

## New simplification

Rotate the same two frozen ordinary polynomials:

    Y_phi=sin(phi)V-cos(phi)Vq,
    Z_phi=cos(phi)V+sin(phi)Vq.

Then the complete response equals, up to the inherited vanishing shift
error, one half of the integral of DH(Y_phi)[Z_phi] for 0<=phi<=pi/2.
The first and third modes remain combined in H, with the original gate.
No unknown covariance, transport phase or inverse needs evaluation in
this formula. This does NOT set K to zero or assert independence; it is
an alternative exact identity, not a favourable signed estimate.

## A zero-safe error bound

The actual joint H has global gradient Lipschitz constant 11/10, verified
by two complete interval polynomial certificates. A rotating pair traces
a centered ellipse or a line, whose direction turns at most pi. This gives
a bounded-variation estimate even if the path passes through zero.

Consequently the whole signed integral can be replaced by M fixed
midpoint currents with asymptotic normalized error at most

    19 E(I0)/M^2.

No inverse-amplitude or fourth-moment estimate is needed. Rare large
values remain fully charged through the second energy.

Using the already-audited E(I0)<45/4 and M=32768, the allowance is

    855/4294967296 < 0.0000002.

The explicit next arithmetic target is the corresponding JOINT fixed-node
signed sum below 0.00092108989 in both original windows as N tends to
infinity. The nodes stay fixed before N. This signed bound is NOT proved
or numerically evaluated here; the stage removes coefficients and continuum
approximation obstacles, not the central arithmetic sign problem.

An auxiliary joint-covariance argument also improves the original transport
weight bound from 52/47 to 72/67. It likewise does not bound the current.

## Verification

- 30,579 incoming replay assertions; four outputs byte-identical.
- 14,533 new exact assertions, including two whole-angle Hessian certificates.
- 60 synthetic cases at 60 and 90 digits; 180 interval overlaps.
- 13 rejection guards.
- Every executed source version is retained and hash-authenticated.

See DERIVATION.md for all endpoint terms, inherited hypotheses, the
distributional midpoint error proof and the remaining signed obligation.
Historical evidence is unchanged. No source scan, provider call, Forge
experiment, GitHub update or external publication occurred.
