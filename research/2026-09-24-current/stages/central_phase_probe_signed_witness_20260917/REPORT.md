# Central phase probes and a joint adverse-overlap certificate

17 September 2026. Private local research. No publication or GitHub update.

## What moved

A mixed-response theorem now reaches the growing central variable X:

    E_W[X_N Psi] -> E_P[Re(e^(i chi) A_P) Psi],
    A_P=product_(p in P)(1-p^(-1/2)e^(i theta_p))^(-1).

P is a FIXED finite prime set and Psi is Lipschitz on that finite phase
torus. The theorem uses the inherited central L2 bound O(sqrt(log N)).
It pays the moving stationary indices and uses a growing Fejer approximation
to a fixed observable; it does not multiply an unchanging approximation
error by a divergent norm. The argument also covers bounded-complexity
families, such as a fixed-prefix minimum with a varying terminal band edge.

Bounded dual witnesses then give one-sided nonlinear-gate lower bounds for
nonnegative fixed observables. These are not a distributional limit, a
Gaussian model, equality for |X|, or independence of the full central sum.
The prime support cannot be enlarged to d^8 using this theorem for free.
See MIXED_PROBE_THEOREM.md for the proof and its hypotheses.

## A simpler sufficient signed condition

Keep the COMPLETE retained minimum m, including the middle region:

    H=min(kappa X,Y)-m, V=H-Y,
    T=E[X_+ H]-E[|X|Y]/2.

There is an exact identity

    T=E[XH]-E[XY]/2+E[X_- V].

Only the overlap X<0 AND V>0 can increase the last term. Put

    m_+=E[V_+], A_+=E[X V_+], B_+=E[X^2 V_+].

Then

    T <= E[XH]-E[XY]/2+(sqrt(m_+ B_+)-A_+)/2.

This formulation needs no explicit central absolute-value gate. Its inputs
are shared-height first/second moments, not independent worst cases for
each band. The positive part V_+ still contains the full source dependence
and has not been declared smooth or asymptotically small. The source,
collar, smoothing, transport and count obligations remain. See SIGNED_WITNESS.md.

## The complete saved-window test passes

All 65,536 children of the ONE saved N=2980 window were kept, with the same
W, X, Y and full minimum. Only the predeclared finite-prime features were
newly evaluated, over their whole cells, at 160 and 224 bits.

The positive adverse-region allowance is:

- Joint cap sqrt(E[V_+] E[X^2 V_+]): <0.16341.
- Global cap sqrt(E[X^2] E[V_+^2]): <0.68112.

All predeclared witness bounds are retained:

| Witness | Certified upper bound for T | Interpretation |
| --- | ---: | --- |
| Constant h=1 | <-0.64598 | Negative bound; best of these certificates |
| Clipped common phase | <-0.62019 | Negative bound |
| Clipped prime-2 predictor | <-0.60206 | Negative bound |
| Clipped primes 2,3,5,7 predictor | <-0.64135 | Negative bound |
| Constant h=-1 control | <+0.89283 | Inconclusive, not a positive actual T |

The prime predictors did NOT improve the constant certificate. No parameters
were tuned after these results. The new sufficient bound is also looser than
the preceding direct/rational enclosure. The gain is a simpler arithmetic
condition, not new verified heights or a stronger physical result.

## What still blocks the proof

The fixed-probe theorem cannot yet be applied to H or V_+, since they contain
the growing full minimum. Fixed-prime witnesses have finite responses, while
a separately estimated adverse cap can still be allowed to grow like
sqrt(log N). Finite negativity does not resolve that scale mismatch.

Keep G_+=E[X_- V_+] and G_-=E[X_- V_-]. The primary next target is

    E[XH]+G_+-G_- < mu_edge+c_*/2

on a defensible unbounded family, with every remaining payment. The simpler
one-sided moment cap is a fallback sufficient condition, not the main
asymptotic target: discarding G_- may discard the needed compensation.
No one-sided gate LOWER bound is used as an upper bound, and no theorem for
a fixed observable is assigned to a growing one.

The complete original-cell reaggregation retains both signs explicitly:

- Positive overlap G_+ lies in [0.08422, 0.10293].
- Negative overlap G_- lies in [0.04412, 0.05705].
- Net signed overlap lies in [0.02828, 0.05794].

These endpoints are outward-rounded. The signed result is still positive;
the negative part offsets some, not all, of the positive cost. The exact
reconstructed T remains negative and overlaps the prior direct enclosure.
This is a diagnostic on the same window, not an unbounded sign assertion.

## Verification and preservation

- 25/25 named tests passed, including the retained negative-credit check.
- 122,880 cross-precision parent-moment intervals agree.
- 196,608 cross-precision feature intervals agree.
- 122,880 independent integer parent-moment checks pass.
- 96 feature intervals checked separately with mpmath.iv at 70/100 digits.
- The constant-witness formula also passes independent expanded aggregation.

Arb runs share a library. The independent integer and mpmath computations
are separate elementary checks. The analytic arguments are written proofs,
not machine-formalized theorems. Original source integrations were not rerun.
The fresh-extraction replay receipt is saved next to the final sealed ZIP.

All 1,760 incoming historical hashes and all 17 snapshotted inputs are
checked again at the final checkpoint. Previous experiments and packages
are unchanged. No provider call, new Forge experiment or publication was
made. The bookmark records this completed stage with no background worker.
