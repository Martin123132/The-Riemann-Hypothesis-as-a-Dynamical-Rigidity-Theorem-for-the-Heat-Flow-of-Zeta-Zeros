# Terminal sag is below the growing central credit

Private local RH research, completed 18 September 2026 local time.
Package date uses 17 September UTC turn start. No publication.

## Result

Let A_N=<X_-^2>_W be the original negative-central energy. On BOTH original
physical cutoff windows, uniformly in every deterministic terminal edge
b_N in [1/5,1/2], the terminal sag satisfies

    B_N^sag=<X_+ [min(C_N,Y_N)-min(C_N,B_N^term)]>_W
             <= C_1+C_2 sqrt(A_N)+o(1)=o(A_N).

This is an analytic asymptotic result using the authenticated source estimates,
not an extrapolation from successful finite windows. No effective onset or
numerical values for its uniform constants are supplied.

The incoming batch's stronger credit rate is accepted after reconciliation:

    liminf_N min_nu A_(N,nu)/(log log N) >=1/3.

The new terminal proof only needs A_N to diverge, which was already established
by the previous sealed stage. Thus its central conclusion does not depend on
improving the old rate to log log N.

## Mechanism

The endpoint-minus-minimum of a curve is at most its positive variation.
For the PAID short comparison Psi_N, that variation has the nonnegative bound

    Q_N=integral_(1/5)^(1/2) (partial_tau Psi_N+1)^2/4 d tau.

Only a quadratic arithmetic expression remains. Its signed central moment
<X_N Q_N> is bounded by the exact product AND divisibility resonances. Its
second moment <Q_N^2> is bounded by a multiplicative fourth-moment calculation.
Every remaining physical-height character is paid; no independent-phase law
is assigned to the source. Both moments are O(1), uniformly in N.

The exact source gap is <=Q_N+2 Delta_N, and <|X| Delta_N>=o(1). Therefore

    <X_+ Q_N>=<X Q_N>+<X_- Q_N>
               <=O(1)+sqrt(A_N <Q_N^2>).

That proves the claimed subcredit result. The factor two pays replacement of
both the endpoint and the minimum. No error estimate is differentiated.

The stronger claim <X_N(Y_N-B_N^term)>=O(1) for the EXACT minimum is still
open. It is not needed here. A counter-control proves why an upper envelope
alone would not justify that signed-order inference.

## Remaining Target

We can choose the terminal edge to be exactly alpha=1/5. Put

    L_N=min_(0<=tau<=1/5) F_N(4tau^2),
    B_N^alpha=min_(1/5<=tau<=1/2) F_N(4tau^2),
    K_N^alpha=(min(C_N,B_N^alpha)-L_N)_+.

The exact collected residual becomes

    B_N=<X_+ K_N^alpha>_W + o(A_N).

Thus the next selected-linear target is

    limsup_N max_nu <X_+ K_N^alpha>_W/(kappa_N A_N)<1.

A simpler sufficient target is

    limsup_N max_nu <X_+ K_N^alpha>_W/(log log N)<1/3.

The inner continuum is retained exactly; the source, common height, gate,
reflected terms and existing other payments are not changed. The next work
should attack this inner competition, not repeat the now-unneeded attempt
to transfer a general exact-minimum operator from a product torus.

This does NOT prove the full positive-gated score. Even its selected linear
component being favourable would not alone prove that score, zero counts or
RH. We do not claim B_sag=o(log log N), since only a LOWER bound for A_N at
that scale has been proved. The correct conclusion is B_sag=o(A_N).

## Checks and Integrity

* 21 named test families passed. They include exact Gaussian-rational
  character convolution, 72 mixed-moment comparisons, eight fourth-moment
  convolution checks, endpoint-error rejection controls, and 3,125 fixed-edge
  reanchoring controls. Imported pure scalar checkers were also replayed.
* Independent 50- and 80-digit mpmath.iv builds give 407 overlapping exported
  interval fields, including 36 fourth-moment factor-pair controls.
* The interval examples use actual finite coefficient formulas. They do not
  evaluate a physical-height source, approximate an integral over a window,
  or turn asymptotic Big-O constants into a finite bound.
* The two interval builds share mpmath.iv. The rational character algebra is
  a separate exact implementation. The analytic proof is written, not
  proof-assistant formalized.
* All 14 supplied hashes and 20 copied dependencies authenticated. Historical
  evidence is rechecked against 2,474 existing pins before final sealing.
* Routine wrapper failures and the successful earlier 20-test iteration are
  preserved in ITERATIONS.md, source snapshots, logs and receipts.

Clean extraction, replay hashes, seal identity and the final historical-pin
audit are in the delivery receipts and FINAL_CHECKPOINT.json. All compute
workers run singly, below normal priority, under 4 GiB and 900-second limits.
No provider call, Forge experiment, repository push or external publication.

## Files

DERIVATION.md supplies the proof. RECONCILIATION.md audits the incoming stages.
VERIFICATION.json and TEST_RESULT_FINAL.json give machine-readable checks.
The package manifest pins every delivered member. Historical records remain
byte-identical; the session bookmark receives only an additive checkpoint and
the updated next-action description.
