# Sharp two-mode control; the actual joint sign is still open

19 September2026. Private local RH research. Historical results unchanged.

## Main advance

We can now retain just the FIRST and THIRD actual angular moments and pay
for every higher angular mode, including its possible arithmetic coupling,
by an explicit scalar inequality. No prime-phase independence is assumed.

Writing J=Q/2+S for the same original weighted hinge, the strongest new
bound is

    J <= (1/2+v)Q + uE + alpha M_1 + beta M_3,

with alpha approximately0.42351467175516, beta approximately0.08128839630441,
u approximately0.01967638717558, and v approximately-0.01487331911601.
Here E and Q are the actual complex energy and real square, respectively.
The negative v is retained; no finite limiting-energy substitution occurs.

The envelope is SHARP at E=1,Q=1/2,M_1=M_3=0: an explicit three-point
probability measure attains its residual cost0.0122397276175733.
The proof is a pair of exact polynomial factorizations plus a matching
extremizer, not an optimization inferred from samples. A separate rational
envelope has a complete16-cell Bernstein certificate as a backup.

The initial uniform minimax cubic is also established. Its error is
(7sqrt(3)-12)/9. The quadratic-aware bound improves its normalized reserve
from0.0008866 to0.0011251.

## Remaining theorem

The inherited central credit and energy estimate leave the sufficient target

    limsup max_original_window
       [alpha M_1,N + beta M_3,N]/log N <0.0011251.

Fixed x-compacts come first, then the N limit, then the authenticated
REAL-source endpoint restoration. The finite complex source is not
integrated over the whole x-half-line before taking N to infinity.

THIS ACTUAL JOINT SIGN BOUND IS NOT PROVED. Both moments being asymptotically
nonpositive would suffice, but neither sign is assumed. This replaces the
unresolved higher-mode correlation ledger with two explicit nonlinear
moments and a paid remainder; it does not turn a reduction into an RH proof.
The original weight leaves eleven physical indices, not one.

## Actual finite checks

All four frozen controls use complete retained Dirichlet polynomials and
whole original height windows at x=1. They do not integrate the profile
variable and are not evaluations of the complete auxiliary function.
Both112- and160-bit runs enclose every height cell.

| N | Window | First moment | Third moment | Peeled third-mode coupling |
| --- | --- | --- | --- | --- |
| 64 | First | Positive | Negative | Positive |
| 64 | Second | Positive | Positive | Negative |
| 1024 | First | Positive | Unresolved at this enclosure width | Negative |
| 1024 | Second | Positive | Negative | Positive |

The fixed-profile two-mode combinations are positive in all four controls.
No finite result is compared to an asymptotic integrated-profile threshold.
The sign changes prevent a universal finite nonpositive-coupling shortcut;
they neither prove nor disprove an eventual asymptotic sign estimate.
Direct-versus-peeled third-mode differences are recorded separately. The
Euler peeling's asymptotic equality is never used as a finite equality.

## Whole-compact transfer

The previously sealed N1024 second-window certificate covers x in[1/8,16].
For R=uE+vQ and T=alpha M_1+beta M_3, the new inequality gives

    |S-T|<=R, hence Q/2+T+R <= J+2R.

Its authenticated finite E,Q,J enclosures therefore certify the new upper
majorant below0.036255, strictly below the historical0.04073 reference.
This is a deductive transfer of existing whole-compact evidence, NOT new
profile integration or separate determination of M_1 and M_3. The old
score is not changed. HISTORICAL_TRANSFER.json retains the exact inputs.

## Verification and scope

-10 exact test families,20181 scalar checks. The grid checks supplement,
 not replace, the analytic polynomial proofs.
-82 output interval comparisons and10240 midpoint-cache comparisons agree
 across precision. All20480 error-payment upper bounds remain nonnegative.
-228 component checks at12 fixed source points agree with independently
 coded direct-complex-power formulas, still using the same Arb library.
-11 guards reject unsupported strengthened claims.
-81920 whole height cells across the eight source builds; no sign-based
 selection or skipped source case.

The raw interval payments are precision-dependent upper bounds, not values
of a single invariant quantity; their equality or overlap is not required.
The first verifier incorrectly requested such overlap and was corrected.
The failed verifier source, log and receipt are retained. Source outputs
were not altered. AUDIT.md gives the details.

No provider calls, subagents, Forge experiments, repository push or external
publication. The user requested night mode; at most two active local workers
were used, below normal priority with4GiB process limits. Clean replay and
final integrity receipts are recorded separately after sealing.
