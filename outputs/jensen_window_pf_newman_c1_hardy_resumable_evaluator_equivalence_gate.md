# Newman C1 Hardy Resumable-Evaluator Equivalence Gate

Date: 2026-08-05

Status: finite computational reproducibility certificate validated at two
calibration heights; no physical-height value or proof error bound is claimed;
not a proof of RH.

## Implementation

The pinned upstream GPL-3.0 source remains unchanged.  The RH-owned derivative
retains the accepted `h=0.01` fifteen-value grid and adds a deterministic
three-stage restart state machine:

```text
stage 0: completed direct/Gaussian block
stage 1: completed Riemann--Siegel tail unit
stage 2: corrected final fifteen-value output
```

Each committed state contains the stage and unit, input and roster metadata,
all loop indices needed to resume, all fifteen direct partial sums, all fifteen
cumulative tail sums, and the latest fifteen tail contributions.  A mandatory
64-hex run identifier binds the Docker image, binary, checkpoint module,
evaluator source, input, MPI roster, and tail chunk policy.

After an atomic unit, rank zero writes a 41-significand-digit snapshot to a
temporary file, flushes and fsyncs it, atomically renames it over the
authoritative snapshot, fsyncs that file, appends the same state to JSONL, and
fsyncs the journal.  Runtime, deterministic-unit, and stop-file requests are
honored only after this commit.

The first executable run exposed a serialization defect: `ES48.40E4` fitted a
positive 41-digit value but not the minus sign of a negative value.  The final
format is the sign-safe `ES49.40E4`, and the source-contract checker forbids
reintroduction of the narrow format.

## Exact fixtures

The fixture-only tail chunk size is `257`, chosen to expose multiple restart
units without changing the mathematical terms.  The independent checker
validates:

| Height | States per path | Gaussian | Tail | Final | Parked invocations |
|---|---:|---:|---:|---:|---:|
| `10^10` | 40 | 36 | 3 | 1 | 4 |
| `10^12` | 56 | 47 | 8 | 1 | 2 |

The uninterrupted and multiply interrupted paths therefore contribute `192`
state records to the exact pairwise comparison.  Every arbitrary-precision
JSON number agrees, the final snapshots agree field by field, and all fifteen
displayed outputs agree exactly.  Both cases reproduce their previously
accepted ten-decimal output fixtures with maximum absolute difference `0`.

At `10^10`, two additional resume attempts alter respectively the run
identifier and the tail chunk size.  Both exit with status `2` and leave the
snapshot and journal unchanged.

## Proof Boundary

This gate proves finite stop/resume equivalence and fail-closed provenance for
the tested one-rank evaluator at `T=10^10` and `T=10^12`.  It does not certify
the external asymptotic approximation, supply a uniform error constant,
control correlated error in the fifteen shifted samples, evaluate the
physical carrier height, or provide interval arithmetic.  It proves no
retained observation, determinant sign or bound, complete-current inequality,
contact exclusion, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.

Validated command:

```text
python work/rh_compute/scripts/check_hardy_resumable_equivalence_fixture.py
```
