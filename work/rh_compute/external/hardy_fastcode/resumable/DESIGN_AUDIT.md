# Resumable Evaluator Design Audit

Status: implemented source contract; executable equivalence fixtures pending a
resource-policy-compliant Docker window.

## State Machine

| Stage | Unit | Meaning | Next unit |
| --- | ---: | --- | --- |
| `0` | `0` | Direct block zero is complete. | Gaussian block 1. |
| `0` | `k` | Gaussian blocks `0..k` are complete. | Gaussian block `k+1`. |
| `1` | `r` | All complete tail rounds `1..r` are accumulated. | Tail round `r+1`. |
| `1` | `R+1` | Complete tail rounds and the final fragment are accumulated. | Endpoint correction. |
| `2` | `0` | Endpoint correction is added and all 15 Hardy values are final. | None. |

For a fixed MPI size, tail round `r` assigns at most one disjoint term range
to each rank. One all-reduction closes the round before rank zero commits it.
The default range length is `100000`; validation uses `257` to expose several
tail boundaries at both calibrated heights.

## Preserved State

Every snapshot contains:

- source-format magic and a required 64-hex provenance hash;
- stage, unit, MPI size, `numbercalc`, `nchalf`, and `totblock`;
- `t`, `et`, `RN1`, `MT`, `MTO`, `aenums`, and `mmax`;
- `NC`, `N_MAX`, and the configured tail chunk size;
- all 15 `zsum`, `rszsumtot`, and latest `rszsum` values;
- three weighted full-precision array checksums.

`t`, `et`, grid metadata, MPI size, source format, provenance hash, and tail
chunk size are checked before a loaded value can affect the resumed loop. The
MPI size is part of provenance because it fixes the reduction tree and thus
the numerical path.

## Commit Protocol

1. Finish one Gaussian block or one tail round on every rank.
2. Complete the MPI reduction, leaving a common accumulated state.
3. Rank zero writes a 40-decimal-place real snapshot to `<path>.tmp`.
4. Flush, close, and `sync -f` the temporary file.
5. Atomically rename it over the authoritative snapshot and `sync -f` it.
6. Append the same state as one JSON object and `sync -f` the journal.
7. Broadcast commit status.
8. Only then inspect the stop file, elapsed-time limit, or deterministic test
   stop and broadcast the parking decision.
9. On a park, close any active PARI phase, close output, finalize MPI, and exit
   with status 75.

No queued unit is launched after a stop is observed. A malformed, mismatched,
or non-finite snapshot exits with status 2 before the stored recurrence state
is installed.

## Equivalence Matrix

The executable harness is configured to establish:

| Height | ZP blocks | RS cutoff | Forced RS units | Parking sequence |
| ---: | ---: | ---: | ---: | --- |
| `10^10` | 35 | 621 | 3 | runtime, block 5, stop file, RS unit 1 |
| `10^12` | 46 | 2034 | 8 | block 23, RS unit 1 |

For each height it compares the complete uninterrupted and resumed JSONL
records as exact arbitrary-precision decimals, compares the final snapshots,
and compares all 15 displayed values. It also checks those displayed values
against the accepted pre-restart fixture. Separate negative resumes alter the
provenance hash and chunk size and must exit without changing either canonical
checkpoint artifact.

## Deliberate Arithmetic Clarification

The upstream `pari_calc` declares its output array `INTENT(OUT)` and then adds
into it without first defining it. The derivative explicitly initializes that
array and accumulates each disjoint tail chunk into `rszsumtot` after a round.
This realizes the source's intended sum while making each tail unit
self-contained. The calibrated accepted runs have `N_MAX=0`, so their original
fixtures do not exercise the ambiguous multi-chunk path; the forced `257`
fixtures do.

This clarification needs both equivalence testing and comparison with an
independently constructed direct Riemann-Siegel main sum before a physical
diagnostic can be authorized.

## Proof Boundary

Restartability can establish deterministic preservation of this executable's
state. It cannot certify the upstream asymptotic error, provide a physical
height error constant, control the shifted-Hardy coefficient amplification,
or convert a diagnostic carrier value into an interval theorem. The physical
run remains barred until executable equivalence and the subsequent joint
correlated-error ledger both close.
