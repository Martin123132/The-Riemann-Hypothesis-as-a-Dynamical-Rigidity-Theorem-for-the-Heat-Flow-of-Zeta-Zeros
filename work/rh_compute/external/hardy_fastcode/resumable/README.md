# Resumable Hardy Multi Evaluator

This directory contains an RH-owned derivative of the pinned GPL-3.0
`zeta14cubicmult.f90` diagnostic evaluator. It is not a proof oracle and must
not be used as an interval evaluator without a separate validated error bound.

## Provenance

- Upstream repository: `https://github.com/dml2391/Hardy-function-fastcodes`
- Pinned commit: `2e16dac3206b707052c3ac4cacdf3d1a2325e636`
- License: GPL-3.0; a local copy is in `COPYING.GPL-3.0.txt`.
- Base source SHA-256 after portability patches 0003 and 0004:
  `5d0699865ab58f6968adfbe11ed969dafaabf5e965d5bd74391206af8b58b620`.
- Retained shift spacing: exactly `0.01`. The rejected `0.02` and `0.04`
  experimental grids are not incorporated.

`zeta14cubicmult_resumable.f90` starts from that patched base. Its arithmetic
changes are limited to restart ownership and an explicit zero initialization
of each `pari_calc` output array. The latter makes the subroutine's declared
`INTENT(OUT)` contract explicit before accumulation.

## Atomic State

`rh_hardy_checkpoint.f90` writes two records after every completed unit:

1. A full-precision snapshot is written to a temporary file, flushed, fsynced,
   and atomically renamed over the previous snapshot.
2. The same state is appended as one JSON object to an fsynced JSONL journal.

Each record preserves all 15 `zsum`, `rszsumtot`, and `rszsum` channels plus
the Gaussian-block recurrence state. A resume fails closed unless the height,
tolerance, number of samples, total block count, communicator size, and source
format match. The Riemann-Siegel tail is partitioned into deterministic MPI
rounds, so it can also resume after a committed chunk.

The sign-safe `ES49.40E4` records carry 41 significand digits. The journal records the
actual radix and model digit count so the checker can verify that this exceeds
the round-trip requirement for the selected `real(kind=selected_real_kind(33))`.

The snapshot is authoritative. The JSONL file is an append-only audit trail;
an interrupted final JSONL line does not invalidate the last atomically
renamed snapshot.

## Controls

The executable reads these environment variables:

| Variable | Meaning |
| --- | --- |
| `RH_HARDY_CHECKPOINT` | Snapshot path; default `rh_hardy_checkpoint.dat`. |
| `RH_HARDY_JOURNAL` | JSONL path; default `<snapshot>.jsonl`. |
| `RH_HARDY_STOP_FILE` | Stop-file path; default `<snapshot>.stop`. |
| `RH_HARDY_RESUME` | `1/true/yes` loads the existing snapshot. |
| `RH_HARDY_RUN_ID` | Required 64-hex hash of image, binary, source, input, MPI size, and chunk policy. |
| `RH_HARDY_MAX_SECONDS` | Fresh wall-clock budget for this invocation; zero disables it. |
| `RH_HARDY_STOP_AFTER_STAGE` | Deterministic test/control stage: `0` for ZP, `1` for RS. |
| `RH_HARDY_STOP_AFTER_UNIT` | Park after this committed unit in the selected stage. |
| `RH_HARDY_RS_CHUNK_SIZE` | Tail terms per rank and round; default `100000`. |

A stop request is honored only after the current atomic block or tail round is
committed. The program then finalizes MPI and exits with status 75. A completed
run writes stage 2 and exits normally. Resume requires the same MPI process
count because reduction order is part of the numerical provenance.

## Build

Inside the pinned Debian Bookworm container image:

```sh
mpif90 -O3 -z noexecstack -fallow-argument-mismatch \
  -ffree-line-length-none \
  rh_hardy_checkpoint.f90 zeta14cubicmult_resumable.f90 \
  -lpari -o zeta14cubicmult.resumable
```

Run from a directory containing the fixed-format `inputs3.nml` and an existing
`zeta14v6resa` output file. Initial validation uses one MPI rank so comparison
with the accepted `h=0.01` fixtures has a fixed reduction order. The resume
fixtures set `RH_HARDY_RS_CHUNK_SIZE=257` to force multiple tail checkpoints at
the calibrated heights; production diagnostics retain the `100000` default.

## Proof Boundary

Restart equivalence establishes only that interruption does not change this
diagnostic program's numerical path for a fixed binary, input, and MPI size.
It does not validate the upstream asymptotic error estimate, supply its missing
constant at physical height, control the shifted-kernel coefficient
amplification, or prove any Riemann-hypothesis implication.
