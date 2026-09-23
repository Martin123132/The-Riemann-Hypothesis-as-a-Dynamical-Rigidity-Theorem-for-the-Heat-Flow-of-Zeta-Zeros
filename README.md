# Riemann Zeta Zero Dynamics Under Heat Flow

> **PROPRIETARY - WRITTEN PERMISSION REQUIRED**
> Copyright (c) 2026 TWO HANDS NETWORK LTD. All rights reserved.
> Project creator: Martin Ollett. Company-owned protected material is not
> open source and has no general non-commercial-use permission. Prior written
> Company permission is required for licensed reuse; attribution alone is not
> permission. Existing valid licences, third-party rights, statutory exceptions
> and hosting-platform rights remain unaffected. See [LICENSE.md](LICENSE.md),
> [NOTICE.md](NOTICE.md), and [permission requests](COMMERCIAL-LICENSE.md).

This repository is a research corpus on the Riemann Hypothesis, the
de Bruijn-Newman heat flow, Jensen-window and determinant methods, and a
current Xi-specific Hardy endpoint programme.

> **Status:** this repository does not contain a proof of the Riemann
> Hypothesis or a proof that `Lambda <= 0`. It contains exact reductions,
> rigorous finite and interval certificates, countermodel gates, computational
> evidence, and explicitly labelled open theorem targets.

The corpus preserves the original manuscripts and numerical experiments, but
the current audited programme is documented in `outputs/` and
`work/rh_compute/`.

## Start Here

- [Current completed research, 24 September 2026](research/2026-09-24-current/README.md)
  contains the full latest derivation and audit, an index of 90 recorded stages,
  and lossless hash-verified evidence through Uniform Filtered First Moment.
  The favourable signed-current bound remains open; this is not a proof of RH.
- [`outputs/Clay_Prize_Readiness_Audit.md`](outputs/Clay_Prize_Readiness_Audit.md)
  gives the honest prize-readiness verdict and the missing theorem chain.
- [`outputs/formal_core.md`](outputs/formal_core.md) is the cumulative formal
  mathematical core.
- [`outputs/proof_claim_ledger.md`](outputs/proof_claim_ledger.md) classifies
  proved statements, finite certificates, diagnostics, countermodels, and open
  targets.
- [`outputs/core_proof_programme_gates.md`](outputs/core_proof_programme_gates.md)
  explains the executable proof-safety gate registry.
- [`outputs/RH_Proof_Programme_Roadmap.md`](outputs/RH_Proof_Programme_Roadmap.md)
  records the broader programme and route history.

## Research Update: 24 September 2026

The [new snapshot](research/2026-09-24-current/README.md) brings the completed
16-23 September research into the repository without changing previous sealed
evidence. It includes complete mathematical files, code and saved verification,
with full-byte hash checks and explicit administrative exclusions. Original
archives are inventoried and their included members are losslessly repackaged;
the original local ZIP containers and seals are unchanged.

The latest result supplies uniform filtered first-moment decay for the literal
rough source under the retained reciprocal-phase estimate, and a squared
vanishing nonlinear-overlap allowance. The remaining signed alignment is not
proved. The separate next-stage GPT work is unfinished and is not included.
See the [current audit](research/2026-09-24-current/latest/AUDIT.md) for precise
dependencies and verification limits. Earlier checkpoints below are historical.

## Research Update: 16 September 2026

The current selection-stage audit and immutable evidence are in [research/2026-09-16-selection/AUDIT.md](research/2026-09-16-selection/AUDIT.md). The four packages pass 142 tests and 17 byte-identical result rebuilds. They establish finite checks and scoped asymptotic deductions, not RH. The unbounded gated coefficient-minimum inequality remains open. The earlier checkpoint below is retained as historical context.

## Archived Checkpoint

The formal core currently runs through Section 11.521. Its newest continuous
certificate proves, for the exact Hardy/Newman endpoint object used there,

```text
Q_K(t) - T(t) < 0
for every t in [10^10 + 12.5, 10^10 + 13].
```

The proof uses an anchored grouped 752-label phase transport with retained
panel and derivative errors. Independent production and reverse-order routes
give strictly negative upper margins. No stationary-point uniqueness is
assumed.

This is one certified phase cell, not the global Newman theorem. The next
local target is the adjacent cell `[10^10+13, 10^10+13.5]`. The programme must
still connect the remaining event-cell cover, prove both wall handoffs, obtain
uniform control across all required heights and changing rosters, and deduce
that no positive Newman boundary can occur.

## Repository Map

```text
outputs/
  Formal notes, theorem gates, route audits, countermodels, and claim ledgers.

work/rh_compute/scripts/
  Generators, rigorous checkers, independent replays, and bounded runners.

work/rh_compute/results/
  Machine-readable certificates, interval outputs, ledgers, and caches.

work/rh_compute/external/
  Pinned external-source metadata, GPL-3.0 Hardy adapters, and exact fixtures
  used by the reproducibility gates.
```

Large certificate ledgers are stored with Git LFS. Virtual environments,
runtime policy, PID files, session bookmarks, scratch files, transient logs,
and two multi-gigabyte resumable ledgers are intentionally not published.

## Core Integrity Checks

From the repository root, the main lightweight checks are:

```powershell
python work/rh_compute/scripts/check_proof_claim_ledger.py
python work/rh_compute/scripts/check_result_language_boundaries.py
python work/rh_compute/scripts/check_output_status_manifest.py
python work/rh_compute/scripts/check_output_reference_integrity.py
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
```

The serial core replay is intentionally resource bounded:

```powershell
python work/rh_compute/scripts/run_core_gates_resource_bounded.py --mode day
```

Passing finite, static, or computational gates does not promote an open target
to a theorem. The proof-claim ledger and dependency graph enforce that
boundary.

## Route History

The original zero-motion, rank-reduction, and finite-flow experiments remain
part of the historical corpus. Later work established that local repulsion by
itself cannot exclude a positive square-root birth at a Newman boundary. An
all-shift signed-Hankel endpoint route was also rejected by rigorous order-ten
counterexamples. Those failures are retained because they prevent circular or
overstrong arguments from re-entering the programme.

The live route is the Xi-specific Hardy/Newman endpoint comparison developed
in the later formal core. Its status should be read from the readiness audit,
not inferred from the size of the computational corpus.

## Licensing

Company-owned protected material is proprietary under
[LICENSE.md](LICENSE.md); no general commercial or non-commercial licence
is offered. Prior written Company permission is required where a licence is
needed. Attribution alone is not permission. See [NOTICE.md](NOTICE.md) and
[written permission requests](COMMERCIAL-LICENSE.md).

[LICENSING-HISTORY.md](LICENSING-HISTORY.md) preserves the earlier licensing
record. Existing valid grants, legal exceptions, hosting-platform rights,
and third-party licences remain unaffected. Separately licensed material,
including `work/rh_compute/external/`, retains its provenance and licence.
This policy does not claim exclusive rights in mathematical ideas or facts.
