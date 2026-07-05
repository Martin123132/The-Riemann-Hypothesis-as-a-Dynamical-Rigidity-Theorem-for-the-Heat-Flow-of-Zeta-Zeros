# Core Proof-Programme Gate Runner

Date: 2026-07-04

Status: reproducibility ledger. This is not a proof of PF-infinity,
Laguerre-Polya membership, RH, or `Lambda <= 0`.

## Purpose

The corpus now has many individual finite-evidence manifests and
countermodel guards. The core runner gives a single command that checks the
main promoted finite certificates and the non-promotion gates used in the
referee-facing notes.

Script:

```text
work/rh_compute/scripts/check_core_proof_programme_gates.py
```

Command:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py
```

Current result:

```text
OK core gate: countermodel proof-safety gates
OK core gate: result-language boundary scan
OK core gate: output artifact status manifest
OK core gate: output reference integrity
OK core gate: proof-claim ledger
OK core gate: signed-Hankel finite certificates
OK core gate: Toeplitz/PF finite certificates
OK core gate: Toeplitz/Jacobi-Trudi reindexing
OK core gate: Hankel sign-consistency reduction point audits
OK core gate: Hankel sign-consistency reduction finite certificates
OK core gate: shifted Hankel sign-consistency finite certificates
OK core gate: Jensen/Hankel bridge algebra gate
OK core gate: signed-Hankel/Jensen bridge target specification
OK core gate: sign-regularity theorem fit matrix
OK core gate: positive Schur-specialization target note
OK core gate: Edrei-log sign diagnostics
OK core gate: Edrei power-Hankel diagnostics
OK core gate: Edrei midpoint frontier non-promotion guard
OK core gate: Edrei power-Hankel boundary repair manifest
OK core gate: Edrei moment-recurrence scout manifest
validated 20/20 core proof-programme gates
```

## Included Gates

```text
countermodel proof-safety gates:
  validates 10 exact/logical traps

result-language boundary scan:
  scans output markdown for obvious unqualified global-proof overclaims

output artifact status manifest:
  requires every output markdown file to declare Date, Status, artifact kind,
  and proof-boundary language near the top

output reference integrity:
  checks markdown references to scripts, output notes, and concrete result
  files; current run has 0 missing required paths and 3 planned deliverables

proof-claim ledger:
  validates 19 classified claims, including 5 open theorem targets that remain
  explicitly unproved

signed-Hankel finite certificates:
  validates 2,500 finite signed-Hankel certificates

Toeplitz/PF finite certificates:
  validates 95 promoted positive summaries and 1 negative control

Toeplitz/Jacobi-Trudi reindexing:
  validates N=10, orders<=5 exact map over 124,129 minors
  records 39,094 nonzero maps and 85,035 structural zeros

Hankel sign-consistency reduction point audits:
  validates 20 exact-rationalized finite reshaped-Hankel point audits
  for the five-lambda grid, orders k=2..5, and N=18

Hankel sign-consistency reduction finite certificates:
  validates 689,795 Arb/enclosure-backed reshaped-Hankel finite certificates
  for the five-lambda grid, orders k=2..7, and N=20

shifted Hankel sign-consistency finite certificates:
  validates 1,322,685 Arb/enclosure-backed shifted reshaped-Hankel finite
  certificates for the five-lambda grid, shifts n=0..20, orders k=2..5,
  and N=18

Jensen/Hankel bridge algebra gate:
  validates the exact degree-2 signed-Hankel/Jensen identity and an exact
  degree-3 finite countermodel blocking low-order finite-sign promotion

signed-Hankel/Jensen bridge target specification:
  validates that the active bridge target is stated as an all-order,
  all-shift open theorem with explicit low-degree obstructions

sign-regularity theorem fit matrix:
  validates that ASW/Edrei, Schur-positive, Schoenberg, Jensen, Hankel, and
  signed-Hankel route notes preserve fit/misfit boundaries and kill gates

positive Schur-specialization target note:
  validates that the coefficient-PF Schur target states the exact
  Jacobi-Trudi equivalence, noncircular sufficient constructions, and
  rejection gates for finite or circular proof attempts

Edrei-log sign diagnostics:
  validates 320 finite Edrei-log sign diagnostics

Edrei power-Hankel diagnostics:
  validates 4,205 finite Edrei power-Hankel diagnostics

Edrei midpoint frontier non-promotion guard:
  validates 5 non-rigorous midpoint frontier scouts

Edrei power-Hankel boundary repair manifest:
  validates 2 retired inconclusive blocker rows and 3 repaired positives

Edrei moment-recurrence scout manifest:
  validates 1 positive Arb recurrence scout and 1 inconclusive frontier scout
```

## Boundary

Passing this runner means the advertised finite ledgers and proof-safety
guards are internally reproducible, and that the exact Toeplitz-to-skew-Schur
reindexing used by the coefficient-PF route still passes its algebraic check.
It also checks for obvious forbidden result-language drift in the markdown
notes and enforces explicit artifact classification metadata. It does not
convert finite certificates, moment diagnostics, midpoint scouts,
counterexamples, reindexing identities, claim classification, clean wording,
clean metadata, or intact file references into an all-order bridge theorem.

The remaining theorem burden is unchanged:

```text
construct a noncircular all-order bridge from the coefficient/sign-regularity
data to PF-infinity, Laguerre-Polya membership, Jensen hyperbolicity for all
degrees/shifts, or directly to Lambda <= 0.
```

## Useful Options

For quick checks that skip gates marked slow:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py --skip-slow
```

Current quick result:

```text
validated 19/19 core proof-programme gates
```

For machine-readable output:

```text
python work/rh_compute/scripts/check_core_proof_programme_gates.py --json
```
