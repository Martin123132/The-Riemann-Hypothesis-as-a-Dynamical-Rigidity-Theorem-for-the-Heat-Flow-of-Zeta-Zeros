# Signed-Hankel/Jensen Dependency Graph

Date: 2026-07-05

Status: dependency hygiene gate. This is not a proof of PF-infinity,
Laguerre-Polya membership, RH, or `Lambda <= 0`; it records how finite
evidence, countermodel gates, and open theorem targets are allowed to depend on
one another.

## Purpose

The signed-Hankel and Jensen-window route now has several exact algebra gates,
finite Arb diagnostics, and countermodel guards. The dependency graph makes
the proof boundary executable: finite evidence may support open theorem
targets, and countermodels may block invalid promotions, but finite evidence
has no direct proving edge to `lambda_le_0_goal`.

Machine-readable graph:

```text
work/rh_compute/results/signed_hankel_jensen_dependency_graph.json
```

Checker:

```text
python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py
```

Current result:

```text
validated signed-Hankel/Jensen dependency graph with 0 issues
```

## Core Invariant

The graph has an explicit conclusion node:

```text
lambda_le_0_goal:
  status `not_proved`
```

The only edges into that node are `would_imply_if_proved` edges from ledger
nodes whose status is still `open_target`, plus documentation edges from this
hygiene gate. The graph also has `blocked_by` edges back from
`lambda_le_0_goal` to the missing theorem targets, so the direction of the
remaining burden stays visible.

## Active Bridge Region

```text
signed_hankel_finite_certificate
hankel_sign_consistency_reduction_finite_certificate
shifted_hankel_sign_consistency_finite_certificate
arb_jensen_window_pf_obligation_diagnostic
arb_jensen_window_sturm_hyperbolicity_diagnostic
arb_jensen_window_sturm_d5_hyperbolicity_diagnostic
jensen_window_sturm_pf_finite_consequence
```

These nodes support:

```text
target_signed_hankel_jensen_bridge
target_jensen_window_pf_bridge
```

Both targets remain open. The graph checker compares their statuses with
`work/rh_compute/results/proof_claim_ledger.json`.

## Non-Promotion Region

```text
jensen_hankel_bridge_algebra_gate
jensen_window_pf_obligation_algebra_gate
countermodel_gates
rejected_finite_prefix_promotion
```

These nodes only sharpen the theorem targets or validate rejection of finite
prefix promotion. They do not supply a theorem proving path to the conclusion.

## Validation Rules

The checker enforces:

```text
all ledger-backed graph nodes exist in the proof-claim ledger
ledger-backed graph roles and statuses match the ledger categories and statuses
all edge endpoints exist
all edge types are in the allowed set
no forbidden edge type appears
lambda_le_0_goal remains not_proved
would_imply_if_proved edges are sourced only from open theorem targets
finite and diagnostic nodes do not point directly to lambda_le_0_goal
the Jensen-window evidence nodes feed only into open bridge targets
the countermodel nodes block finite-prefix promotion
```

## Boundary

Passing this checker means the dependency ledger has not silently promoted a
finite signed-Hankel, finite Jensen-window, Sturm, PF, countermodel, or hygiene
artifact into the missing Newman-direction theorem. It is a proof-programme
firewall, not a proof of the firewall's target.
