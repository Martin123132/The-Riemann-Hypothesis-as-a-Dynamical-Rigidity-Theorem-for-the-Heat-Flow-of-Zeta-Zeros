#!/usr/bin/env python3
"""Validate the signed-Hankel/Jensen dependency graph.

The graph is a proof-safety artifact.  It checks that finite evidence and
countermodel gates feed into open theorem targets without recording any direct
proof edge to the Newman-direction conclusion Lambda <= 0.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_GRAPH = REPO_ROOT / "work/rh_compute/results/signed_hankel_jensen_dependency_graph.json"
DEFAULT_LEDGER = REPO_ROOT / "work/rh_compute/results/proof_claim_ledger.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/signed_hankel_jensen_dependency_graph.md"


ALLOWED_EDGE_TYPES = {
    "supports",
    "sharpens",
    "requires",
    "blocks_promotion",
    "validates_rejection",
    "would_imply_if_proved",
    "blocked_by",
    "documents",
}

FORBIDDEN_EDGE_TYPES = {"proves", "establishes", "solves", "implies"}

NON_PROVING_ROLES = {
    "exact_lemma",
    "finite_certificate",
    "diagnostic",
    "algebraic_reindexing",
    "countermodel_gate",
    "forbidden_promotion",
    "hygiene_gate",
}

REQUIRED_LEDGER_NODES = {
    "exact_heat_equation",
    "signed_hankel_finite_certificate",
    "hankel_sign_consistency_reduction_finite_certificate",
    "shifted_hankel_sign_consistency_finite_certificate",
    "jensen_hankel_bridge_algebra_gate",
    "jensen_window_pf_obligation_algebra_gate",
    "arb_jensen_window_pf_obligation_diagnostic",
    "arb_jensen_window_sturm_hyperbolicity_diagnostic",
    "arb_jensen_window_sturm_d5_hyperbolicity_diagnostic",
    "jensen_window_sturm_pf_finite_consequence",
    "countermodel_gates",
    "rejected_finite_prefix_promotion",
    "target_signed_hankel_jensen_bridge",
    "target_jensen_window_pf_bridge",
    "target_direct_coefficient_pf",
    "target_schur_positive_specialization",
    "target_positive_determinant_integral",
    "target_edrei_log_power_representation",
    "signed_hankel_jensen_dependency_graph",
    "core_reproducibility_gates",
}

OPEN_TARGET_NODES = {
    "target_signed_hankel_jensen_bridge",
    "target_jensen_window_pf_bridge",
    "target_direct_coefficient_pf",
    "target_schur_positive_specialization",
    "target_positive_determinant_integral",
    "target_edrei_log_power_representation",
}

FINITE_JENSEN_EVIDENCE = {
    "arb_jensen_window_pf_obligation_diagnostic",
    "arb_jensen_window_sturm_hyperbolicity_diagnostic",
    "arb_jensen_window_sturm_d5_hyperbolicity_diagnostic",
    "jensen_window_sturm_pf_finite_consequence",
}

REQUIRED_EDGES = {
    ("signed_hankel_finite_certificate", "target_signed_hankel_jensen_bridge", "supports"),
    ("shifted_hankel_sign_consistency_finite_certificate", "target_jensen_window_pf_bridge", "supports"),
    ("jensen_hankel_bridge_algebra_gate", "rejected_finite_prefix_promotion", "blocks_promotion"),
    ("jensen_window_pf_obligation_algebra_gate", "rejected_finite_prefix_promotion", "blocks_promotion"),
    ("countermodel_gates", "rejected_finite_prefix_promotion", "validates_rejection"),
    ("rejected_finite_prefix_promotion", "target_signed_hankel_jensen_bridge", "blocks_promotion"),
    ("rejected_finite_prefix_promotion", "target_jensen_window_pf_bridge", "blocks_promotion"),
    ("target_signed_hankel_jensen_bridge", "lambda_le_0_goal", "would_imply_if_proved"),
    ("target_jensen_window_pf_bridge", "lambda_le_0_goal", "would_imply_if_proved"),
    ("lambda_le_0_goal", "target_signed_hankel_jensen_bridge", "blocked_by"),
    ("lambda_le_0_goal", "target_jensen_window_pf_bridge", "blocked_by"),
}

REQUIRED_NOTE_STRINGS = (
    "# Signed-Hankel/Jensen Dependency Graph",
    "Status: dependency hygiene gate",
    "This is not a proof of PF-infinity",
    "work/rh_compute/results/signed_hankel_jensen_dependency_graph.json",
    "python work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py",
    "validated signed-Hankel/Jensen dependency graph with 0 issues",
    "lambda_le_0_goal",
    "status `not_proved`",
    "has no direct proving edge to `lambda_le_0_goal`",
    "would_imply_if_proved",
    "blocked_by",
    "work/rh_compute/results/proof_claim_ledger.json",
)


@dataclass(frozen=True)
class GraphIssue:
    section: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def ledger_claims(path: Path) -> dict[str, dict]:
    ledger = load_json(path)
    claims = ledger.get("claims", [])
    if not isinstance(claims, list):
        return {}
    return {str(claim.get("id")): claim for claim in claims if isinstance(claim, dict)}


def issue(section: str, name: str, detail: str) -> GraphIssue:
    return GraphIssue(section=section, issue=name, detail=detail)


def validate_top_level(graph: dict) -> list[GraphIssue]:
    issues: list[GraphIssue] = []
    if graph.get("kind") != "signed_hankel_jensen_dependency_graph":
        issues.append(issue("<graph>", "bad-kind", repr(graph.get("kind"))))
    boundary = str(graph.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<graph>", "weak-proof-boundary", graph.get("proof_boundary", "")))

    allowed = set(graph.get("allowed_edge_types", []))
    forbidden = set(graph.get("forbidden_edge_types", []))
    if not ALLOWED_EDGE_TYPES.issubset(allowed):
        issues.append(issue("<graph>", "missing-allowed-edge-types", repr(sorted(ALLOWED_EDGE_TYPES - allowed))))
    if not FORBIDDEN_EDGE_TYPES.issubset(forbidden):
        issues.append(issue("<graph>", "missing-forbidden-edge-types", repr(sorted(FORBIDDEN_EDGE_TYPES - forbidden))))
    overlap = allowed & forbidden
    if overlap:
        issues.append(issue("<graph>", "edge-type-overlap", repr(sorted(overlap))))
    return issues


def validate_nodes(graph: dict, ledger: dict[str, dict]) -> tuple[dict[str, dict], list[GraphIssue]]:
    issues: list[GraphIssue] = []
    nodes = graph.get("nodes", [])
    if not isinstance(nodes, list) or not nodes:
        return {}, [issue("nodes", "missing-nodes", "nodes must be a nonempty list")]

    by_id: dict[str, dict] = {}
    for node in nodes:
        if not isinstance(node, dict):
            issues.append(issue("nodes", "bad-node", repr(node)))
            continue
        node_id = str(node.get("id", ""))
        if not node_id:
            issues.append(issue("nodes", "missing-id", repr(node)))
            continue
        if node_id in by_id:
            issues.append(issue(node_id, "duplicate-node", node_id))
        by_id[node_id] = node

        role = node.get("role")
        status = node.get("status")
        if not role:
            issues.append(issue(node_id, "missing-role", repr(node)))
        if not status:
            issues.append(issue(node_id, "missing-status", repr(node)))

        if node.get("ledger_claim") is True:
            ledger_claim = ledger.get(node_id)
            if ledger_claim is None:
                issues.append(issue(node_id, "missing-ledger-claim", node_id))
                continue
            if role != ledger_claim.get("category"):
                issues.append(
                    issue(
                        node_id,
                        "ledger-role-mismatch",
                        f"graph {role!r} != ledger {ledger_claim.get('category')!r}",
                    )
                )
            if status != ledger_claim.get("status"):
                issues.append(
                    issue(
                        node_id,
                        "ledger-status-mismatch",
                        f"graph {status!r} != ledger {ledger_claim.get('status')!r}",
                    )
                )

    missing_required = REQUIRED_LEDGER_NODES - set(by_id)
    for node_id in sorted(missing_required):
        issues.append(issue(node_id, "missing-required-node", node_id))

    lambda_node = by_id.get("lambda_le_0_goal")
    if not lambda_node:
        issues.append(issue("lambda_le_0_goal", "missing-required-node", "lambda_le_0_goal"))
    else:
        if lambda_node.get("ledger_claim") is not False:
            issues.append(issue("lambda_le_0_goal", "bad-ledger-flag", repr(lambda_node.get("ledger_claim"))))
        if lambda_node.get("role") != "target_conclusion":
            issues.append(issue("lambda_le_0_goal", "bad-role", repr(lambda_node.get("role"))))
        if lambda_node.get("status") != "not_proved":
            issues.append(issue("lambda_le_0_goal", "bad-status", repr(lambda_node.get("status"))))
        boundary = f"{lambda_node.get('summary', '')} {lambda_node.get('proof_boundary', '')}".lower()
        if "not proved" not in boundary and "not_proved" not in boundary:
            issues.append(issue("lambda_le_0_goal", "weak-boundary", boundary))

    for node_id in sorted(OPEN_TARGET_NODES):
        node = by_id.get(node_id)
        ledger_claim = ledger.get(node_id)
        if node and node.get("status") != "open_target":
            issues.append(issue(node_id, "target-not-open-in-graph", repr(node.get("status"))))
        if ledger_claim and ledger_claim.get("status") != "open_target":
            issues.append(issue(node_id, "target-not-open-in-ledger", repr(ledger_claim.get("status"))))
    return by_id, issues


def validate_edges(graph: dict, nodes: dict[str, dict]) -> list[GraphIssue]:
    issues: list[GraphIssue] = []
    edges = graph.get("edges", [])
    if not isinstance(edges, list) or not edges:
        return [issue("edges", "missing-edges", "edges must be a nonempty list")]

    seen_edges: set[tuple[str, str, str]] = set()
    incoming: dict[str, list[dict]] = {}
    outgoing: dict[str, list[dict]] = {}

    for edge in edges:
        if not isinstance(edge, dict):
            issues.append(issue("edges", "bad-edge", repr(edge)))
            continue
        src = str(edge.get("from", ""))
        dst = str(edge.get("to", ""))
        edge_type = str(edge.get("type", ""))
        key = (src, dst, edge_type)
        seen_edges.add(key)
        incoming.setdefault(dst, []).append(edge)
        outgoing.setdefault(src, []).append(edge)

        if src not in nodes:
            issues.append(issue("edges", "missing-source", repr(edge)))
        if dst not in nodes:
            issues.append(issue("edges", "missing-target", repr(edge)))
        if edge_type not in ALLOWED_EDGE_TYPES:
            issues.append(issue("edges", "bad-edge-type", repr(edge)))
        if edge_type in FORBIDDEN_EDGE_TYPES:
            issues.append(issue("edges", "forbidden-edge-type", repr(edge)))
        if edge_type == "would_imply_if_proved":
            src_node = nodes.get(src, {})
            if src_node.get("role") != "theorem_target" or src_node.get("status") != "open_target":
                issues.append(issue(src, "conditional-edge-from-nontarget", repr(edge)))

        boundary = str(edge.get("boundary", "")).lower()
        if not boundary:
            issues.append(issue("edges", "missing-boundary", repr(edge)))
        if dst == "lambda_le_0_goal":
            src_role = nodes.get(src, {}).get("role")
            if edge_type not in {"would_imply_if_proved", "documents"}:
                issues.append(issue(src, "bad-edge-to-goal", repr(edge)))
            if src_role in NON_PROVING_ROLES and edge_type != "documents":
                issues.append(issue(src, "nonproving-node-points-to-goal", repr(edge)))
            if edge_type == "documents" and src != "signed_hankel_jensen_dependency_graph":
                issues.append(issue(src, "unexpected-documenting-edge-to-goal", repr(edge)))

    for required in sorted(REQUIRED_EDGES):
        if required not in seen_edges:
            issues.append(issue("edges", "missing-required-edge", repr(required)))

    for node_id in sorted(FINITE_JENSEN_EVIDENCE):
        targets = {
            edge.get("to")
            for edge in outgoing.get(node_id, [])
            if edge.get("type") == "supports"
        }
        if targets != {"target_jensen_window_pf_bridge"}:
            issues.append(issue(node_id, "bad-jensen-evidence-targets", repr(sorted(targets))))

    for target in sorted(OPEN_TARGET_NODES):
        conditional_edges = [
            edge for edge in outgoing.get(target, []) if edge.get("type") == "would_imply_if_proved"
        ]
        if not conditional_edges:
            issues.append(issue(target, "missing-conditional-route-edge", target))
        blocked_by_edges = [
            edge
            for edge in incoming.get(target, [])
            if edge.get("from") == "lambda_le_0_goal" and edge.get("type") == "blocked_by"
        ]
        if not blocked_by_edges:
            issues.append(issue(target, "missing-goal-blocked-by-edge", target))

    return issues


def validate_note(path: Path) -> list[GraphIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[GraphIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(graph_path: Path, ledger_path: Path, note_path: Path) -> list[GraphIssue]:
    issues: list[GraphIssue] = []
    graph = load_json(graph_path)
    ledger = ledger_claims(ledger_path)
    issues.extend(validate_top_level(graph))
    nodes, node_issues = validate_nodes(graph, ledger)
    issues.extend(node_issues)
    issues.extend(validate_edges(graph, nodes))
    issues.extend(validate_note(note_path))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(args.graph, args.ledger, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"DEPENDENCY-GRAPH {item.section} [{item.issue}] {item.detail}")
        print(f"validated signed-Hankel/Jensen dependency graph with {len(issues)} issues")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
