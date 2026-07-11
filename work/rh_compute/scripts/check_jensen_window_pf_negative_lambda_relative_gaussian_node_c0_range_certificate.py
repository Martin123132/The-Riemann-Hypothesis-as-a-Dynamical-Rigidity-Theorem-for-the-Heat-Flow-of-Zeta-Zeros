#!/usr/bin/env python3
"""Validate the relative-Gaussian node-range and Phi(0) lower certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate import (  # noqa: E402
    DEFAULT_GRID_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TAIL_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgnc_01_laguerre_jacobi_reduction",
    "nlrgnc_02_gershgorin_node_range",
    "nlrgnc_03_phi0_n1_lower_certificate",
    "nlrgnc_04_phi_tail_conditions_handoff",
    "nlrgnc_05_weight_and_quadrature_promotion_rejected",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "exact_bound",
    "arb_lower_bound",
    "conditional_handoff",
    "rejected_route",
}

ALLOWED_READINESS = {
    "available_exact",
    "available_arb",
    "available_for_tail_scout",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Node-C0 Range Certificate",
    "Status: node-range and Phi0 lower certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: 5 rows, 0 issues, 16 Laguerre bound rows, 2 certified side conditions, 0 ready-to-apply rows",
    "worst Laguerre node upper bound: 809",
    "worst x^2 upper bound: 809/1156",
    "node range x<=1 certified: True",
    "certified c0 lower: 0.44",
    "c0 lower certified by n=1 term: True",
    "2*pi^2*n^4 - 3*pi*n^2 = pi*n^2*(2*pi*n^2-3)>0 for n>=1, using pi>3",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class NodeC0Issue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> NodeC0Issue:
    return NodeC0Issue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("E+00", "E+0"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[NodeC0Issue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def validate_top_level(artifact: dict) -> list[NodeC0Issue]:
    issues: list[NodeC0Issue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "node-range and Phi0 lower certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_cancellation_reduced_grid_scout",
        "source_cancellation_reduced_grid_json",
        "source_phi_tail_bound_scout",
        "source_phi_tail_bound_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "node-range and phi(0) lower certificate only",
        "does not provide individual laguerre node or weight intervals",
        "does not bound quadrature",
        "uniform collar",
        "scaled-curvature",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict, grid_path: Path, phi_tail_path: Path, interval_path: Path) -> list[NodeC0Issue]:
    try:
        recomputed = build_artifact(grid_path, phi_tail_path, interval_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[NodeC0Issue] = []
    for key in ("laguerre_bound_rows", "matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_laguerre_rows(artifact: dict) -> list[NodeC0Issue]:
    rows = artifact.get("laguerre_bound_rows", [])
    issues: list[NodeC0Issue] = []
    if not isinstance(rows, list):
        return [issue("laguerre_bound_rows", "bad-rows", repr(type(rows)))]
    if len(rows) != 16:
        issues.append(issue("laguerre_bound_rows", "bad-row-count", str(len(rows))))
    expected_orders = [64, 96, 128, 192]
    expected_indices = [21, 22, 23, 24]
    seen = {(row.get("quadrature_order"), row.get("index")) for row in rows if isinstance(row, dict)}
    if seen != {(order, index) for order in expected_orders for index in expected_indices}:
        issues.append(issue("laguerre_bound_rows", "bad-order-index-grid", repr(sorted(seen))))
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("laguerre_bound_rows", "bad-row", repr(row)))
            continue
        row_id = f"N={row.get('quadrature_order')},i={row.get('index')}"
        for key in (
            "alpha",
            "laguerre_largest_node_upper_bound",
            "min_T",
            "x_square_upper_bound",
            "x_less_than_one_certified",
            "proof_boundary",
        ):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("min_T") != 1156:
            issues.append(issue(row_id, "bad-min-T", repr(row.get("min_T"))))
        if row.get("x_less_than_one_certified") is not True:
            issues.append(issue(row_id, "node-range-not-certified", repr(row)))
        if dec(row.get("x_square_upper_bound_decimal")) >= Decimal(1):
            issues.append(issue(row_id, "x-square-bound-too-large", repr(row.get("x_square_upper_bound_decimal"))))
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("gershgorin", "laguerre nodes only", "does not enclose"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-row-boundary", required))
    worst = max(rows, key=lambda row: dec(row["x_square_upper_bound_decimal"]))
    expected_worst = {
        "quadrature_order": 192,
        "index": 24,
        "alpha": "47/2",
        "laguerre_largest_node_upper_bound": "809",
        "x_square_upper_bound": "809/1156",
    }
    for key, value in expected_worst.items():
        if worst.get(key) != value:
            issues.append(issue("laguerre_bound_rows", f"bad-worst-{key}", f"{worst.get(key)!r} != {value!r}"))
    return issues


def validate_summary(artifact: dict) -> list[NodeC0Issue]:
    summary = artifact.get("summary", {})
    issues: list[NodeC0Issue] = []
    expected = {
        "matrix_rows": 5,
        "laguerre_bound_rows": 16,
        "certified_side_conditions": 2,
        "ready_to_apply_rows": 0,
        "source_grid_rows": 20,
        "source_phi_tail_rows": 6,
        "source_interval_obligations": 8,
        "min_T": 1156,
        "max_quadrature_order": 192,
        "max_index": 24,
        "worst_laguerre_node_upper_bound": "809",
        "worst_x_square_upper_bound": "809/1156",
        "worst_x_square_upper_bound_decimal": "6.998269896193771711E-01",
        "worst_x_square_slack_to_one": "347/1156",
        "node_range_x_le_1_certified": True,
        "certified_c0_lower": "0.44",
        "c0_lower_certified_by_n1_term": True,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if not (Decimal("0.69") < dec(summary.get("worst_x_square_upper_bound_decimal")) < Decimal("0.70")):
        issues.append(issue("summary", "bad-worst-x-square-range", repr(summary.get("worst_x_square_upper_bound_decimal"))))
    if not dec(summary.get("n1_phi0_term_lower")) > Decimal("0.445"):
        issues.append(issue("summary", "n1-lower-too-small", repr(summary.get("n1_phi0_term_lower"))))
    if not dec(summary.get("n1_minus_c0_lower_lower")) > Decimal("0.005"):
        issues.append(issue("summary", "c0-margin-too-small", repr(summary.get("n1_minus_c0_lower_lower"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("gershgorin", "809<t_min=1156", "x<=1", "arb certifies", "0.44", "weights", "grid-to-collar"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("x<=1", "phi(0)>=0.44", "laguerre weights", "uniform collar", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_rows(artifact: dict) -> list[NodeC0Issue]:
    rows = artifact.get("matrix_rows", [])
    issues: list[NodeC0Issue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))]
    if len(rows) != 5:
        issues.append(issue("matrix_rows", "bad-row-count", str(len(rows))))
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") not in ALLOWED_READINESS:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    c0_diag = rows_by_id.get("nlrgnc_03_phi0_n1_lower_certificate", {}).get("diagnostics", {})
    try:
        n1_term = flint.arb(str(c0_diag.get("n1_phi0_term_ball")))
        margin = flint.arb(str(c0_diag.get("n1_minus_c0_lower_ball")))
    except Exception as exc:
        issues.append(issue("nlrgnc_03_phi0_n1_lower_certificate", "bad-arb-ball", f"{type(exc).__name__}: {exc}"))
    else:
        if not arb_positive(n1_term - flint.arb("0.44")):
            issues.append(issue("nlrgnc_03_phi0_n1_lower_certificate", "n1-not-above-c0-lower", str(n1_term)))
        if not arb_positive(margin):
            issues.append(issue("nlrgnc_03_phi0_n1_lower_certificate", "margin-not-positive", str(margin)))
    handoff = rows_by_id.get("nlrgnc_04_phi_tail_conditions_handoff", {}).get("diagnostics", {})
    if handoff.get("certified_side_conditions") != ["node-induced x<=1", "Phi(0)>=0.44"]:
        issues.append(issue("nlrgnc_04_phi_tail_conditions_handoff", "bad-certified-side-conditions", repr(handoff)))
    if len(handoff.get("remaining_intervalization_obligations", [])) != 5:
        issues.append(issue("nlrgnc_04_phi_tail_conditions_handoff", "bad-remaining-obligations", repr(handoff)))
    return issues


def validate_note(path: Path) -> list[NodeC0Issue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[NodeC0Issue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "finite-grid interval certificate is complete",
        "uniform residual estimate is proved",
        "quadrature-remainder theorem is proved",
        "laguerre weight intervals are certified",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(
    target_path: Path,
    note_path: Path,
    grid_path: Path,
    phi_tail_path: Path,
    interval_path: Path,
) -> tuple[list[NodeC0Issue], dict]:
    artifact = load_json(target_path)
    issues: list[NodeC0Issue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, grid_path, phi_tail_path, interval_path))
    issues.extend(validate_laguerre_rows(artifact))
    issues.extend(validate_summary(artifact))
    issues.extend(validate_rows(artifact))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    issues, summary = validate(target, note, grid_path, phi_tail_path, interval_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-NODE-C0 {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('laguerre_bound_rows')} Laguerre bound rows, "
            f"{summary.get('certified_side_conditions')} certified side conditions, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
