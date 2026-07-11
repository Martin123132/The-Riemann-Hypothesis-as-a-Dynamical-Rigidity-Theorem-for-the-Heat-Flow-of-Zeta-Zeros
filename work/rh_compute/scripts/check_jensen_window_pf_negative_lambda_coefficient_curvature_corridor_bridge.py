#!/usr/bin/env python3
"""Validate the negative-lambda coefficient-curvature corridor bridge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
    exact_curvature_witness,
)


REQUIRED_ROW_IDS = {
    "nlcccb_01_curvature_coordinate",
    "nlcccb_02_exact_curvature_corridor",
    "nlcccb_03_monotone_side_identification",
    "nlcccb_04_repaired_k300_curvature_corridor",
    "nlcccb_05_monotone_curvature_shortcut_blocked",
    "nlcccb_06_raw_wall_shortcut_blocked",
    "nlcccb_07_live_lower_barrier_route",
    "nlcccb_08_conditional_closure_route",
    "nlcccb_09_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_stress",
    "exact_counterexample",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Coefficient-Curvature Corridor Bridge",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",
    "validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: 9 rows, 0 issues, 894 curvature-corridor rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows",
    "B_k = -log(x_k)",
    "p_k = log(R_k) = log((2*k+1)/(2*k-1))-B_k",
    "log((2*k+3)/(2+(2*k+1)*exp(-B_k))) <= B_(k+1)",
    "B_(k+1) <= B_k",
    "B positive rows: 897 / 897",
    "B upper-wall rows: 897 / 897",
    "monotone-curvature rows: 894 / 894",
    "lower-barrier rows: 894 / 894",
    "curvature-corridor rows: 894 / 894",
    "curvature-width rows: 894 / 894",
    "R_1=2, R_2=3/2: raw walls and B_(k+1)<=B_k hold, but the lower curvature barrier fails",
    "R_1=2, R_2=1: raw walls hold, but B_(k+1)<=B_k fails",
    "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
    "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
    "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
)


@dataclass(frozen=True)
class CurvatureBridgeIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> CurvatureBridgeIssue:
    return CurvatureBridgeIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[CurvatureBridgeIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[CurvatureBridgeIssue]:
    k, x, y = sp.symbols("k x y")
    raw = (2 * k + 1) * x / (2 * k - 1)
    raw_next = (2 * k + 3) * y / (2 * k + 1)
    raw_lower = (2 * k - 1) * (2 * k + 3) * raw / (2 * k + 1) ** 2
    raw_upper = (2 + (2 * k - 1) * raw) / (2 * k + 1)
    checks = {
        "lower-raw-to-monotone-x": sp.simplify((raw_next - raw_lower) - ((2 * k + 3) / (2 * k + 1)) * (y - x)),
        "upper-raw-to-scaled-x": sp.simplify((raw_upper - raw_next) - (2 + (2 * k + 1) * x - (2 * k + 3) * y) / (2 * k + 1)),
    }
    issues: list[CurvatureBridgeIssue] = []
    for name, value in checks.items():
        if sp.simplify(value) != 0:
            issues.append(issue("symbolic", f"bad-{name}", str(value)))
    return issues


def validate_top_level(artifact: dict) -> list[CurvatureBridgeIssue]:
    issues: list[CurvatureBridgeIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_raw_log_bridge",
        "source_k300_precision_repair",
        "source_raw_corridor_target",
        "source_log_curvature_bridge",
        "source_monotone_contraction_target",
        "source_raw_obstruction",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "raw-corridor", "cone entry", "jwpf_06", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[CurvatureBridgeIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_artifact([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[CurvatureBridgeIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_witnesses(rows_by_id: dict[str, dict]) -> list[CurvatureBridgeIssue]:
    issues: list[CurvatureBridgeIssue] = []
    expected = {
        "nlcccb_05_monotone_curvature_shortcut_blocked": exact_curvature_witness(Fraction(2), Fraction(3, 2)),
        "nlcccb_06_raw_wall_shortcut_blocked": exact_curvature_witness(Fraction(2), Fraction(1)),
    }
    for row_id, witness in expected.items():
        row = rows_by_id.get(row_id, {})
        if row.get("witness") != witness:
            issues.append(issue(row_id, "bad-witness", repr(row.get("witness"))))
        if not witness["raw_wall_holds_at_k"] or not witness["raw_wall_holds_at_next"]:
            issues.append(issue(row_id, "bad-raw-wall-witness", repr(witness)))
    first = expected["nlcccb_05_monotone_curvature_shortcut_blocked"]
    second = expected["nlcccb_06_raw_wall_shortcut_blocked"]
    if not first["B_monotone_holds"] or first["B_lower_barrier_holds"]:
        issues.append(issue("witness", "first-does-not-isolate-lower-barrier", repr(first)))
    if second["B_monotone_holds"]:
        issues.append(issue("witness", "second-does-not-fail-monotone-side", repr(second)))
    return issues


def validate_rows(artifact: dict) -> tuple[list[CurvatureBridgeIssue], int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[CurvatureBridgeIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))

    exact_available = 0
    exact_counterexamples = 0
    live_routes = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        readiness = row.get("readiness")
        if row.get("role") == "exact_reduction":
            exact_available += 1
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "exact_counterexample":
            exact_counterexamples += 1
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "open", "live", "conditional", "hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    issues.extend(validate_witnesses(rows_by_id))
    return issues, len(rows), exact_available, exact_counterexamples, live_routes if ready_to_apply == 0 else -live_routes


def validate_summary(artifact: dict, row_count: int, exact_available: int, exact_counterexamples: int, live_routes: int) -> list[CurvatureBridgeIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 9,
        "exact_available_rows": 3,
        "finite_stress_rows": 1,
        "exact_counterexample_rows": 2,
        "live_routes": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "coefficient_cap": 300,
        "b_total_rows": 897,
        "adjacent_total_rows": 894,
        "b_positive_rows": 897,
        "b_upper_wall_rows": 897,
        "b_monotone_rows": 894,
        "b_lower_barrier_rows": 894,
        "b_corridor_rows": 894,
        "b_width_rows": 894,
    }
    issues: list[CurvatureBridgeIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_available != 3:
        issues.append(issue("summary", "bad-exact-available-count", str(exact_available)))
    if exact_counterexamples != 2:
        issues.append(issue("summary", "bad-counterexample-count", str(exact_counterexamples)))
    if live_routes != 2:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "coefficient-curvature corridor",
        "b_(k+1) <= b_k",
        "897/897",
        "894/894",
        "monotone curvature alone",
        "exact cone counterexample",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "open zeta-specific", "finite stress", "monotone curvature", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[CurvatureBridgeIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[CurvatureBridgeIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "raw-corridor theorem is proved",
        "curvature corridor is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[CurvatureBridgeIssue], dict]:
    artifact = load_json(target_path)
    issues: list[CurvatureBridgeIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, exact_available, exact_counterexamples, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_available, exact_counterexamples, live_routes))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    issues, summary = validate(target, note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-CURV-CORRIDOR {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('b_corridor_rows')} curvature-corridor rows, "
            f"{summary.get('b_monotone_rows')} monotone-curvature rows, "
            f"{summary.get('exact_counterexample_rows')} exact counterexamples, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
