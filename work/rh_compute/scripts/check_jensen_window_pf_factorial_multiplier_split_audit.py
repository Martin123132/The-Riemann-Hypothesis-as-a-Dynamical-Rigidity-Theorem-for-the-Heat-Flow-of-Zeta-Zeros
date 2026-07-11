#!/usr/bin/env python3
"""Validate the Jensen-window PF factorial multiplier split audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_AUDIT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_factorial_multiplier_split_audit.md"


EXPECTED_ROW_IDS = (
    "fm_01_factorial_multiplier_sequence",
    "fm_02_shifted_window_preserver",
    "fm_03_raw_moment_degree2_anti_alignment",
    "fm_04_normalized_degree2_threshold",
    "fm_05_route_verdict",
)

REQUIRED_PATHS = (
    "outputs/coefficient_pf_bridge_obstruction.md",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md",
    "work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json",
    "work/rh_compute/scripts/jensen_window_pf_factorial_multiplier_split_audit.py",
    "work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py",
)

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Factorial Multiplier Split Audit",
    "Status: finite theorem-search diagnostic",
    "Artifact kind: `jensen_window_pf_factorial_multiplier_split_audit`.",
    "work/rh_compute/results/jensen_window_pf_factorial_multiplier_split_audit.json",
    "python work/rh_compute/scripts/jensen_window_pf_factorial_multiplier_split_audit.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py",
    "validated Jensen-window PF factorial multiplier split audit: 5 exact rows, 315 raw degree-2 anti-hyperbolic rows, 315 normalized degree-2 positive rows, 0 ready-to-apply rows, 0 issues",
    "gamma_k = k!/(2*k)!",
    "sum gamma_k*z^k/k! = cosh(sqrt(z))",
    "Delta_2(M_{2,n}) = 4*(mu_{2*n+2}^2 - mu_{2*n}*mu_{2*n+4}) <= 0",
    "A_{n+1}^2 >= A_n*A_{n+2}",
    "iff mu_{2*n+2}^2/(mu_{2*n}*mu_{2*n+4}) >= (2*n+1)/(2*n+3)",
    "raw degree-2 anti-hyperbolic rows: 315",
    "normalized degree-2 positive rows: 315",
    "Raw moment positivity is not enough.",
)

FORBIDDEN_TEXT = (
    "therefore rh",
    "we have proved lambda <= 0",
    "proves lambda <= 0",
    "jensen-window pf-infinity is proved",
    "the sign-regular-to-jensen transfer theorem is proved",
    "raw moment positivity proves",
)


@dataclass(frozen=True)
class AuditIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require(condition: bool, issues: list[AuditIssue], row_id: str, issue: str, detail: str) -> None:
    if not condition:
        issues.append(AuditIssue(row_id, issue, detail))


def validate_paths(audit: dict) -> list[AuditIssue]:
    issues: list[AuditIssue] = []
    for path in REQUIRED_PATHS:
        require((REPO_ROOT / path).exists(), issues, "<paths>", "missing-path", path)
    for path in audit.get("enclosure_jsonl", []):
        require((REPO_ROOT / path).exists(), issues, "<paths>", "missing-enclosure-jsonl", str(path))
    return issues


def validate_rows(audit: dict) -> tuple[list[AuditIssue], int]:
    issues: list[AuditIssue] = []
    rows = audit.get("exact_rows", [])
    if not isinstance(rows, list):
        return [AuditIssue("<rows>", "bad-rows", repr(type(rows)))], 0
    ids = [row.get("id") for row in rows if isinstance(row, dict)]
    require(tuple(ids) == EXPECTED_ROW_IDS, issues, "<rows>", "bad-row-ids", repr(ids))
    ready_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(AuditIssue("<rows>", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        require(row.get("readiness") != "ready_to_apply", issues, row_id, "row-ready-to-apply", repr(row))
        for key in ("role", "formula", "proof_boundary"):
            require(bool(str(row.get(key, "")).strip()), issues, row_id, f"missing-{key}", repr(row.get(key)))
        combined = " ".join(str(row.get(key, "")) for key in ("formula", "reason", "conditional_claim", "proof_boundary")).lower()
        for forbidden in FORBIDDEN_TEXT:
            if forbidden in combined:
                issues.append(AuditIssue(row_id, "forbidden-overclaim", forbidden))
    return issues, ready_count


def validate_note(path: Path) -> list[AuditIssue]:
    if not path.exists():
        return [AuditIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[AuditIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        require(required in text, issues, "note", "missing-text", required)
    lowered = text.lower()
    for forbidden in FORBIDDEN_TEXT:
        if forbidden in lowered:
            issues.append(AuditIssue("note", "forbidden-text", forbidden))
    return issues


def validate(audit_path: Path, note_path: Path) -> tuple[list[AuditIssue], dict]:
    audit = load_json(audit_path)
    issues: list[AuditIssue] = []
    require(
        audit.get("kind") == "jensen_window_pf_factorial_multiplier_split_audit",
        issues,
        "<audit>",
        "bad-kind",
        repr(audit.get("kind")),
    )
    require(
        audit.get("status") == "finite_theorem_search_diagnostic",
        issues,
        "<audit>",
        "bad-status",
        repr(audit.get("status")),
    )
    boundary = str(audit.get("proof_boundary", "")).lower()
    for required in ("conditional", "does not prove", "lambda <= 0"):
        require(required in boundary, issues, "<audit>", "weak-proof-boundary", audit.get("proof_boundary", ""))

    issues.extend(validate_paths(audit))
    row_issues, ready_count = validate_rows(audit)
    issues.extend(row_issues)

    diagnostics = audit.get("finite_diagnostics", {})
    expected_diagnostics = {
        "max_coefficient_index": 64,
        "raw_degree2_rows": 315,
        "raw_degree2_negative_rows": 315,
        "normalized_degree2_rows": 315,
        "normalized_degree2_positive_rows": 315,
        "raw_ratio_below_one_rows": 315,
        "normalized_threshold_rows": 315,
        "normalized_threshold_positive_rows": 315,
    }
    for key, value in expected_diagnostics.items():
        require(diagnostics.get(key) == value, issues, "finite_diagnostics", f"bad-{key}", repr(diagnostics.get(key)))
    require(
        diagnostics.get("lambdas") == ["0.0", "0.000001", "0.0001", "0.01", "0.1"],
        issues,
        "finite_diagnostics",
        "bad-lambdas",
        repr(diagnostics.get("lambdas")),
    )
    for key in ("max_raw_discriminant", "min_normalized_discriminant", "min_normalized_threshold_margin"):
        value = diagnostics.get(key, {})
        require(isinstance(value, dict) and {"sample", "lam", "n"}.issubset(value), issues, "finite_diagnostics", f"bad-{key}", repr(value))

    summary = audit.get("summary", {})
    expected_summary = {
        "exact_rows": 5,
        "raw_degree2_rows": 315,
        "raw_degree2_negative_rows": 315,
        "normalized_degree2_rows": 315,
        "normalized_degree2_positive_rows": 315,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        require(summary.get(key) == value, issues, "summary", f"bad-{key}", repr(summary.get(key)))
    require(summary.get("ready_to_apply_rows") == ready_count, issues, "summary", "ready-count-mismatch", repr(ready_count))

    invariants = "\n".join(str(item) for item in audit.get("invariants", []))
    for required in (
        "No row is ready_to_apply.",
        "The shifted multiplier sequence is used only conditionally.",
        "Raw moment positivity is not promoted to Jensen-window hyperbolicity.",
        "Degree-2 normalized positivity is not promoted to Jensen-window PF-infinity.",
        "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
    ):
        require(required in invariants, issues, "invariants", "missing-invariant", required)

    issues.extend(validate_note(note_path))
    return issues, audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, audit = validate(args.audit, args.note)
    summary = audit.get("summary", {})
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-WINDOW-PF-FACTORIAL-MULTIPLIER-SPLIT {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated Jensen-window PF factorial multiplier split audit: "
            f"{summary.get('exact_rows', 0)} exact rows, "
            f"{summary.get('raw_degree2_negative_rows', 0)} raw degree-2 anti-hyperbolic rows, "
            f"{summary.get('normalized_degree2_positive_rows', 0)} normalized degree-2 positive rows, "
            f"{summary.get('ready_to_apply_rows', 0)} ready-to-apply rows, "
            f"{len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
