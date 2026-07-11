#!/usr/bin/env python3
"""Validate the negative-lambda k300 precision-repair audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_k300_precision_repair_audit import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlk300pra_01_broad_k300_input",
    "nlk300pra_02_broad_run_failure_gate",
    "nlk300pra_03_local_precision_repair",
    "nlk300pra_04_repaired_raw_wall_stress",
    "nlk300pra_05_repaired_decrement_corridor_stress",
    "nlk300pra_06_repaired_theta_shape_stress",
    "nlk300pra_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "finite_input",
    "precision_gate",
    "finite_stress",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda k300 Precision-Repair Audit",
    "Status: finite precision-repair theorem-search diagnostic",
    "This is not",
    "Artifact kind: `jensen_window_pf_negative_lambda_k300_precision_repair_audit`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_k300_precision_repair_audit.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",
    "validated Jensen-window PF negative-lambda k300 precision-repair audit: 7 rows, 0 issues, 894 repaired decrement-corridor rows, 891 repaired theta-k-monotone rows, 0 ready-to-apply rows",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k300.jsonl",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k220_k250_dps220.jsonl",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
    "first raw decrease failure: k=252",
    "first lower decrement failure: k=252",
    "first upper decrement failure: k=253",
    "first raw upper-wall failure: k=271",
    "first theta-k monotone failure: k=233",
    "raw lower wall rows: 897 / 897",
    "raw upper wall rows: 897 / 897",
    "raw decrease rows: 894 / 894",
    "lower decrement rows: 894 / 894",
    "upper decrement rows: 894 / 894",
    "decrement-corridor rows: 894 / 894",
    "theta unit rows: 894 / 894",
    "theta-k monotone rows: 891 / 891",
    "theta lambda-order rows: 596 / 596",
    "precision alarms, not mathematical",
    "finite stress evidence, not an all-k theorem",
)


@dataclass(frozen=True)
class K300AuditIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> K300AuditIssue:
    return K300AuditIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[K300AuditIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[K300AuditIssue]:
    issues: list[K300AuditIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_k300_precision_repair_audit":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite precision-repair theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "base_enclosure_jsonl",
        "base_summary",
        "source_decrement_scout",
        "source_raw_corridor_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for key in ("repair_enclosure_jsonl", "repair_summaries"):
        values = artifact.get(key, [])
        if not isinstance(values, list) or not values:
            issues.append(issue("<artifact>", f"missing-{key}", repr(values)))
            continue
        for ref in values:
            issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite", "precision", "does not prove", "raw-corridor", "cone entry", "jwpf_06", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[K300AuditIssue]:
    recomputed = build_artifact()
    issues: list[K300AuditIssue] = []
    for key in ("audit_rows", "base_diagnostics", "repaired_diagnostics", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[K300AuditIssue], int, int]:
    rows = artifact.get("audit_rows", [])
    issues: list[K300AuditIssue] = []
    if not isinstance(rows, list):
        return [issue("audit_rows", "bad-rows", repr(type(rows)))], 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))

    ready_to_apply = 0
    precision_gates = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("audit_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "precision_gate":
            precision_gates += 1
        for key in ("source", "summary"):
            if key in row:
                issues.extend(validate_ref(row_id, row[key]))
        for key in ("sources", "summaries"):
            for ref in row.get(key, []):
                issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "not", "only", "precision", "hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), ready_to_apply


def validate_summary(artifact: dict, row_count: int, ready_to_apply: int) -> list[K300AuditIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "audit_rows": 7,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "base_rows": 903,
        "base_k_max": 300,
        "repair_rows": 107,
        "repair_ranges": [[220, 250], [245, 320]],
        "repaired_raw_total_rows": 897,
        "repaired_adjacent_total_rows": 894,
        "repaired_theta_k_total_rows": 891,
        "repaired_raw_lower_rows": 897,
        "repaired_raw_upper_rows": 897,
        "repaired_raw_decrease_rows": 894,
        "repaired_lower_decrement_rows": 894,
        "repaired_upper_decrement_rows": 894,
        "repaired_decrement_corridor_rows": 894,
        "repaired_theta_unit_rows": 894,
        "repaired_theta_k_monotone_rows": 891,
        "repaired_theta_lambda_order_rows": 596,
        "repaired_theta_lambda_order_total_rows": 596,
        "open_theorem_target": False,
    }
    issues: list[K300AuditIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 7:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_to_apply != 0:
        issues.append(issue("summary", "bad-ready-to-apply-count", str(ready_to_apply)))
    firsts = summary.get("base_lambda_minus_100_first_failures", {})
    expected_firsts = {
        "raw_decrease": 252,
        "lower_decrement": 252,
        "upper_decrement": 253,
        "raw_upper": 271,
        "raw_lower": 272,
        "theta_k_monotone": 233,
    }
    for key, value in expected_firsts.items():
        if firsts.get(key) != value:
            issues.append(issue("summary", f"bad-first-{key}", repr(firsts)))
    repaired = artifact.get("repaired_diagnostics", {})
    for label, row in repaired.get("per_lambda", {}).items():
        if row.get("first_failures") != {}:
            issues.append(issue("repaired", "unexpected-first-failure", f"{label}: {row.get('first_failures')}"))
    if repaired.get("theta_lambda_first_failure") is not None:
        issues.append(issue("repaired", "unexpected-theta-lambda-failure", repr(repaired.get("theta_lambda_first_failure"))))
    base = artifact.get("base_diagnostics", {})
    if base.get("theta_lambda_first_failure") != {"left": "-50.0", "right": "-100.0", "k": 249}:
        issues.append(issue("base", "bad-theta-lambda-first-failure", repr(base.get("theta_lambda_first_failure"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "k300 stress extension",
        "local precision repair",
        "broad dps160",
        "lambda=-100",
        "897/897 raw wall",
        "894/894 decrement-corridor",
        "891/891 theta-k",
        "596/596 theta lambda-order",
        "finite stress evidence",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "precision alarms", "finite evidence", "open zeta-specific", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[K300AuditIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[K300AuditIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "raw-corridor theorem is proved",
        "decrement route is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[K300AuditIssue], dict]:
    artifact = load_json(target_path)
    issues: list[K300AuditIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, ready_to_apply = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, ready_to_apply))
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
            print(f"JWPF-NEG-LAMBDA-K300-PRECISION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda k300 precision-repair audit: "
            f"{summary.get('audit_rows')} rows, {len(issues)} issues, "
            f"{summary.get('repaired_decrement_corridor_rows')} repaired decrement-corridor rows, "
            f"{summary.get('repaired_theta_k_monotone_rows')} repaired theta-k-monotone rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
