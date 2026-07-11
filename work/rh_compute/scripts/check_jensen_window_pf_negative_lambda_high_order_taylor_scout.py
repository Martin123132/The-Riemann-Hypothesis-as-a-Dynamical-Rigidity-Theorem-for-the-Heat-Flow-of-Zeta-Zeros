#!/usr/bin/env python3
"""Validate the negative-lambda high-order Taylor truncation scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_high_order_taylor_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
    generate_polynomials,
)


REQUIRED_ROW_IDS = {
    "nhts_01_formal_polynomial_generator",
    "nhts_02_even_coefficient_certification",
    "nhts_03_high_order_truncation_matrix",
    "nhts_04_higher_order_promotion_rejected",
    "nhts_05_remainder_theorem_handoff",
}

ALLOWED_ROLES = {
    "exact_algebra",
    "finite_certificate",
    "finite_diagnostic",
    "rejected_route",
    "conditional_handoff",
}

EXPECTED_COEFFICIENT_SIGNS = {
    0: "positive",
    2: "negative",
    4: "positive",
    6: "negative",
    8: "positive",
    10: "positive",
    12: "negative",
    14: "positive",
}

EXPECTED_STATUSES = {
    6: [
        "invalid_normalizer",
        "invalid_normalizer",
        "invalid_normalizer",
        "invalid_normalizer",
        "overbound",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
    ],
    8: [
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "overbound",
        "invalid_normalizer",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
    ],
    10: [
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "upper_wall_violation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
    ],
    12: [
        "invalid_normalizer",
        "invalid_normalizer",
        "invalid_normalizer",
        "invalid_normalizer",
        "overbound",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
    ],
    14: [
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
        "upper_wall_violation",
        "target_satisfying_truncation",
        "target_satisfying_truncation",
    ],
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda High-Order Taylor Scout",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_high_order_taylor_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_high_order_taylor_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_high_order_taylor_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py",
    "validated Jensen-window PF negative-lambda high-order Taylor scout: 8 coefficient rows, 35 truncation rows, 9 invalid normalizers, 2 upper-wall violations, 3 overbound rows, 0 ready-to-apply rows, 0 issues",
    "(2*q^2*exp(9u)-3*q*exp(5u))*exp(-q*(exp(4u)-1))",
    "c14",
    "truncation degrees 6, 8, 10, 12, and 14",
    "invalid normalizers: 9",
    "upper-wall violations: 2",
    "overbound rows: 3",
    "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
    "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
)


def issue(section: str, name: str, detail: str) -> dict:
    return {"section": section, "issue": name, "detail": detail}


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[dict]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_polynomial_generator() -> list[dict]:
    polynomials = generate_polynomials(14)
    expected = {
        0: {1: Fraction(-3), 2: Fraction(2)},
        2: {1: Fraction(-75, 2), 2: Fraction(165), 3: Fraction(-112), 4: Fraction(16)},
        6: {
            1: Fraction(-3125, 48),
            2: Fraction(87011, 24),
            3: Fraction(-891254, 45),
            4: Fraction(1289582, 45),
            5: Fraction(-46288, 3),
            6: Fraction(31712, 9),
            7: Fraction(-1024, 3),
            8: Fraction(512, 45),
        },
    }
    issues: list[dict] = []
    for degree, expected_poly in expected.items():
        if polynomials[degree] != expected_poly:
            issues.append(issue("polynomial_generator", f"bad-P{degree}", repr(polynomials[degree])))
    return issues


def validate_top_level(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_high_order_taylor_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite_theorem_search_diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_phi_taylor_sign_scout",
        "source_taylor_moment_budget",
        "source_uniform_remainder_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite theorem-search diagnostic",
        "does not bound the infinite taylor remainder",
        "does not prove bounded log-curvature",
        "does not prove cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    params = artifact.get("parameters", {})
    if params.get("max_taylor_degree") != 14:
        issues.append(issue("parameters", "bad-max-degree", repr(params.get("max_taylor_degree"))))
    if params.get("tail_cutoff_n", 0) < 16:
        issues.append(issue("parameters", "weak-tail-cutoff", repr(params.get("tail_cutoff_n"))))
    if params.get("tail_start_k") != 22:
        issues.append(issue("parameters", "bad-tail-start", repr(params.get("tail_start_k"))))
    if params.get("sample_T_values") != [25, 50, 100, 200, 500, 1000, 2000]:
        issues.append(issue("parameters", "bad-t-values", repr(params.get("sample_T_values"))))
    return issues


def validate_recomputed(artifact: dict) -> list[dict]:
    params = artifact.get("parameters", {})
    try:
        diagnostics = asdict(
            build_diagnostics(
                int(params.get("max_taylor_degree", 14)),
                int(params.get("tail_cutoff_n", 16)),
                int(params.get("precision_bits", 256)),
                int(params.get("tail_start_k", 22)),
                list(params.get("sample_T_values", [25, 50, 100, 200, 500, 1000, 2000])),
            )
        )
    except Exception as exc:
        return [issue("diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = artifact.get("diagnostics", {})
    issues: list[dict] = []
    for key, value in diagnostics.items():
        if recorded.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))

    coefficient_signs = {
        int(row.get("degree")): row.get("sign")
        for row in recorded.get("coefficient_rows", [])
        if isinstance(row, dict) and row.get("degree") is not None
    }
    if coefficient_signs != EXPECTED_COEFFICIENT_SIGNS:
        issues.append(issue("coefficient_rows", "bad-signs", repr(coefficient_signs)))

    by_degree: dict[int, list[str]] = {}
    for row in recorded.get("truncation_rows", []):
        if isinstance(row, dict):
            by_degree.setdefault(int(row.get("truncation_degree")), []).append(str(row.get("status")))
    if by_degree != EXPECTED_STATUSES:
        issues.append(issue("truncation_rows", "bad-status-matrix", repr(by_degree)))
    return issues


def validate_rows(artifact: dict) -> tuple[list[dict], int, int]:
    rows = artifact.get("scout_rows", [])
    issues: list[dict] = []
    if not isinstance(rows, list):
        return [issue("scout_rows", "bad-rows", repr(type(rows)))], 0, 0
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        issues.append(issue("scout_rows", "bad-row-ids", repr(sorted(ids))))
    ready_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("scout_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
            issues.append(issue(row_id, "ready-to-apply-overclaim", row_id))
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        text = " ".join(str(row.get(key, "")) for key in ("claim", "proof_boundary")).lower()
        if row.get("role") == "rejected_route" and "reject" not in text:
            issues.append(issue(row_id, "missing-rejection-language", text))
        if row.get("role") == "conditional_handoff" and "not" not in text:
            issues.append(issue(row_id, "weak-boundary", text))
    return issues, len(rows), ready_count


def validate_summary(artifact: dict, row_count: int, ready_count: int) -> list[dict]:
    summary = artifact.get("summary", {})
    issues: list[dict] = []
    expected = {
        "scout_rows": 5,
        "coefficient_rows": 8,
        "truncation_rows": 35,
        "invalid_normalizer_rows": 9,
        "upper_wall_violation_rows": 2,
        "overbound_rows": 3,
        "target_satisfying_rows": 21,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 5:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("high-order", "degree 14", "invalid normalizers", "upper-wall", "infinite-tail remainder"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("no row", "finite taylor truncations", "actual zeta finite prefix", "remainder theorem remains open", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[dict]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[dict] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "bounded log-curvature theorem is proved",
        "cone entry is proved",
        "infinite taylor-tail theorem is proved",
        "higher finite taylor truncation proves",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(artifact_path: Path, note_path: Path) -> tuple[list[dict], dict]:
    artifact = load_json(artifact_path)
    issues: list[dict] = []
    issues.extend(validate_polynomial_generator())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, ready_count = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, ready_count))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, summary = validate(args.artifact, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": issues}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-HIGH-TAYLOR {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda high-order Taylor scout: "
            f"{summary.get('coefficient_rows')} coefficient rows, "
            f"{summary.get('truncation_rows')} truncation rows, "
            f"{summary.get('invalid_normalizer_rows')} invalid normalizers, "
            f"{summary.get('upper_wall_violation_rows')} upper-wall violations, "
            f"{summary.get('overbound_rows')} overbound rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
