#!/usr/bin/env python3
"""Validate the worst-row Laguerre root-bracket certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate import (  # noqa: E402
    DEFAULT_ARB_PRECISION_BITS,
    DEFAULT_BISECTION_STEPS,
    DEFAULT_DECIMAL_PRECISION,
    DEFAULT_GRID_JSON,
    DEFAULT_INDEX,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_ORDER,
    DEFAULT_OUT_JSON,
    DEFAULT_QUADRATURE_JSON,
    DEFAULT_T,
    REPO_ROOT,
    build_artifact,
    laguerre_value,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwrbr_01_worst_row_setup",
    "nlrgwrbr_02_arb_root_brackets",
    "nlrgwrbr_03_node_interval_handoff",
    "nlrgwrbr_04_weight_underflow_diagnostic",
    "nlrgwrbr_05_weight_interval_target",
    "nlrgwrbr_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "arb_node_certificate",
    "node_interval_handoff",
    "floating_diagnostic",
    "open_requirement",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_exact",
    "available_for_intervalization",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Laguerre Root-Bracket Certificate",
    "Status: worst-row Laguerre root-bracket certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: 6 rows, 0 issues, 320 root brackets, 30 zero floating weights, 2 intervalization rows, 0 ready-to-apply rows",
    "widest bracket: 4.145787935554796563E-17 at root 320",
    "narrowest bracket: 1.710100983788340277E-19 at root 2",
    "zero floating weights: 30",
    "first zero floating weight index: 291",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class RootBracketIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> RootBracketIssue:
    return RootBracketIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[RootBracketIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def arb_from_decimal_text(text: str) -> flint.arb:
    return flint.arb(text)


def sign_name(value: flint.arb) -> str:
    if bool(value > 0 and not value.contains(0)):
        return "positive"
    if bool(value < 0 and not value.contains(0)):
        return "negative"
    return "unresolved"


def validate_top_level(artifact: dict) -> list[RootBracketIssue]:
    issues: list[RootBracketIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row Laguerre root-bracket certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_cancellation_reduced_grid_scout",
        "source_cancellation_reduced_grid_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_quadrature_ladder_scout",
        "source_quadrature_ladder_json",
        "source_node_c0_certificate",
        "source_phi_tail_grid_certificate",
        "source_coefficient_core_certificate",
        "source_first_omitted_denominator_certificate",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "root-bracket certificate only",
        "does not certify christoffel",
        "does not evaluate phi",
        "quadrature-remainder",
        "all recorded rows",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[RootBracketIssue]:
    try:
        recomputed = build_artifact(DEFAULT_GRID_JSON, DEFAULT_INTERVAL_JSON, DEFAULT_QUADRATURE_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[RootBracketIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[RootBracketIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[RootBracketIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    available = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row or row[key] in (None, ""):
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") not in ALLOWED_READINESS:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "available_for_intervalization":
            available += 1
        if row.get("readiness") == "ready_to_apply":
            issues.append(issue(row_id, "forbidden-ready-to-apply", row_id))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "open")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    return issues, len(rows), available


def validate_root_brackets(artifact: dict) -> list[RootBracketIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    diagnostics = rows_by_id.get("nlrgwrbr_02_arb_root_brackets", {}).get("diagnostics", {})
    root_summary = diagnostics.get("root_bracket_summary", {})
    root_rows = diagnostics.get("root_bracket_rows", [])
    issues: list[RootBracketIssue] = []
    expected_summary = {
        "quadrature_order": DEFAULT_ORDER,
        "index": DEFAULT_INDEX,
        "T": DEFAULT_T,
        "alpha": "41/2",
        "root_bracket_rows": DEFAULT_ORDER,
        "bisection_steps": DEFAULT_BISECTION_STEPS,
        "decimal_precision": DEFAULT_DECIMAL_PRECISION,
        "arb_precision_bits": DEFAULT_ARB_PRECISION_BITS,
        "all_endpoint_signs_certified": True,
        "all_brackets_sign_changing": True,
        "all_brackets_ordered_and_disjoint": True,
        "widest_bracket_root_index": 320,
        "widest_bracket_width": "4.145787935554796563E-17",
        "narrowest_bracket_root_index": 2,
        "narrowest_bracket_width": "1.710100983788340277E-19",
        "node_brackets_available_for_worst_row": True,
    }
    for key, expected in expected_summary.items():
        if root_summary.get(key) != expected:
            issues.append(issue("root_summary", f"bad-{key}", f"{root_summary.get(key)!r} != {expected!r}"))
    if len(root_rows) != DEFAULT_ORDER:
        issues.append(issue("root_rows", "bad-row-count", str(len(root_rows))))
        return issues

    flint.ctx.prec = DEFAULT_ARB_PRECISION_BITS
    alpha = flint.arb(2 * DEFAULT_INDEX - 1) / flint.arb(2)
    previous_right: Decimal | None = None
    max_width = Decimal(0)
    min_width: Decimal | None = None
    for expected_index, row in enumerate(root_rows, start=1):
        if not isinstance(row, dict):
            issues.append(issue("root_rows", "bad-row", repr(row)))
            continue
        row_id = f"root_{expected_index}"
        if row.get("root_index") != expected_index:
            issues.append(issue(row_id, "bad-root-index", repr(row.get("root_index"))))
        left = dec(row.get("left_endpoint", "0"))
        right = dec(row.get("right_endpoint", "0"))
        width = dec(row.get("width", "0"))
        if not left < right:
            issues.append(issue(row_id, "bad-order", f"{left} !< {right}"))
        if width != right - left:
            issues.append(issue(row_id, "bad-width", f"{width} != {right-left}"))
        if previous_right is not None and not previous_right < left:
            issues.append(issue(row_id, "brackets-not-disjoint", f"{previous_right} !< {left}"))
        previous_right = right
        max_width = max(max_width, width)
        min_width = width if min_width is None else min(min_width, width)
        if width >= Decimal("5E-17"):
            issues.append(issue(row_id, "width-too-large", str(width)))
        left_value = laguerre_value(DEFAULT_ORDER, alpha, arb_from_decimal_text(row["left_endpoint"]))
        right_value = laguerre_value(DEFAULT_ORDER, alpha, arb_from_decimal_text(row["right_endpoint"]))
        left_sign = sign_name(left_value)
        right_sign = sign_name(right_value)
        if left_sign != row.get("left_sign"):
            issues.append(issue(row_id, "bad-left-sign", f"{left_sign} != {row.get('left_sign')!r}"))
        if right_sign != row.get("right_sign"):
            issues.append(issue(row_id, "bad-right-sign", f"{right_sign} != {row.get('right_sign')!r}"))
        if left_sign == "unresolved" or right_sign == "unresolved" or left_sign == right_sign:
            issues.append(issue(row_id, "not-sign-changing", f"{left_sign}, {right_sign}"))
        if dec(row.get("v_interval_left_for_T_10000", "0")) != left / Decimal(DEFAULT_T):
            issues.append(issue(row_id, "bad-v-left", repr(row.get("v_interval_left_for_T_10000"))))
        if dec(row.get("v_interval_right_for_T_10000", "0")) != right / Decimal(DEFAULT_T):
            issues.append(issue(row_id, "bad-v-right", repr(row.get("v_interval_right_for_T_10000"))))
    if max_width != Decimal("4.14578793555479656329776361189942690543830394744873046875E-17"):
        issues.append(issue("root_rows", "unexpected-max-width", str(max_width)))
    if min_width != Decimal("1.7101009837883402765756901686700075515545904636383056640625E-19"):
        issues.append(issue("root_rows", "unexpected-min-width", str(min_width)))
    return issues


def validate_weight_underflow(artifact: dict) -> list[RootBracketIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    diagnostics = rows_by_id.get("nlrgwrbr_04_weight_underflow_diagnostic", {}).get("diagnostics", {})
    issues: list[RootBracketIssue] = []
    expected = {
        "quadrature_order": 320,
        "alpha_float": 20.5,
        "floating_weight_count": 320,
        "zero_weight_count": 30,
        "first_zero_weight_index": 291,
        "last_zero_weight_index": 320,
        "last_positive_weight_index": 290,
        "last_positive_weight_float": "9.303e-321",
        "first_zero_weight_node_float": "886.9174722103623",
    }
    for key, value in expected.items():
        if diagnostics.get(key) != value:
            issues.append(issue("weight_underflow", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    boundary = str(diagnostics.get("proof_boundary", "")).lower()
    for required in ("floating underflow", "not interval"):
        if required not in boundary:
            issues.append(issue("weight_underflow", "weak-proof-boundary", required))
    return issues


def validate_summary(artifact: dict, row_count: int, available: int) -> list[RootBracketIssue]:
    summary = artifact.get("summary", {})
    issues: list[RootBracketIssue] = []
    expected = {
        "matrix_rows": 6,
        "root_bracket_rows": 320,
        "quadrature_order": 320,
        "index": 21,
        "T": 10000,
        "alpha": "41/2",
        "bisection_steps": 60,
        "widest_bracket_width": "4.145787935554796563E-17",
        "widest_bracket_root_index": 320,
        "narrowest_bracket_width": "1.710100983788340277E-19",
        "narrowest_bracket_root_index": 2,
        "zero_float_weight_count": 30,
        "first_zero_float_weight_index": 291,
        "available_for_intervalization_rows": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if available != summary.get("available_for_intervalization_rows"):
        issues.append(
            issue("summary", "available-count-mismatch", f"{available} != {summary.get('available_for_intervalization_rows')!r}")
        )
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("certified arb brackets", "320 nodes", "underflow", "christoffel weights"):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    return issues


def validate_note(path: Path) -> list[RootBracketIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RootBracketIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    forbidden_phrases = (
        "christoffel-weight interval certificate is proved",
        "quadrature-remainder theorem is proved",
        "finite-grid interval certificate is complete",
        "rh is proved",
        "lambda <= 0 is proved",
    )
    lowered = text.lower()
    for phrase in forbidden_phrases:
        if phrase in lowered:
            issues.append(issue("note", "forbidden-promotion-language", phrase))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    getcontext().prec = DEFAULT_DECIMAL_PRECISION
    args = build_parser().parse_args()
    artifact_path = args.artifact if args.artifact.is_absolute() else REPO_ROOT / args.artifact
    note_path = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = load_json(artifact_path)
    issues: list[RootBracketIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, available = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_root_brackets(artifact))
    issues.extend(validate_weight_underflow(artifact))
    issues.extend(validate_summary(artifact, row_count, available))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"LAGUERRE-ROOT {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
