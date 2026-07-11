#!/usr/bin/env python3
"""Validate the finite-collar contract for the negative-lambda cone prefix."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from jensen_window_pf_negative_lambda_finite_collar_contract import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(section: str, name: str, detail: str) -> dict:
    return {"section": section, "issue": name, "detail": detail}


def validate_ref(section: str, ref: object) -> list[dict]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def repo_ref(path: Path) -> str:
    abs_path = path if path.is_absolute() else REPO_ROOT / path
    return abs_path.relative_to(REPO_ROOT).as_posix()


def validate_top_level(contract: dict) -> list[dict]:
    issues: list[dict] = []
    if contract.get("kind") != "jensen_window_pf_negative_lambda_finite_collar_contract":
        issues.append(issue("<contract>", "bad-kind", repr(contract.get("kind"))))
    if contract.get("status") != "finite_collar_contract":
        issues.append(issue("<contract>", "bad-status", repr(contract.get("status"))))
    for key in ("source_prefix_scout", "source_ratio_cone_lemma", "source_cone_entry_target", "prefix_json", "generator", "checker"):
        issues.extend(validate_ref("<contract>", contract.get(key)))
    for ref in contract.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<contract>", ref))
    for key in ("result_json", "note"):
        if key in contract:
            issues.extend(validate_ref("<contract>", contract.get(key)))
    boundary = str(contract.get("proof_boundary", "")).lower()
    for required in (
        "finite-collar accounting",
        "active finite depth",
        "does not prove an infinite",
        "tail theorem",
        "collared finite flow theorem beyond the stated depth",
        "jwpf_06",
        "lambda <= 0 unsettled",
    ):
        if required not in boundary:
            issues.append(issue("<contract>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(contract: dict) -> list[dict]:
    issues: list[dict] = []
    paths = [REPO_ROOT / ref for ref in contract.get("enclosure_jsonl", [])]
    try:
        diagnostics = asdict(build_diagnostics(paths))
    except Exception as exc:
        return [issue("finite_diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = contract.get("finite_diagnostics", {})
    for key, value in diagnostics.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))

    for row_key, positive_key in (
        ("active_lower_rows", "active_lower_positive_rows"),
        ("active_upper_rows", "active_upper_positive_rows"),
        ("active_monotone_rows", "active_monotone_positive_rows"),
        ("collar_lower_upper_rows", "collar_lower_upper_positive_rows"),
        ("collar_monotone_rows", "collar_monotone_positive_rows"),
    ):
        if recorded.get(positive_key) != recorded.get(row_key):
            issues.append(
                issue(
                    "finite_diagnostics",
                    f"not-all-positive-{positive_key}",
                    f"{recorded.get(positive_key)!r} != {recorded.get(row_key)!r}",
                )
            )
    return issues


def validate_rows_and_summary(contract: dict) -> list[dict]:
    issues: list[dict] = []
    rows = contract.get("contract_rows", [])
    if not isinstance(rows, list) or len(rows) != 4:
        issues.append(issue("contract_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else type(rows))))
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            issues.append(issue("contract_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("requirement", "finite", "only", "target", "not")):
            issues.append(issue(row_id, "weak-row-boundary", boundary))

    summary = contract.get("summary", {})
    diagnostics = contract.get("finite_diagnostics", {})
    expected_summary = {
        "certified_active_depth": diagnostics.get("certified_active_depth"),
        "first_collar_x": diagnostics.get("first_collar_x"),
        "second_collar_x": diagnostics.get("second_collar_x"),
        "active_lower_rows": diagnostics.get("active_lower_rows"),
        "active_upper_rows": diagnostics.get("active_upper_rows"),
        "active_monotone_rows": diagnostics.get("active_monotone_rows"),
        "collar_lower_upper_rows": diagnostics.get("collar_lower_upper_rows"),
        "collar_monotone_rows": diagnostics.get("collar_monotone_rows"),
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if summary.get("ready_to_apply_rows") != 0:
        issues.append(issue("summary", "bad-ready_to_apply_rows", repr(summary.get("ready_to_apply_rows"))))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "bad-target_closing", repr(summary.get("target_closing"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        f"k={diagnostics.get('certified_active_depth')}",
        f"x_{diagnostics.get('first_collar_x')}",
        f"x_{diagnostics.get('second_collar_x')}",
        "local seed",
        "all-k tail theorem",
        "finite-collar flow theorem",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in contract.get("invariants", [])).lower()
    for required in ("finite-collar requirement", "finite", "ready_to_apply", "tail", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path, contract_path: Path, contract: dict) -> list[dict]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[dict] = []
    diagnostics = contract.get("finite_diagnostics", {})
    summary = contract.get("summary", {})
    required_strings = [
        "# Jensen-Window PF Negative-Lambda Finite-Collar Contract",
        "Status: finite-collar accounting contract",
        "This is not a proof",
        "Artifact kind: `jensen_window_pf_negative_lambda_finite_collar_contract`",
        contract.get("result_json") or repo_ref(contract_path),
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_finite_collar_contract.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
        (
            "validated Jensen-window PF negative-lambda finite-collar contract: "
            f"active depth K={summary.get('certified_active_depth')}, "
            f"{summary.get('active_lower_rows')} active lower rows, "
            f"{summary.get('active_upper_rows')} active upper rows, "
            f"{summary.get('active_monotone_rows')} active monotone rows, "
            f"{summary.get('collar_lower_upper_rows')} collar lower/upper rows, "
            f"{summary.get('collar_monotone_rows')} collar monotone rows, 0 issues"
        ),
        "checking x_1..x_K requires collar variables x_{K+1}",
        "and x_{K+2} for monotone wall tests",
        f"certified active depth: K={diagnostics.get('certified_active_depth')}",
        f"first collar: x_{diagnostics.get('first_collar_x')}",
        f"second collar: x_{diagnostics.get('second_collar_x')}",
        f"available coefficient range: A_0..A_{diagnostics.get('coefficient_k_max')}",
        f"available ratio range: x_1..x_{diagnostics.get('available_x_max')}",
        f"lower/upper range: x_1..x_{diagnostics.get('lower_upper_k_max')}",
        f"monotone-gap range: 1 <= k <= {diagnostics.get('monotone_gap_k_max')}",
        f"active lower rows: {diagnostics.get('active_lower_positive_rows')} / {diagnostics.get('active_lower_rows')}",
        f"active upper rows: {diagnostics.get('active_upper_positive_rows')} / {diagnostics.get('active_upper_rows')}",
        f"active monotone rows: {diagnostics.get('active_monotone_positive_rows')} / {diagnostics.get('active_monotone_rows')}",
        f"collar lower/upper rows: {diagnostics.get('collar_lower_upper_positive_rows')} / {diagnostics.get('collar_lower_upper_rows')}",
        f"collar monotone rows: {diagnostics.get('collar_monotone_positive_rows')} / {diagnostics.get('collar_monotone_rows')}",
        "It does not seed",
        "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
        "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
    ]
    for required in required_strings:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "cone entry is proved",
        "monotone contractions are proved for all zeta windows",
        "analytic monotone-contraction theorem is proved",
        "closed differential inequality is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(contract_path: Path, note_path: Path) -> tuple[list[dict], dict]:
    contract = load_json(contract_path)
    issues: list[dict] = []
    issues.extend(validate_top_level(contract))
    issues.extend(validate_recomputed(contract))
    issues.extend(validate_rows_and_summary(contract))
    issues.extend(validate_note(note_path, contract_path, contract))
    return issues, contract.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, summary = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "summary": summary,
                    "issues": issues,
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-FINITE-COLLAR {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda finite-collar contract: "
            f"active depth K={summary.get('certified_active_depth')}, "
            f"{summary.get('active_lower_rows')} active lower rows, "
            f"{summary.get('active_upper_rows')} active upper rows, "
            f"{summary.get('active_monotone_rows')} active monotone rows, "
            f"{summary.get('collar_lower_upper_rows')} collar lower/upper rows, "
            f"{summary.get('collar_monotone_rows')} collar monotone rows, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
