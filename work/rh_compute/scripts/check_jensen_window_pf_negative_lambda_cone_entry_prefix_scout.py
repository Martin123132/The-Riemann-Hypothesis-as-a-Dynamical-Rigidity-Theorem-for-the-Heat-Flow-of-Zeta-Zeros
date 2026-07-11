#!/usr/bin/env python3
"""Validate the finite large-negative-lambda ratio-cone prefix scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from jensen_window_pf_negative_lambda_cone_entry_prefix_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    finite_diagnostics,
)


@dataclass(frozen=True)
class NegativeLambdaIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> NegativeLambdaIssue:
    return NegativeLambdaIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[NegativeLambdaIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def repo_ref(path: Path) -> str:
    abs_path = path if path.is_absolute() else REPO_ROOT / path
    return abs_path.relative_to(REPO_ROOT).as_posix()


def validate_top_level(target: dict) -> list[NegativeLambdaIssue]:
    issues: list[NegativeLambdaIssue] = []
    if target.get("kind") != "jensen_window_pf_negative_lambda_cone_entry_prefix_scout":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "finite_negative_lambda_prefix_certificate":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    for key in ("source_target", "source_sign_scout", "enclosure_summary", "generator", "checker"):
        issues.extend(validate_ref("<target>", target.get(key)))
    for ref in target.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<target>", ref))
    for key in ("result_json", "note"):
        if key in target:
            issues.extend(validate_ref("<target>", target.get(key)))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in (
        "finite large-negative-lambda prefix",
        "does not prove an infinite",
        "collared finite",
        "all k",
        "jwpf_06",
        "lambda <= 0 unsettled",
    ):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    return issues


def validate_recomputed_diagnostics(target: dict) -> list[NegativeLambdaIssue]:
    issues: list[NegativeLambdaIssue] = []
    paths = [REPO_ROOT / ref for ref in target.get("enclosure_jsonl", [])]
    prefix = target.get("finite_diagnostics", {})
    lower_upper_k_max = int(prefix.get("lower_upper_k_max", -1))
    monotone_k_max = int(prefix.get("monotone_k_max", -1))
    try:
        diagnostics = finite_diagnostics(paths, lower_upper_k_max, monotone_k_max)
    except Exception as exc:
        return [issue("finite_diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]

    recomputed = asdict(diagnostics)
    for key in (
        "lambdas",
        "coefficient_rows",
        "coefficient_k_min",
        "coefficient_k_max",
        "lower_upper_k_max",
        "monotone_k_max",
        "lower_wall_rows",
        "lower_wall_positive_rows",
        "upper_wall_rows",
        "upper_wall_positive_rows",
        "monotone_gap_rows",
        "monotone_gap_positive_rows",
        "min_lower_wall_margin",
        "min_upper_wall_margin",
        "min_monotone_gap",
    ):
        if prefix.get(key) != recomputed.get(key):
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{prefix.get(key)!r} != {recomputed.get(key)!r}"))

    for row_key, positive_key in (
        ("lower_wall_rows", "lower_wall_positive_rows"),
        ("upper_wall_rows", "upper_wall_positive_rows"),
        ("monotone_gap_rows", "monotone_gap_positive_rows"),
    ):
        if prefix.get(positive_key) != prefix.get(row_key):
            issues.append(
                issue(
                    "finite_diagnostics",
                    f"not-all-positive-{positive_key}",
                    f"{prefix.get(positive_key)!r} != {prefix.get(row_key)!r}",
                )
            )
    return issues


def validate_rows_and_summary(target: dict) -> list[NegativeLambdaIssue]:
    issues: list[NegativeLambdaIssue] = []
    rows = target.get("finite_rows", [])
    if not isinstance(rows, list) or len(rows) != 4:
        issues.append(issue("finite_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else type(rows))))
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            issues.append(issue("finite_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "claim", "rows", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if "positive_rows" in row and row.get("positive_rows") != row.get("rows"):
            issues.append(issue(row_id, "not-all-positive", repr(row)))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "only", "no all-k", "no collared")):
            issues.append(issue(row_id, "weak-row-boundary", boundary))

    summary = target.get("summary", {})
    diagnostics = target.get("finite_diagnostics", {})
    expected_summary = {
        "coefficient_rows": diagnostics.get("coefficient_rows"),
        "lower_wall_rows": diagnostics.get("lower_wall_rows"),
        "lower_wall_positive_rows": diagnostics.get("lower_wall_positive_rows"),
        "upper_wall_rows": diagnostics.get("upper_wall_rows"),
        "upper_wall_positive_rows": diagnostics.get("upper_wall_positive_rows"),
        "monotone_gap_rows": diagnostics.get("monotone_gap_rows"),
        "monotone_gap_positive_rows": diagnostics.get("monotone_gap_positive_rows"),
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if summary.get("ready_to_apply_rows") != 0:
        issues.append(issue("summary", "bad-ready_to_apply_rows", repr(summary.get("ready_to_apply_rows"))))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "bad-target_closing", repr(summary.get("target_closing"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("actual acb-enclosed", "finite", "large-negative-lambda", "infinite or collared finite", "tail control"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("actual a_k", "finite in lambda and k", "ready_to_apply", "target remains open", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path, target_path: Path, target: dict) -> list[NegativeLambdaIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[NegativeLambdaIssue] = []
    diagnostics = target.get("finite_diagnostics", {})
    summary = target.get("summary", {})
    required_strings = [
        "# Jensen-Window PF Negative-Lambda Cone-Entry Prefix Scout",
        "Status: finite negative-lambda prefix certificate",
        "This is not a proof",
        "Artifact kind: `jensen_window_pf_negative_lambda_cone_entry_prefix_scout`",
        target.get("result_json") or repo_ref(target_path),
        str(target.get("enclosure_summary", "")),
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
        (
            "validated Jensen-window PF negative-lambda cone-entry prefix scout: "
            f"{summary.get('coefficient_rows')} coefficient rows, "
            f"{summary.get('lower_wall_rows')} lower-wall rows, "
            f"{summary.get('upper_wall_rows')} upper-wall rows, "
            f"{summary.get('monotone_gap_rows')} monotone-gap rows, 0 issues"
        ),
        f"lambdas: {', '.join(str(item) for item in diagnostics.get('lambdas', []))}",
        (
            "coefficient range: "
            f"A_{diagnostics.get('coefficient_k_min')}..A_{diagnostics.get('coefficient_k_max')}"
        ),
        f"lower and upper walls: 1 <= k <= {diagnostics.get('lower_upper_k_max')}",
        f"monotone wall: 1 <= k <= {diagnostics.get('monotone_k_max')}",
        f"lower wall rows: {diagnostics.get('lower_wall_positive_rows')} / {diagnostics.get('lower_wall_rows')}",
        f"upper wall rows: {diagnostics.get('upper_wall_positive_rows')} / {diagnostics.get('upper_wall_rows')}",
        f"monotone gap rows: {diagnostics.get('monotone_gap_positive_rows')} / {diagnostics.get('monotone_gap_rows')}",
        "finite interval-certified large-negative-lambda prefix",
        "finite-prefix collar theorem",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
    ]
    required_strings.extend(str(ref) for ref in target.get("enclosure_jsonl", []))
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


def validate(target_path: Path, note_path: Path) -> tuple[list[NegativeLambdaIssue], int, int, int, int]:
    target = load_json(target_path)
    issues: list[NegativeLambdaIssue] = []
    issues.extend(validate_top_level(target))
    issues.extend(validate_recomputed_diagnostics(target))
    issues.extend(validate_rows_and_summary(target))
    issues.extend(validate_note(note_path, target_path, target))
    summary = target.get("summary", {})
    return (
        issues,
        int(summary.get("coefficient_rows", 0)),
        int(summary.get("lower_wall_rows", 0)),
        int(summary.get("upper_wall_rows", 0)),
        int(summary.get("monotone_gap_rows", 0)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, coeffs, lower, upper, mono = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "coefficient_rows": coeffs,
                    "lower_wall_rows": lower,
                    "upper_wall_rows": upper,
                    "monotone_gap_rows": mono,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-CONE-PREFIX {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda cone-entry prefix scout: "
            f"{coeffs} coefficient rows, {lower} lower-wall rows, "
            f"{upper} upper-wall rows, {mono} monotone-gap rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
