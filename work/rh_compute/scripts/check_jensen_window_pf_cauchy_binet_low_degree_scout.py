#!/usr/bin/env python3
"""Validate the Jensen-window PF Cauchy-Binet low-degree scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json"
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md"


REQUIRED_FORMULA_COUNT = 15
REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Cauchy-Binet Low-Degree Scout",
    "Status: symbolic theorem-search scout",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_cauchy_binet_low_degree_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py",
    "validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found",
    "adjacent log-concavity",
    "nonnegative Bernstein coefficients",
    "cauchy_binet_identity_found=false",
    "d3_first_negative_contiguous_toeplitz_minor_size=8",
    "d4_first_negative_contiguous_toeplitz_minor_size=6",
    "Frontier Extension",
    "outputs/jensen_window_pf_log_concavity_frontier_scout.md",
    "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py",
    "validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues",
    "degree 3 size 6",
    "degree 4 size 5",
)


@dataclass(frozen=True)
class ScoutIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(row_id: str, name: str, detail: str) -> ScoutIssue:
    return ScoutIssue(row_id=row_id, issue=name, detail=detail)


def algebra_formula_set(algebra: dict) -> set[str]:
    formulas = {str(algebra["degree2"]["jensen_discriminant"])}
    for degree_key in ("degree3", "degree4"):
        for row in algebra[degree_key]["selected_toeplitz_minors"]:
            formulas.add(str(row["determinant"]))
    return formulas


def validate_top_level(scout: dict) -> list[ScoutIssue]:
    issues: list[ScoutIssue] = []
    if scout.get("kind") != "jensen_window_pf_cauchy_binet_low_degree_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    if scout.get("target_ansatz") != "ansatz_02_positive_cauchy_binet_kernel":
        issues.append(issue("<scout>", "bad-target-ansatz", repr(scout.get("target_ansatz"))))
    if scout.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<scout>", "bad-target-obligation", repr(scout.get("target_obligation"))))
    boundary = str(scout.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<scout>", "weak-proof-boundary", scout.get("proof_boundary", "")))

    summary = scout.get("summary", {})
    if summary.get("formula_rows") != REQUIRED_FORMULA_COUNT:
        issues.append(issue("<summary>", "bad-formula-count", repr(summary.get("formula_rows"))))
    if summary.get("bernstein_nonnegative_rows") != REQUIRED_FORMULA_COUNT:
        issues.append(issue("<summary>", "not-all-bernstein-certified", repr(summary.get("bernstein_nonnegative_rows"))))
    if summary.get("cauchy_binet_identity_found") is not False:
        issues.append(issue("<summary>", "cauchy-binet-overclaim", repr(summary.get("cauchy_binet_identity_found"))))
    if summary.get("kernel_identity_found") is not False:
        issues.append(issue("<summary>", "kernel-overclaim", repr(summary.get("kernel_identity_found"))))
    if summary.get("target_closing") is not False:
        issues.append(issue("<summary>", "target-closing-overclaim", repr(summary.get("target_closing"))))
    if "too weak" not in str(summary.get("main_finding", "")).lower():
        issues.append(issue("<summary>", "missing-weakness-finding", str(summary.get("main_finding", ""))))
    return issues


def validate_formula_rows(scout: dict, algebra: dict) -> list[ScoutIssue]:
    rows = scout.get("formula_rows", [])
    issues: list[ScoutIssue] = []
    if not isinstance(rows, list):
        return [issue("formula_rows", "bad-rows", repr(type(rows)))]
    if len(rows) != REQUIRED_FORMULA_COUNT:
        issues.append(issue("formula_rows", "bad-row-count", str(len(rows))))

    seen_formulas: set[str] = set()
    degree2_zero_seen = False
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("formula_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in (
            "determinant",
            "ratio_factorization",
            "positive_monomial_factor",
            "normalized_ratio_polynomial",
            "bernstein_coefficients",
            "bernstein_min_coefficient",
            "proof_boundary",
        ):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        seen_formulas.add(str(row.get("determinant", "")))
        if row.get("bernstein_coefficients_nonnegative") is not True:
            issues.append(issue(row_id, "not-bernstein-nonnegative", repr(row.get("bernstein_coefficients_nonnegative"))))
        if row.get("low_degree_log_concavity_certificate") is not True:
            issues.append(issue(row_id, "missing-log-concavity-certificate", repr(row.get("low_degree_log_concavity_certificate"))))
        if row_id == "degree2_jensen_discriminant":
            if row.get("bernstein_min_coefficient") != "0":
                issues.append(issue(row_id, "degree2-min-should-be-zero", str(row.get("bernstein_min_coefficient"))))
            degree2_zero_seen = True
        elif row.get("bernstein_coefficients_strictly_positive") is not True:
            issues.append(issue(row_id, "expected-strict-positive-coefficients", repr(row.get("bernstein_coefficients_strictly_positive"))))
        boundary = str(row.get("proof_boundary", "")).lower()
        if "not a cauchy-binet kernel identity" not in boundary:
            issues.append(issue(row_id, "weak-row-boundary", str(row.get("proof_boundary", ""))))

    missing_formulas = algebra_formula_set(algebra) - seen_formulas
    if missing_formulas:
        issues.append(issue("formula_rows", "missing-algebra-formulas", repr(sorted(missing_formulas))))
    if not degree2_zero_seen:
        issues.append(issue("degree2_jensen_discriminant", "missing-degree2-row", "degree2 row absent"))
    return issues


def validate_countermodel(scout: dict, algebra: dict) -> list[ScoutIssue]:
    issues: list[ScoutIssue] = []
    row = scout.get("countermodel_ratio_check", {})
    if row.get("ratio_contractions_in_unit_interval") is not True:
        issues.append(issue("countermodel", "ratio-contractors-not-in-box", repr(row.get("ratio_contractions"))))
    if row.get("adjacent_log_concavity_gaps_nonnegative") is not True:
        issues.append(issue("countermodel", "not-log-concave", repr(row.get("adjacent_log_concavity_gaps"))))
    if row.get("selected_low_degree_minors_positive") is not True:
        issues.append(issue("countermodel", "selected-minors-not-positive", repr(row.get("selected_low_degree_minors_positive"))))
    if row.get("d3_first_negative_contiguous_toeplitz_minor_size") != 8:
        issues.append(issue("countermodel", "bad-d3-negative-size", repr(row.get("d3_first_negative_contiguous_toeplitz_minor_size"))))
    if row.get("d4_first_negative_contiguous_toeplitz_minor_size") != 6:
        issues.append(issue("countermodel", "bad-d4-negative-size", repr(row.get("d4_first_negative_contiguous_toeplitz_minor_size"))))
    algebra_counter = algebra.get("finite_countermodel", {})
    if algebra_counter.get("d3_first_negative_contiguous_toeplitz_minor", {}).get("size") != 8:
        issues.append(issue("countermodel", "algebra-d3-size-changed", repr(algebra_counter.get("d3_first_negative_contiguous_toeplitz_minor"))))
    if algebra_counter.get("d4_first_negative_contiguous_toeplitz_minor", {}).get("size") != 6:
        issues.append(issue("countermodel", "algebra-d4-size-changed", repr(algebra_counter.get("d4_first_negative_contiguous_toeplitz_minor"))))
    if "larger contiguous" not in str(row.get("interpretation", "")).lower():
        issues.append(issue("countermodel", "missing-larger-minor-warning", str(row.get("interpretation", ""))))
    return issues


def validate_refs(scout: dict) -> list[ScoutIssue]:
    issues: list[ScoutIssue] = []
    source = scout.get("source_algebra")
    if not isinstance(source, str) or not (REPO_ROOT / source).exists():
        issues.append(issue("<scout>", "missing-source-algebra", repr(source)))
    return issues


def validate_note(path: Path) -> list[ScoutIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ScoutIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "kernel identity is proved", "cauchy-binet identity is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, algebra_path: Path, note_path: Path) -> tuple[list[ScoutIssue], int]:
    scout = load_json(scout_path)
    algebra = load_json(algebra_path)
    issues: list[ScoutIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_refs(scout))
    issues.extend(validate_formula_rows(scout, algebra))
    issues.extend(validate_countermodel(scout, algebra))
    issues.extend(validate_note(note_path))
    return issues, len(scout.get("formula_rows", []))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_SCOUT)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count = validate(args.scout, args.algebra, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "formula_rows": row_count,
                    "kernel_identities_found": 0,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-CB-SCOUT {item.row_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF Cauchy-Binet low-degree scout: "
            f"{row_count} formula rows, {len(issues)} issues, 0 kernel identities found"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
