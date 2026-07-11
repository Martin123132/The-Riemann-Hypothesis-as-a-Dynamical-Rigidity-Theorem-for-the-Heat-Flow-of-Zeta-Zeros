#!/usr/bin/env python3
"""Validate the Jensen-window PF monotone-contraction frontier scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_frontier_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md"

EXPECTED_EXACT_ROWS = {
    "mcfs_01_d3_m8_monotone_contractions": {
        "degree": 3,
        "minor_size": 8,
        "source_frontier_row": "d3_column_recurrence_m8",
        "bernstein_multidegree": [7, 2],
        "bernstein_coefficient_count": 24,
        "bernstein_min_coefficient": "45",
    },
    "mcfs_02_d4_m6_monotone_contractions": {
        "degree": 4,
        "minor_size": 6,
        "source_frontier_row": "d4_column_recurrence_m6",
        "bernstein_multidegree": [7, 3, 1],
        "bernstein_coefficient_count": 64,
        "bernstein_min_coefficient": "84",
    },
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Monotone-Contraction Frontier Scout",
    "Status: symbolic frontier scout",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_monotone_contraction_frontier_scout`",
    "work/rh_compute/results/jensen_window_pf_monotone_contraction_frontier_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_frontier_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_frontier_scout.py",
    "validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues",
    "x1 <= x2 <= x3",
    "d3_column_recurrence_m8",
    "d4_column_recurrence_m6",
    "Bernstein minimum 45",
    "Bernstein minimum 84",
    "65/66",
    "44/65",
    "50/77",
    "210 finite zeta rows",
    "outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md",
    "outputs/jensen_window_pf_column_recurrence_contract.md",
)


@dataclass(frozen=True)
class MonotoneScoutIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> MonotoneScoutIssue:
    return MonotoneScoutIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[MonotoneScoutIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if ref.startswith(("http://", "https://")):
        return []
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(scout: dict) -> list[MonotoneScoutIssue]:
    issues: list[MonotoneScoutIssue] = []
    if scout.get("kind") != "jensen_window_pf_monotone_contraction_frontier_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    if scout.get("status") != "symbolic_frontier_scout":
        issues.append(issue("<scout>", "bad-status", repr(scout.get("status"))))
    if scout.get("target_frontier_rows") != ["d3_column_recurrence_m8", "d4_column_recurrence_m6"]:
        issues.append(issue("<scout>", "bad-target-frontier-rows", repr(scout.get("target_frontier_rows"))))
    boundary = str(scout.get("proof_boundary", "")).lower()
    for required in ("frontier scout", "does not prove", "all-shape", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<scout>", "weak-proof-boundary", required))
    for key in ("source_column_recurrence_contract", "source_frontier_matrix", "source_algebra"):
        issues.extend(validate_ref("<scout>", scout.get(key)))
    definition = scout.get("ratio_contraction_definition", {})
    if definition.get("extra_sufficient_region") != "0 <= x1 <= x2 <= x3 <= 1, truncated to available degree":
        issues.append(issue("ratio_contraction_definition", "bad-extra-region", repr(definition.get("extra_sufficient_region"))))
    return issues


def validate_exact_rows(scout: dict) -> list[MonotoneScoutIssue]:
    rows = scout.get("exact_frontier_rows", [])
    issues: list[MonotoneScoutIssue] = []
    if not isinstance(rows, list):
        return [issue("exact_frontier_rows", "bad-rows", repr(type(rows)))]
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for row_id, expected in EXPECTED_EXACT_ROWS.items():
        row = by_id.get(row_id)
        if row is None:
            issues.append(issue(row_id, "missing-row", row_id))
            continue
        for key, value in expected.items():
            if row.get(key) != value:
                issues.append(issue(row_id, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        if row.get("bernstein_negative_coefficient_count") != 0:
            issues.append(issue(row_id, "negative-bernstein-coefficient", repr(row.get("bernstein_negative_coefficient_count"))))
        coefficients = row.get("bernstein_coefficients", [])
        if not isinstance(coefficients, list) or len(coefficients) != expected["bernstein_coefficient_count"]:
            issues.append(issue(row_id, "bad-coefficient-list", repr(len(coefficients) if isinstance(coefficients, list) else type(coefficients))))
        text = " ".join(str(row.get(key, "")) for key in ("certificate", "proof_boundary", "monotone_substitution")).lower()
        for required in ("bernstein", "monotone", "sufficient", "does not prove"):
            if required not in text:
                issues.append(issue(row_id, "missing-row-boundary-text", required))
    return issues


def validate_countermodel(scout: dict) -> list[MonotoneScoutIssue]:
    counter = scout.get("countermodel_check", {})
    issues: list[MonotoneScoutIssue] = []
    expected_contractions = {"x1": "65/66", "x2": "44/65", "x3": "50/77"}
    if counter.get("ratio_contractions") != expected_contractions:
        issues.append(issue("countermodel", "bad-contractions", repr(counter.get("ratio_contractions"))))
    if counter.get("monotone_contractions") is not False:
        issues.append(issue("countermodel", "countermodel-should-violate", repr(counter.get("monotone_contractions"))))
    for required in ("x2 < x1", "x3 < x2"):
        if required not in counter.get("violations", []):
            issues.append(issue("countermodel", "missing-violation", required))
    if counter.get("d3_m8_normalized_value") != "-122505653/6324912":
        issues.append(issue("countermodel", "bad-d3-value", repr(counter.get("d3_m8_normalized_value"))))
    if counter.get("d4_m6_normalized_value") != "-19180303/754677":
        issues.append(issue("countermodel", "bad-d4-value", repr(counter.get("d4_m6_normalized_value"))))
    return issues


def validate_finite_zeta(scout: dict) -> list[MonotoneScoutIssue]:
    rows = scout.get("finite_zeta_arb_diagnostics", [])
    issues: list[MonotoneScoutIssue] = []
    if not isinstance(rows, list):
        return [issue("finite_zeta_arb_diagnostics", "bad-rows", repr(type(rows)))]
    by_degree = {row.get("degree"): row for row in rows if isinstance(row, dict)}
    for degree in (3, 4):
        row = by_degree.get(degree)
        if row is None:
            issues.append(issue("finite_zeta", "missing-degree", str(degree)))
            continue
        for key in ("checked_rows", "monotone_contraction_rows", "positive_hard_value_rows"):
            if row.get(key) != 105:
                issues.append(issue(f"finite_zeta_d{degree}", f"bad-{key}", repr(row.get(key))))
        if row.get("failed_or_inconclusive_rows") != 0:
            issues.append(issue(f"finite_zeta_d{degree}", "failed-or-inconclusive", repr(row.get("failed_or_inconclusive_rows"))))
        if "x1" not in row.get("contraction_ranges", {}):
            issues.append(issue(f"finite_zeta_d{degree}", "missing-contraction-range", repr(row.get("contraction_ranges"))))
        for key in ("min_gap_sample", "min_normalized_hard_value_sample"):
            if not str(row.get(key, "")).strip():
                issues.append(issue(f"finite_zeta_d{degree}", f"missing-{key}", repr(row.get(key))))
    return issues


def validate_summary(scout: dict) -> list[MonotoneScoutIssue]:
    summary = scout.get("summary", {})
    issues: list[MonotoneScoutIssue] = []
    expected = {
        "exact_certificate_rows": 2,
        "total_bernstein_coefficients": 88,
        "negative_bernstein_coefficients": 0,
        "finite_zeta_checked_rows": 210,
        "finite_zeta_monotone_rows": 210,
        "finite_zeta_positive_hard_rows": 210,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("monotone-contraction", "countermodel", "zeta-grid", "remaining theorem target"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in scout.get("invariants", [])).lower()
    for required in ("positive bernstein", "countermodel", "finite zeta", "target_closing"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[MonotoneScoutIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[MonotoneScoutIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "monotone contractions are proved for all zeta windows",
        "all-shape positivity is proved",
        "cauchy-binet identity is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, note_path: Path) -> tuple[list[MonotoneScoutIssue], int, int, int]:
    scout = load_json(scout_path)
    issues: list[MonotoneScoutIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_exact_rows(scout))
    issues.extend(validate_countermodel(scout))
    issues.extend(validate_finite_zeta(scout))
    issues.extend(validate_summary(scout))
    issues.extend(validate_note(note_path))
    summary = scout.get("summary", {})
    return (
        issues,
        int(summary.get("exact_certificate_rows", 0)),
        int(summary.get("total_bernstein_coefficients", 0)),
        int(summary.get("finite_zeta_checked_rows", 0)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_SCOUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, exact_rows, bernstein_count, finite_rows = validate(args.scout, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "exact_rows": exact_rows,
                    "bernstein_coefficients": bernstein_count,
                    "finite_zeta_rows": finite_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-MONOTONE-CONTRACTION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF monotone contraction frontier scout: "
            f"{exact_rows} exact rows, {bernstein_count} Bernstein coefficients, "
            f"{finite_rows} finite zeta rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
