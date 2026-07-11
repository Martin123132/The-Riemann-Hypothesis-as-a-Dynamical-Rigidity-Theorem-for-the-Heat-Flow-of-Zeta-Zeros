#!/usr/bin/env python3
"""Validate the Jensen-window PF contraction-log-concavity scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_contraction_log_concavity_scout.md"

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Contraction-Log-Concavity Scout",
    "Status: ratio-condition rejection diagnostic",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_contraction_log_concavity_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py",
    "validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows",
    "rho = 33/40",
    "x1 = 65/66",
    "x2 = 44/65",
    "x3 = 1/3",
    "x2^2 - x1*x3 = 108703/836550",
    "degree 3 size 8",
    "degree 4 size 6",
    "rc_07_contraction_log_concavity",
)


@dataclass(frozen=True)
class ContractionIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ContractionIssue:
    return ContractionIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_refs(scout: dict) -> list[ContractionIssue]:
    ref = scout.get("source_algebra")
    if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
        return [issue("<scout>", "missing-source-algebra", repr(ref))]
    return []


def validate_top_level(scout: dict) -> list[ContractionIssue]:
    issues: list[ContractionIssue] = []
    if scout.get("kind") != "jensen_window_pf_contraction_log_concavity_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    if scout.get("target_ratio_condition") != "rc_07_contraction_log_concavity":
        issues.append(issue("<scout>", "bad-target-ratio-condition", repr(scout.get("target_ratio_condition"))))
    boundary = str(scout.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<scout>", "weak-proof-boundary", str(scout.get("proof_boundary", ""))))

    summary = scout.get("summary", {})
    expected = {
        "candidate_rows_tested": 1,
        "rejected_by_constructed_extension_rows": 1,
        "negative_frontier_rows": 2,
        "target_closing_rows": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("<summary>", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if "still negative" not in str(summary.get("main_finding", "")).lower():
        issues.append(issue("<summary>", "missing-negative-finding", str(summary.get("main_finding", ""))))
    return issues


def validate_constructed_extension(scout: dict) -> list[ContractionIssue]:
    issues: list[ContractionIssue] = []
    extension = scout.get("constructed_extension", {})
    point = extension.get("ratio_point", {})
    expected_point = {"rho": "33/40", "x1": "65/66", "x2": "44/65", "x3": "1/3"}
    for key, value in expected_point.items():
        if point.get(key) != value:
            issues.append(issue("constructed_extension", f"bad-{key}", f"{point.get(key)!r} != {value!r}"))

    expected_sequence = ["1", "33/40", "429/640", "4719/12800", "17303/256000"]
    if extension.get("sequence_A0_to_A4") != expected_sequence:
        issues.append(issue("constructed_extension", "bad-sequence", repr(extension.get("sequence_A0_to_A4"))))

    values = extension.get("condition_values", {})
    expected_values = {
        "1 - x1": "1/66",
        "1 - x2": "21/65",
        "1 - x3": "2/3",
        "x1 - x2": "1321/4290",
        "x2 - x3": "67/195",
        "(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)": "3793/38025",
        "x2^2 - x1*x3": "108703/836550",
    }
    for key, value in expected_values.items():
        if values.get(key) != value:
            issues.append(issue("constructed_extension", f"bad-condition-{key}", repr(values.get(key))))

    satisfies = extension.get("satisfies", {})
    for key in (
        "adjacent_log_concavity_box",
        "decreasing_ratio_contractions",
        "second_order_log_concavity",
        "contraction_log_concavity",
    ):
        if satisfies.get(key) is not True:
            issues.append(issue("constructed_extension", f"condition-not-satisfied-{key}", repr(satisfies.get(key))))
    return issues


def validate_frontiers(scout: dict) -> list[ContractionIssue]:
    issues: list[ContractionIssue] = []
    frontiers = scout.get("frontier_failures", {})
    rows = frontiers.get("rows", [])
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    expected = {
        "d3_contiguous_m8": "-435846079534239/104857600000000",
        "d4_contiguous_m6": "-26359418151/4096000000",
    }
    for row_id, value in expected.items():
        row = by_id.get(row_id)
        if row is None:
            issues.append(issue("frontier_failures", "missing-row", row_id))
            continue
        if row.get("value") != value:
            issues.append(issue(row_id, "bad-value", f"{row.get('value')!r} != {value!r}"))
        if row.get("sign") != -1 or row.get("positive") is not False:
            issues.append(issue(row_id, "not-negative", repr(row)))
    if frontiers.get("negative_rows") != 2:
        issues.append(issue("frontier_failures", "bad-negative-row-count", repr(frontiers.get("negative_rows"))))

    comparison = scout.get("comparison_to_original_countermodel", {})
    if comparison.get("same_A0_to_A3_prefix") is not True:
        issues.append(issue("comparison_to_original_countermodel", "prefix-not-preserved", repr(comparison)))
    if comparison.get("d3_m8_matches_original_countermodel") is not True:
        issues.append(issue("comparison_to_original_countermodel", "d3-not-matched", repr(comparison)))
    return issues


def validate_note(path: Path) -> list[ContractionIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ContractionIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "kernel identity is proved", "cauchy-binet identity is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, note_path: Path) -> tuple[list[ContractionIssue], int, int]:
    scout = load_json(scout_path)
    issues: list[ContractionIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_refs(scout))
    issues.extend(validate_constructed_extension(scout))
    issues.extend(validate_frontiers(scout))
    issues.extend(validate_note(note_path))
    summary = scout.get("summary", {})
    return (
        issues,
        int(summary.get("rejected_by_constructed_extension_rows", 0)),
        int(summary.get("negative_frontier_rows", 0)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_SCOUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rejected_count, negative_count = validate(args.scout, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "rejected_by_constructed_extension_rows": rejected_count,
                    "negative_frontier_rows": negative_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-CONTRACTION-LOG-CONCAVITY {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF contraction-log-concavity scout: "
            f"{rejected_count} rejected by construction, {len(issues)} issues, "
            f"{negative_count} negative frontier rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
