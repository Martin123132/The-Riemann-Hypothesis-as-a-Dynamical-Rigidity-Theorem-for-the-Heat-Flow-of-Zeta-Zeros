#!/usr/bin/env python3
"""Validate the Gaussian-cumulant target for the remaining order-four ray."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order4_gaussian_cumulant_ray_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


EXPECTED_SUMMARY = (
    "validated order-four Gaussian cumulant ray target: 10 rows, 0 issues, "
    "4 exact formal rows, 2 proved formal-corridor rows, "
    "7 positive conditional collars, 3 open analytic rows"
)


@dataclass(frozen=True)
class Finding:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> Finding:
    return Finding(section=section, issue=issue, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    try:
        artifact = load_json(artifact_path)
    except Exception as exc:
        return [finding("artifact", "load-failed", repr(exc))]
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))

    expected_summary = {
        "rows": 10,
        "exact_rows": 4,
        "formal_theorem_rows": 2,
        "open_analytic_rows": 3,
        "conditional_scout_rows": 1,
        "positive_conditional_collars": 7,
        "ready_to_apply_rows": 6,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    row_ids = {row.get("id") for row in artifact.get("rows", [])}
    expected_ids = {
        "co4gcrt_01_epsilon_grading",
        "co4gcrt_02_partition_recurrence",
        "co4gcrt_03_cumulant_polynomials",
        "co4gcrt_04_factorial_signature",
        "co4gcrt_05_candidate_corridor",
        "co4gcrt_06_conditional_collar_scout",
        "co4gcrt_06a_formal_finite_corridor",
        "co4gcrt_06b_formal_asymptotic_corridor",
        "co4gcrt_07_central_remainder_target",
        "co4gcrt_08_tail_and_composition_target",
    }
    if row_ids != expected_ids:
        findings.append(finding("rows", "bad-ids", sorted(str(value) for value in row_ids)))
    scouts = artifact.get("conditional_collar_scout", [])
    if len(scouts) != 7 or not all(row.get("passed_conditionally") for row in scouts):
        findings.append(finding("scouts", "failed", scouts))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    generator = artifact.get("generator")
    if not isinstance(generator, str) or not (REPO_ROOT / generator).exists():
        findings.append(finding("artifact", "missing-generator", generator))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact formal cumulant algebra, global formal-corridor theorem, and",
        "kappa_8~720*epsilon^6",
        "2/5<=q*(kappa_2-1)<=4/5",
        "These are theorem targets, not inferred from the formal series.",
        "conditional finite compatibility check, not a continuum ray",
        "for every `u>=2`",
        "not a proof of the exact-density",
    )
    for marker in required:
        if marker not in text:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-four Gaussian cumulant ray target: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
