#!/usr/bin/env python3
"""Validate the global order-nine first-summand curvature certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order9_first_summand_curvature_certificate as global9  # noqa: E402


@dataclass(frozen=True)
class CurvatureIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CurvatureIssue:
    return CurvatureIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[CurvatureIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[CurvatureIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order9_first_summand_curvature_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("theorem") != global9.GLOBAL_THEOREM:
        issues.append(issue("artifact", "bad-theorem", artifact.get("theorem")))
    try:
        largest = Decimal(
            artifact.get("largest_scaled_curvature_upper", "Infinity")
        )
    except Exception as exc:
        issues.append(issue("artifact", "bad-largest", exc))
    else:
        if largest >= Decimal(4200):
            issues.append(issue("artifact", "curvature-failure", largest))
    expected_summary = {
        "rows": 4,
        "ready_rows": 4,
        "open_rows": 0,
        "lower_bridge_theorems": 1,
        "upper_range_compositions": 1,
        "global_first_summand_curvature_theorems": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    sources = artifact.get("source_contract", {}).get("sources", [])
    if len(sources) != 4:
        issues.append(issue("sources", "bad-count", len(sources)))
    for source in sources:
        try:
            source_path = global9.REPO_ROOT / source["path"]
        except Exception as exc:
            issues.append(issue("sources", "bad-path", exc))
            continue
        if not source_path.exists():
            issues.append(issue("sources", "missing-source", source_path))
        elif global9.sha256(source_path) != source.get("sha256"):
            issues.append(issue("sources", "hash-mismatch", source_path))
    try:
        canonical = global9.build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[CurvatureIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous global first-summand theorem on `t>=1250`",
        "This is not a proof",
        global9.LOWER_THEOREM,
        global9.COMPACT_THEOREM,
        global9.FINITE_RAY_THEOREM,
        global9.ASYMPTOTIC_THEOREM,
        global9.UPPER_COMPOSITION,
        global9.GLOBAL_THEOREM,
        "full-kernel",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=global9.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=global9.DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-GLOBAL-CURVATURE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"global order-nine first-summand curvature: {len(issues)} issues")
        return 1
    print(
        "validated global order-nine first-summand curvature: 0 issues, "
        f"largest scaled upper {artifact['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
