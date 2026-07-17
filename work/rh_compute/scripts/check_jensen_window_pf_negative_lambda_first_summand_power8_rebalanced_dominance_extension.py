#!/usr/bin/env python3
"""Validate the rebalanced inverse-eighth-power dominance extension."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class Issue:
    section: str
    code: str
    detail: str


def validate(artifact_path: Path, note_path: Path) -> tuple[dict, list[Issue]]:
    if not artifact_path.exists():
        return {}, [Issue("artifact", "missing", str(artifact_path))]
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    issues: list[Issue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_negative_lambda_first_summand_"
        "power8_rebalanced_dominance_extension"
    ):
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    expected = {
        "rows": 7,
        "ready_to_apply_rows": 7,
        "positive_analytic_gates": 14,
        "full_tail_power": 8,
        "tail_start_k": 300,
        "dominance_theorems": 1,
        "rebalanced_splits": 1,
    }
    for key, value in expected.items():
        actual = artifact.get("summary", {}).get(key)
        if actual != value:
            issues.append(Issue("summary", key, str(actual)))
    diagnostics = artifact.get("diagnostics", {})
    if not diagnostics.get("all_positive_gates"):
        issues.append(
            Issue("diagnostics", "gates", str(diagnostics.get("positive_gates")))
        )
    if "2/k^8" not in diagnostics.get("full_tail_relative_bound", ""):
        issues.append(
            Issue(
                "diagnostics",
                "bound",
                str(diagnostics.get("full_tail_relative_bound")),
            )
        )
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(Issue("rebuild", "exception", repr(exc)))
    else:
        if canonical != artifact:
            issues.append(Issue("rebuild", "drift", str(artifact_path)))
    if not note_path.exists():
        issues.append(Issue("note", "missing", str(note_path)))
    else:
        text = note_path.read_text(encoding="utf-8")
        for marker in (
            "Status: rigorous inverse-eighth-power",
            "a(k)=log(k)/10",
            "2/k^8",
            "k>=300",
            "All fourteen endpoint and derivative gates",
            "This is not a proof",
            "outputs/formal_core.md",
        ):
            if marker not in text:
                issues.append(Issue("note", "marker", marker))
    return artifact, issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate(args.artifact, args.note)
    for item in issues:
        print(f"POWER8-DOMINANCE {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"power-eight rebalanced first-summand dominance: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated power-eight rebalanced first-summand dominance: "
        f"{summary['rows']} rows, {summary['positive_analytic_gates']} positive gates, "
        f"tail k>={summary['tail_start_k']}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
