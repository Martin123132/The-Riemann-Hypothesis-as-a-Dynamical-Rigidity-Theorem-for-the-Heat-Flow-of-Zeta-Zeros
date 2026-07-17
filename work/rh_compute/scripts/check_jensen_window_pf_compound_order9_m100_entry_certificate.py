#!/usr/bin/env python3
"""Validate all-shift signed order-nine entry at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order9_m100_entry_certificate as entry9  # noqa: E402


@dataclass(frozen=True)
class EntryIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> EntryIssue:
    return EntryIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[EntryIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[EntryIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order9_m100_entry_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "rigorous all-shift signed order-nine entry at lambda=-100"
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    exact = artifact.get("exact", {})
    expected_exact = {
        "continuous_first_summand": entry9.CONTINUOUS_THEOREM,
        "first_discrete_ceiling": entry9.FIRST_DISCRETE_THEOREM,
        "full_kernel_transfer": entry9.FULL_TRANSFER_THEOREM,
        "full_curvature_ceiling": entry9.FULL_CEILING_THEOREM,
        "analytic_tail": entry9.TAIL_THEOREM,
        "finite_prefix": entry9.FINITE_THEOREM,
        "all_shift_entry": entry9.ENTRY_THEOREM,
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    expected_summary = {
        "rows": 6,
        "ready_rows": 6,
        "open_rows": 0,
        "continuous_curvature_inputs": 1,
        "full_kernel_transfer_applications": 1,
        "analytic_tail_theorems": 1,
        "finite_prefix_theorems": 1,
        "all_shift_m100_entry_theorems": 1,
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
            source_path = entry9.REPO_ROOT / source["path"]
        except Exception as exc:
            issues.append(issue("sources", "bad-path", exc))
            continue
        if not source_path.exists():
            issues.append(issue("sources", "missing-source", source_path))
        elif entry9.sha256(source_path) != source.get("sha256"):
            issues.append(issue("sources", "hash-mismatch", source_path))
    try:
        canonical = entry9.build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[EntryIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous all-shift signed order-nine entry",
        "This is not a proof",
        entry9.CONTINUOUS_THEOREM,
        entry9.FIRST_DISCRETE_THEOREM,
        entry9.FULL_TRANSFER_THEOREM,
        entry9.FULL_CEILING_THEOREM,
        entry9.TAIL_THEOREM,
        entry9.FINITE_THEOREM,
        entry9.ENTRY_THEOREM,
        "-100<=lambda<=0",
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
    parser.add_argument("--artifact", type=Path, default=entry9.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=entry9.DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-M100-ENTRY {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine lambda=-100 entry: {len(issues)} issues")
        return 1
    print(
        "validated all-shift signed order-nine entry at lambda=-100: "
        "0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
