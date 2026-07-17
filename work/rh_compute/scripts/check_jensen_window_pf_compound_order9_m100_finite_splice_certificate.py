#!/usr/bin/env python3
"""Validate the lambda=-100 order-nine finite splice certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order9_m100_finite_splice_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class SpliceIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> SpliceIssue:
    return SpliceIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[SpliceIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[SpliceIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order9_m100_finite_splice_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "two-index order-nine finite splice" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 5,
        "ready_rows": 5,
        "open_rows": 0,
        "coefficients": 1259,
        "new_coefficient_rows": 2,
        "positive_Q8_rows": 1245,
        "combined_positive_Q9_rows": 1243,
        "finite_splice_rows": 2,
        "finite_splice_theorems": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    exact = artifact.get("exact", {})
    required = {
        "splice": "Q_(9,n)(-100)>0 for n=1241,1242",
        "combined_prefix": "Q_(9,n)(-100)>0 for every 0<=n<=1242",
        "analytic_tail": "Q_(9,n)(-100)>0 for every n>=1243",
    }
    for key, expected in required.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    finite = artifact.get("finite", {})
    if finite.get("n_range") != [0, 1242] or finite.get("coefficient_range") != [0, 1258]:
        issues.append(issue("finite", "bad-ranges", finite))
    splice_rows = finite.get("splice_rows", [])
    if [row.get("n") for row in splice_rows] != [1241, 1242]:
        issues.append(issue("finite", "bad-splice-rows", splice_rows))
    if any(row.get("Q9_sign") != "positive_by_signed_condensation" for row in splice_rows):
        issues.append(issue("finite", "bad-splice-sign", splice_rows))
    try:
        minimum = Decimal(finite.get("minimum_relative_lower", "nan"))
    except Exception as exc:
        issues.append(issue("finite", "bad-minimum", exc))
    else:
        if finite.get("minimum_relative_n") != 1242:
            issues.append(issue("finite", "bad-minimum-index", finite.get("minimum_relative_n")))
        if not (Decimal("0.004") < minimum < Decimal("0.0041")):
            issues.append(issue("finite", "bad-minimum", minimum))
    rows = artifact.get("rows", [])
    if len(rows) != 5 or len({row.get("id") for row in rows}) != 5:
        issues.append(issue("rows", "bad-rows", rows))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "bad-readiness", rows))
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[SpliceIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous two-index endpoint splice",
        "This is not a proof of",
        "A_1257",
        "A_1258",
        "Q_(9,n)(-100)>0 for n=1241,1242",
        "Q_(9,n)(-100)>0 for every 0<=n<=1242",
        "meets the analytic bridge with no index gap",
        "w_1''(t)<=4200/t^2",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-nine entry is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-FINITE-SPLICE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine finite splice: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-nine finite splice: "
        f"{summary['coefficients']} coefficients, 0 issues, "
        f"{summary['finite_splice_rows']} splice rows, "
        f"{summary['combined_positive_Q9_rows']} combined positive signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
