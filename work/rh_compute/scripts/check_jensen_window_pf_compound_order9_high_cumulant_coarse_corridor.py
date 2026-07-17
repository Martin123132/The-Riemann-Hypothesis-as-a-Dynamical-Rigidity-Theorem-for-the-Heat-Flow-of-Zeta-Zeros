#!/usr/bin/env python3
"""Validate the order-nine high-cumulant coarse corridor."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order9_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class CorridorIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CorridorIssue:
    return CorridorIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[CorridorIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[CorridorIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order9_high_cumulant_coarse_corridor":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 4,
        "ready_rows": 4,
        "formal_orders": 2,
        "formal_terms": 0,
        "cauchy_extensions": 1,
        "global_coarse_corridors": 2,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    exact = artifact.get("exact", {})
    if exact.get("formal_orders") != [15, 16]:
        issues.append(issue("exact", "bad-formal-orders", exact.get("formal_orders")))
    if exact.get("cauchy_factor") != 240:
        issues.append(issue("exact", "bad-cauchy-factor", exact.get("cauchy_factor")))
    if Fraction(exact.get("finite_scaled_residual_bound", "1")) != Fraction(152, 4125):
        issues.append(issue("exact", "bad-finite-residual", exact.get("finite_scaled_residual_bound")))
    if exact.get("ray_scaled_residual_bound") != "1099/41250/u":
        issues.append(issue("exact", "bad-ray-residual", exact.get("ray_scaled_residual_bound")))
    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[CorridorIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: global coarse exact fifteenth- and sixteenth-cumulant corridor",
        "This is not a proof",
        "scaled kappa_15^[10]=scaled kappa_16^[10]=0",
        "16*15=240",
        "finite residual < 152/4125",
        "ray residual < 1099/41250/u",
        "r=15,16, u>=2",
        "outputs/formal_core.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
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
                f"ORDER9-HIGH-CUMULANT {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine high-cumulant corridor: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-nine high-cumulant corridor: "
        f"{summary['formal_terms']} formal terms, 0 issues, "
        f"{summary['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
