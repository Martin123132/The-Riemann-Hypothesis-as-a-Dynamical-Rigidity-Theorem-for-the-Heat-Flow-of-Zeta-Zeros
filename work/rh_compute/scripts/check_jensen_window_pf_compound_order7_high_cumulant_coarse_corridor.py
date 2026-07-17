#!/usr/bin/env python3
"""Validate the order-seven high-cumulant coarse corridor."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order7_high_cumulant_coarse_corridor import (  # noqa: E402
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
        "jensen_window_pf_compound_order7_high_cumulant_coarse_corridor"
    ):
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    summary = artifact.get("summary", {})
    exact = artifact.get("exact", {})
    for key, value in {
        "rows": 5,
        "ready_rows": 5,
        "formal_orders": 2,
        "formal_terms": 72,
        "cauchy_extensions": 1,
        "global_coarse_corridors": 2,
    }.items():
        if summary.get(key) != value:
            issues.append(Issue("summary", key, str(summary.get(key))))
    if exact.get("formal_orders") != [11, 12]:
        issues.append(Issue("exact", "orders", str(exact.get("formal_orders"))))
    if exact.get("cauchy_factor") != 132:
        issues.append(Issue("exact", "cauchy-factor", str(exact.get("cauchy_factor"))))
    if exact.get("finite_exact_corridor_cap") != 14001:
        issues.append(Issue("exact", "finite-cap", str(exact.get("finite_exact_corridor_cap"))))
    if exact.get("ray_exact_corridor_cap") != 700001:
        issues.append(Issue("exact", "ray-cap", str(exact.get("ray_exact_corridor_cap"))))
    if exact.get("exact_corridor_cap") != 1000000:
        issues.append(Issue("exact", "corridor-cap", str(exact.get("exact_corridor_cap"))))
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
            "Status: global coarse exact eleventh- and twelfth-cumulant",
            "L_3,...,L_12",
            "12*11=132",
            "1000000, r=11,12",
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
        print(f"ORDER7-HIGH-CUMULANT {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-seven high-cumulant corridor: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-seven high-cumulant corridor: "
        f"{summary['formal_terms']} formal terms, 0 issues, "
        f"{summary['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
