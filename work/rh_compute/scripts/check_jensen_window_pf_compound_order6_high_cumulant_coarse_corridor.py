#!/usr/bin/env python3
"""Validate the order-six coarse ninth/tenth cumulant corridor."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order6_high_cumulant_coarse_corridor import (  # noqa: E402
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
    if artifact.get("kind") != "jensen_window_pf_compound_order6_high_cumulant_coarse_corridor":
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    expected = {
        "rows": 5,
        "ready_rows": 5,
        "formal_orders": 2,
        "formal_terms": 109,
        "cauchy_extensions": 1,
        "global_coarse_corridors": 2,
    }
    for key, value in expected.items():
        if artifact.get("summary", {}).get(key) != value:
            issues.append(Issue("summary", key, str(artifact.get("summary", {}).get(key))))
    exact = artifact.get("exact", {})
    if exact.get("exact_corridor_cap") != 50000:
        issues.append(Issue("exact", "cap", str(exact.get("exact_corridor_cap"))))
    for order in ("9", "10"):
        row = exact.get("formal_rows", {}).get(order, {})
        if not row or row.get("finite_cap") != 1600 or row.get("ray_cap") != 36000:
            issues.append(Issue("exact", f"formal-{order}", str(row)))
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
            "Status: global coarse exact ninth- and tenth-cumulant",
            "r*(r-1)<=90",
            "<50000",
            "No sign is asserted",
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
        print(f"ORDER6-HIGH-CUMULANT {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-six high-cumulant corridor: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-six high-cumulant corridor: "
        f"{summary['formal_orders']} formal orders, "
        f"{summary['formal_terms']} terms, 0 issues, "
        f"{summary['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
