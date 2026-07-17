#!/usr/bin/env python3
"""Validate the order-eight high-cumulant coarse corridor."""

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

from jensen_window_pf_compound_order8_high_cumulant_coarse_corridor import (  # noqa: E402
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
        "jensen_window_pf_compound_order8_high_cumulant_coarse_corridor"
    ):
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    summary = artifact.get("summary", {})
    exact = artifact.get("exact", {})
    for key, value in {
        "rows": 4,
        "ready_rows": 4,
        "formal_orders": 2,
        "formal_terms": 0,
        "cauchy_extensions": 1,
        "global_coarse_corridors": 2,
    }.items():
        if summary.get(key) != value:
            issues.append(Issue("summary", key, str(summary.get(key))))
    if exact.get("formal_orders") != [13, 14]:
        issues.append(Issue("exact", "orders", str(exact.get("formal_orders"))))
    if exact.get("cauchy_factor") != 182:
        issues.append(Issue("exact", "cauchy-factor", str(exact.get("cauchy_factor"))))
    if exact.get("finite_scaled_residual_bound") != str(Fraction(1729, 61875)):
        issues.append(
            Issue(
                "exact",
                "finite-residual",
                str(exact.get("finite_scaled_residual_bound")),
            )
        )
    if exact.get("ray_scaled_residual_bound") != f"{Fraction(100009, 4950000)}/u":
        issues.append(
            Issue(
                "exact",
                "ray-residual",
                str(exact.get("ray_scaled_residual_bound")),
            )
        )
    for order in ("13", "14"):
        if exact.get("formal_rows", {}).get(order) != {
            "scaled_expression": "0",
            "term_count": 0,
        }:
            issues.append(Issue("exact", f"formal-{order}", str(exact.get("formal_rows"))))
    if exact.get("exact_corridor_cap") != 1:
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
            "Status: global coarse exact thirteenth- and fourteenth-cumulant",
            "scaled kappa_13^[10]=scaled kappa_14^[10]=0",
            "14*13=182",
            "<1, r=13,14, u>=2",
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
        print(f"ORDER8-HIGH-CUMULANT {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-eight high-cumulant corridor: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight high-cumulant corridor: "
        f"{summary['formal_terms']} formal terms, 0 issues, "
        f"{summary['global_coarse_corridors']} exact corridors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
