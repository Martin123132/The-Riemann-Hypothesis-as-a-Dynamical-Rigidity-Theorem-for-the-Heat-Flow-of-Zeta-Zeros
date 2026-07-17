#!/usr/bin/env python3
"""Validate the order-seven first/full curvature bridge."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order7_first_summand_curvature_bridge import (  # noqa: E402
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
        "jensen_window_pf_compound_order7_first_summand_curvature_bridge"
    ):
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    expected = {
        "rows": 9,
        "ready_rows": 8,
        "open_rows": 1,
        "fourth_gap_floor_theorems": 1,
        "full_kernel_transfer_theorems": 1,
        "conditional_rows": 3,
        "open_continuous_targets": 1,
        "transfer_polynomial_degree": 102,
        "positive_transfer_coefficients": 103,
    }
    for key, value in expected.items():
        actual = artifact.get("summary", {}).get(key)
        if actual != value:
            issues.append(Issue("summary", key, str(actual)))
    exact = artifact.get("exact", {})
    if len(exact.get("finite_T_floor_rows", [])) != 2:
        issues.append(
            Issue("exact", "finite-floor-count", str(exact.get("finite_T_floor_rows")))
        )
    if any(
        float(row.get("first_margin_decimal", "-1")) <= 0
        for row in exact.get("finite_T_floor_rows", [])
    ):
        issues.append(Issue("exact", "finite-floor-sign", "nonpositive"))
    transfer = exact.get("transfer_polynomial", {})
    if len(transfer.get("shifted_coefficients", [])) != 103:
        issues.append(Issue("exact", "transfer-coefficients", str(transfer)))
    if not exact.get("full_transfer", "").endswith("k>=321"):
        issues.append(Issue("exact", "full-transfer", str(exact.get("full_transfer"))))
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
            "Status: exact first/full-kernel transfer",
            "min(T_j,T_j^(1))>=3/(2*j), j>=320",
            "2/k^8",
            "<262/k^2",
            "degree-102 shifted numerator",
            "r_1''(t)<=600/t^2",
            "R_k<863/k^2<900/k^2",
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
        print(f"ORDER7-FIRST-FULL {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-seven first/full curvature bridge: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-seven first/full curvature bridge: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['fourth_gap_floor_theorems']} fourth-gap floor theorem, "
        f"{summary['full_kernel_transfer_theorems']} full transfer, "
        f"{summary['conditional_rows']} conditional theorems, "
        f"{summary['open_continuous_targets']} open continuous target, "
        f"degree {summary['transfer_polynomial_degree']}, "
        f"{summary['positive_transfer_coefficients']} positive coefficients"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
