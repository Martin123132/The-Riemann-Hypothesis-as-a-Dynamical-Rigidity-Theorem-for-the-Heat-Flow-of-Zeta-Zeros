#!/usr/bin/env python3
"""Validate the order-six first/full curvature bridge."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order6_first_summand_curvature_bridge import (  # noqa: E402
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
    if artifact.get("kind") != "jensen_window_pf_compound_order6_first_summand_curvature_bridge":
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    expected = {
        "rows": 8,
        "ready_rows": 7,
        "open_rows": 1,
        "third_gap_floor_theorems": 1,
        "full_kernel_transfer_theorems": 1,
        "conditional_rows": 2,
        "open_continuous_targets": 1,
    }
    for key, value in expected.items():
        if artifact.get("summary", {}).get(key) != value:
            issues.append(Issue("summary", key, str(artifact.get("summary", {}).get(key))))
    exact = artifact.get("exact", {})
    if exact.get("stable_derivative_residuals") != ["0", "0"]:
        issues.append(Issue("exact", "derivative-residuals", str(exact.get("stable_derivative_residuals"))))
    transfer = exact.get("transfer_polynomial", {})
    if transfer.get("degree") != 75 or transfer.get("coefficient_count") != 76:
        issues.append(Issue("exact", "transfer-polynomial", str(transfer)))
    if float(exact.get("endpoint_scaled_transfer_decimal", "100")) >= 100:
        issues.append(Issue("exact", "endpoint-transfer", str(exact.get("endpoint_scaled_transfer_decimal"))))
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
            "min(S_j,S_j^(1))>=3/(2*j)",
            "|P_k-P_k^(1)|",
            "p_1''(t)<=200/t^2",
            "degree-75 shifted numerator",
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
        print(f"ORDER6-FIRST-BRIDGE {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-six first/full curvature bridge: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-six first/full curvature bridge: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['full_kernel_transfer_theorems']} full transfer, "
        f"{summary['open_continuous_targets']} open continuous target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
