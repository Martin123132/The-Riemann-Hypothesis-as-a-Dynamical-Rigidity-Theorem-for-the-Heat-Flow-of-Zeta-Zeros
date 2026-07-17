#!/usr/bin/env python3
"""Validate the order-eight first/full-kernel curvature bridge."""

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

from jensen_window_pf_compound_order8_first_summand_curvature_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class BridgeIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> BridgeIssue:
    return BridgeIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[BridgeIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[BridgeIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order8_first_summand_curvature_bridge":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "one open continuous theorem" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 9,
        "ready_rows": 8,
        "open_rows": 1,
        "fifth_gap_floor_theorems": 1,
        "full_kernel_transfer_theorems": 1,
        "conditional_rows": 3,
        "open_continuous_targets": 1,
        "transfer_polynomial_degree": 133,
        "positive_transfer_coefficients": 134,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    exact = artifact.get("exact", {})
    required = {
        "fifth_gap_floor": "min(U_j,U_j^(1))>=3/(2*j), j>=1249",
        "full_transfer": "|W_k-W_k^(1)|<=C_(k-1)+2*C_k+C_(k+1)<190/k^2, k>=1250",
        "continuous_target": "s_1''(t)<=4000/t^2 for every real t>=999",
        "conditional_endpoint_tail": (
            "s_1''(t)<=4000/t^2 on t>=999 => "
            "Q_(8,n)(-100)>0 for every n>=1243"
        ),
    }
    for key, expected in required.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    if exact.get("first_U_floor_polynomial", {}).get("shifted_coefficients") != [
        "6", "12581", "6352461"
    ]:
        issues.append(issue("exact", "bad-first-floor", exact.get("first_U_floor_polynomial")))
    if exact.get("full_U_floor_polynomial", {}).get("shifted_coefficients") != [
        "756", "1456613", "639733131"
    ]:
        issues.append(issue("exact", "bad-full-floor", exact.get("full_U_floor_polynomial")))
    transfer = exact.get("transfer_polynomial", {})
    if transfer.get("degree") != 133 or transfer.get("coefficient_count") != 134:
        issues.append(issue("exact", "bad-transfer-polynomial", transfer))
    try:
        endpoint = Decimal(exact.get("endpoint_scaled_transfer_decimal", "nan"))
        reserve = Decimal(exact.get("endpoint_scaled_reserve_below_190", "nan"))
    except Exception as exc:
        issues.append(issue("exact", "bad-endpoint-decimal", exc))
    else:
        if not (Decimal("177") < endpoint < Decimal("179")):
            issues.append(issue("exact", "bad-endpoint-transfer", endpoint))
        if not (Decimal("12") < reserve < Decimal("13")):
            issues.append(issue("exact", "bad-endpoint-reserve", reserve))
    rows = artifact.get("rows", [])
    if len(rows) != 9 or len({row.get("id") for row in rows}) != 9:
        issues.append(issue("rows", "bad-rows", rows))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 8:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        issues.append(issue("rows", "bad-open-count", rows))
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[BridgeIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact first/full-kernel transfer with one open continuous",
        "This is not a proof of order-eight entry",
        "min(U_j,U_j^(1))>=3/(2*j), j>=1249",
        "6*m**2 + 12581*m + 6352461>0",
        "756*m**2 + 1456613*m + 639733131>0",
        "<190/k^2, k>=1250",
        "degree-133 shifted numerator has 134 strictly positive",
        "s_1''(t)<=4000/t^2 for every real t>=999",
        "4191/k^2<4300/k^2",
        "outputs/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-eight entry is proved",
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
                f"ORDER8-FIRST-FULL {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-eight first/full curvature bridge: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight first/full curvature bridge: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['fifth_gap_floor_theorems']} fifth-gap floor theorem, "
        f"{summary['positive_transfer_coefficients']} positive transfer coefficients, "
        f"{summary['full_kernel_transfer_theorems']} full-kernel transfer theorem, "
        f"{summary['open_continuous_targets']} open continuous target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
