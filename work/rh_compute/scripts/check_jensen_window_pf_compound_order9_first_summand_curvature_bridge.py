#!/usr/bin/env python3
"""Validate the order-nine first/full-kernel curvature bridge."""

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

from jensen_window_pf_compound_order9_first_summand_curvature_bridge import (  # noqa: E402
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
    if artifact.get("kind") != "jensen_window_pf_compound_order9_first_summand_curvature_bridge":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "two-index finite splice" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 10,
        "ready_rows": 8,
        "open_rows": 2,
        "sixth_gap_floor_theorems": 1,
        "full_kernel_transfer_theorems": 1,
        "conditional_rows": 3,
        "open_continuous_targets": 1,
        "open_finite_splices": 1,
        "finite_splice_indices": 2,
        "transfer_polynomial_degree": 168,
        "positive_transfer_coefficients": 169,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    exact = artifact.get("exact", {})
    required = {
        "sixth_gap_floor": "min(V_j,V_j^(1))>=4/(3*j), j>=1250",
        "full_transfer": "|Y_k-Y_k^(1)|<=F_(k-1)+2*F_k+F_(k+1)<550/k^2, k>=1251",
        "continuous_target": "w_1''(t)<=4200/t^2 for every real t>=1250",
        "conditional_endpoint_tail": (
            "w_1''(t)<=4200/t^2 on t>=1250 => "
            "Q_(9,n)(-100)>0 for every n>=1243"
        ),
        "finite_splice": "prove Q_(9,n)(-100)>0 for n=1241,1242",
    }
    for key, expected in required.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    if exact.get("first_V_floor_polynomial", {}).get("shifted_coefficients") != [
        "13",
        "17490",
        "1542497",
    ]:
        issues.append(issue("exact", "bad-first-floor", exact.get("first_V_floor_polynomial")))
    if exact.get("full_V_floor_polynomial", {}).get("shifted_coefficients") != [
        "3271",
        "4140000",
        "62044250",
    ]:
        issues.append(issue("exact", "bad-full-floor", exact.get("full_V_floor_polynomial")))
    transfer = exact.get("transfer_polynomial", {})
    if transfer.get("degree") != 168 or transfer.get("coefficient_count") != 169:
        issues.append(issue("exact", "bad-transfer-polynomial", transfer))
    if transfer.get("minimum_coefficient") != "103950":
        issues.append(issue("exact", "bad-minimum-coefficient", transfer.get("minimum_coefficient")))
    try:
        endpoint = Decimal(exact.get("endpoint_scaled_transfer_decimal", "nan"))
        reserve = Decimal(exact.get("endpoint_scaled_reserve_below_550", "nan"))
    except Exception as exc:
        issues.append(issue("exact", "bad-endpoint-decimal", exc))
    else:
        if not (Decimal("533") < endpoint < Decimal("534")):
            issues.append(issue("exact", "bad-endpoint-transfer", endpoint))
        if not (Decimal("16") < reserve < Decimal("17")):
            issues.append(issue("exact", "bad-endpoint-reserve", reserve))
    rows = artifact.get("rows", [])
    if len(rows) != 10 or len({row.get("id") for row in rows}) != 10:
        issues.append(issue("rows", "bad-rows", rows))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 8:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 2:
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
        "This is not a",
        "min(V_j,V_j^(1))>=4/(3*j), j>=1250",
        "13*m**2 + 17490*m + 1542497>0",
        "3271*m**2 + 4140000*m + 62044250>0",
        "<550/k^2, k>=1251",
        "degree-168 shifted numerator has 169 strictly positive",
        "w_1''(t)<=4200/t^2 for every real t>=1250",
        "4751/k^2<4900/k^2",
        "prove Q_(9,n)(-100)>0 for n=1241,1242",
        "outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md",
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
                f"ORDER9-FIRST-FULL {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine first/full curvature bridge: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-nine first/full curvature bridge: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['sixth_gap_floor_theorems']} sixth-gap floor theorem, "
        f"{summary['positive_transfer_coefficients']} positive transfer coefficients, "
        f"{summary['full_kernel_transfer_theorems']} full-kernel transfer theorem, "
        f"{summary['open_continuous_targets']} open continuous target, "
        f"{summary['finite_splice_indices']} finite splice indices"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
