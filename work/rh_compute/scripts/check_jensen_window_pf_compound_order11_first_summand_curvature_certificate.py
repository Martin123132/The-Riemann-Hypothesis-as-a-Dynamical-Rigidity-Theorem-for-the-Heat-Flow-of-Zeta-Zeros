#!/usr/bin/env python3
"""Validate the global order-eleven first-summand curvature certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_first_summand_curvature_certificate as certificate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=certificate.DEFAULT_OUT)
    args = parser.parse_args()
    issues: list[str] = []
    if not args.artifact.exists():
        print(f"missing artifact: {args.artifact}")
        return 1
    try:
        actual = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = certificate.build_artifact()
    except (OSError, json.JSONDecodeError, RuntimeError, KeyError, ValueError) as exc:
        print(f"order-eleven first-summand certificate: source failure: {exc}")
        return 1
    if actual != expected:
        issues.append("artifact differs from deterministic live-source rebuild")
    if (
        actual.get("kind")
        != "jensen_window_pf_compound_order11_first_summand_curvature_certificate"
        or actual.get("status")
        != (
            "rigorous global order-eleven first-summand curvature theorem "
            "on t>=1252"
        )
        or actual.get("theorem") != certificate.GLOBAL_THEOREM
        or actual.get("generator") != certificate.GENERATOR_PATH
        or actual.get("checker") != certificate.CHECKER_PATH
    ):
        issues.append("artifact identity, status, or theorem changed")
    rows = actual.get("rows", [])
    if (
        len(rows) != 5
        or any(row.get("readiness") != "ready_to_apply" for row in rows)
        or rows[-1].get("formula") != certificate.GLOBAL_THEOREM
    ):
        issues.append("composition ledger changed")
    summary = actual.get("summary", {})
    if (
        summary.get("global_first_summand_curvature_theorems") != 1
        or summary.get("full_kernel_theorems") != 0
        or summary.get("heat_forward_theorems") != 0
        or summary.get("rh_claims") != 0
        or summary.get("open_rows") != 0
    ):
        issues.append("claim-boundary summary changed")
    try:
        largest = Decimal(actual["largest_scaled_curvature_upper"])
        overlap = actual["source_contract"]["compact_saddle_overlap"]
        transition = Decimal(overlap.split("<=", 1)[1].split("<", 1)[0])
        scaled = actual["source_contract"]["scaled_uppers"]
        if not largest < Decimal(certificate.GLOBAL_CURVATURE_CONSTANT):
            issues.append("largest scaled upper is not below 6000")
        if not transition < Decimal(38020):
            issues.append("compact-saddle overlap is not strict")
        if set(scaled) != {"lower", "compact", "finite_ray", "asymptotic_ray"}:
            issues.append("four-range scaled summary changed")
        elif largest != max(Decimal(value) for value in scaled.values()):
            issues.append("global largest scaled upper is inconsistent")
    except (KeyError, InvalidOperation, IndexError, TypeError, ValueError) as exc:
        issues.append(f"numeric summary parse failure: {exc}")
    sources = actual.get("source_contract", {}).get("sources", [])
    if len(sources) != 4 or len({row.get("sha256") for row in sources}) != 4:
        issues.append("four-source hash contract changed")
    if issues:
        print(f"order-eleven first-summand curvature certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated global order-eleven first-summand curvature certificate: "
        "1 half-line first-summand theorem, 0 full-kernel claims, 0 RH claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
