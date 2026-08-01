#!/usr/bin/env python3
"""Validate the complete order-eleven compact adaptive-H23 certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_adaptive_h23_certificate as certificate  # noqa: E402
import jensen_window_pf_compound_order11_compact_adaptive_h23_segments as source  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=certificate.DEFAULT_CACHE)
    parser.add_argument(
        "--run-contract",
        type=Path,
        default=certificate.DEFAULT_RUN_CONTRACT,
    )
    parser.add_argument("--artifact", type=Path, default=certificate.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=certificate.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    try:
        actual = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = certificate.build_artifact(args.cache, args.run_contract)
    except (OSError, json.JSONDecodeError, RuntimeError, KeyError, ValueError) as exc:
        print(f"order-eleven compact certificate: source failure: {exc}")
        return 1
    if actual != expected:
        issues.append("artifact differs from deterministic cache/source rebuild")
    summary = actual.get("summary", {})
    if (
        actual.get("kind")
        != "jensen_window_pf_compound_order11_compact_adaptive_h23_certificate"
        or actual.get("status")
        != (
            "rigorous order-eleven first-summand curvature theorem on "
            "5700<=t<=38020"
        )
        or actual.get("theorem") != certificate.THEOREM
        or actual.get("generator") != certificate.GENERATOR_PATH
        or actual.get("checker") != certificate.CHECKER_PATH
    ):
        issues.append("artifact identity, status, or theorem changed")
    rows = actual.get("rows", [])
    if (
        len(rows) != 4
        or any(row.get("readiness") != "ready_to_apply" for row in rows)
        or rows[2].get("formula") != certificate.THEOREM
        or rows[3].get("formula") != "V'(2001/1000)<38020"
    ):
        issues.append("compact theorem ledger changed")
    if (
        summary.get("segments") != source.EXPECTED_SEGMENTS
        or summary.get("quarter_blocks") != source.EXPECTED_QUARTER_BLOCKS
        or summary.get("compact_first_summand_theorems") != 1
        or summary.get("global_first_summand_theorems") != 0
        or summary.get("full_kernel_theorems") != 0
        or summary.get("rh_claims") != 0
    ):
        issues.append("coverage or claim-boundary summary changed")
    try:
        largest = Decimal(summary["largest_scaled_curvature_upper"])
        smallest = Decimal(summary["smallest_margin_lower"])
        theorem_margin = Decimal(summary["theorem_margin_lower"])
        transition = Decimal(summary["saddle_transition_upper"])
        if not largest < Decimal(source.CURVATURE_CONSTANT):
            issues.append("largest scaled upper is not below 6000")
        if not smallest > 0 or not theorem_margin > 0:
            issues.append("compact curvature margin is not positive")
        if theorem_margin != Decimal(source.CURVATURE_CONSTANT) - largest:
            issues.append("theorem margin is not the exact global difference")
        compact_end = Decimal(source.END_T.numerator) / Decimal(
            source.END_T.denominator
        )
        if not transition < compact_end:
            issues.append("compact-saddle overlap is not strict")
    except (KeyError, InvalidOperation, TypeError, ValueError) as exc:
        issues.append(f"numeric summary parse failure: {exc}")
    cache_contract = actual.get("source_contract", {}).get("segment_cache", {})
    if (
        cache_contract.get("sha256") != summary.get("sha256")
        or cache_contract.get("bytes") != summary.get("bytes")
        or cache_contract.get("segments") != source.EXPECTED_SEGMENTS
        or cache_contract.get("quarter_blocks")
        != source.EXPECTED_QUARTER_BLOCKS
    ):
        issues.append("cache binding changed")
    if not args.note.exists():
        issues.append("missing note")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in ("2020", "129280", "V'(2001/1000)", "not a proof of RH"):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-eleven compact certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated complete order-eleven compact theorem: "
        "2020 segments, 129280 quarter blocks, 0 full-kernel claims, "
        "0 RH claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
