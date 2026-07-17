#!/usr/bin/env python3
"""Validate the order-ten endpoint-tail curvature reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_m100_tail_curvature_reduction as target  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []

    if not args.artifact.exists():
        issues.append(f"missing artifact: {args.artifact}")
    else:
        artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = target.build_artifact()
        if artifact != expected:
            issues.append("artifact does not match a fresh deterministic build")
        summary = artifact.get("summary", {})
        expected_summary = {
            "rows": 10,
            "ready_to_apply_rows": 9,
            "open_rows": 1,
            "exact_factorizations": 1,
            "exact_curvature_reductions": 1,
            "coefficient_positive_comparisons": 1,
            "conditional_tail_theorems": 1,
            "conditional_finite_splices": 1,
            "conditional_heat_handoffs": 1,
            "open_curvature_targets": 1,
            "preserved_negative_endpoint_shifts": 4,
        }
        if summary != expected_summary:
            issues.append(f"bad summary: {summary!r}")
        exact = artifact.get("exact", {})
        if exact.get("shifted_coefficients") != ["2259", "2906536", "96616536"]:
            issues.append("tail comparison coefficients changed")
        if exact.get("sufficient_ceiling") != "Z_k<=5500/k^2 for every integer k>=1252":
            issues.append("curvature target changed")
        if "n>=4" not in str(exact.get("conditional_heat_handoff", "")):
            issues.append("heat handoff lost the n>=4 boundary")
        rows = artifact.get("rows", [])
        open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
        if len(open_rows) != 1 or open_rows[0].get("id") != "co10m100tcr_10_open_curvature":
            issues.append("unexpected open-row contract")

    if not args.note.exists():
        issues.append(f"missing note: {args.note}")
    else:
        note = args.note.read_text(encoding="utf-8")
        for required in (
            "5500/k^2",
            "`4<=n<=1242`",
            "`n=0,1,2,3` remain negative",
            "not a proof of RH",
        ):
            if required not in note:
                issues.append(f"note missing required text: {required!r}")

    if issues:
        print(f"order-ten endpoint-tail curvature reduction: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-ten endpoint-tail curvature reduction: "
        "10 rows, 0 issues, 1 exact factorization, 1 positive comparison, "
        "1 conditional n>=4 heat handoff, 1 open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
