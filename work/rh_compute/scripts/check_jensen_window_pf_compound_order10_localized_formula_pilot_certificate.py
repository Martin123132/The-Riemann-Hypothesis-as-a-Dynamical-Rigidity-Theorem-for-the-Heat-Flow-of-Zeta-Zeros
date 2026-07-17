#!/usr/bin/env python3
"""Validate the order-ten localized formula pilot certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_localized_formula_pilot_certificate as target  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    issues: list[str] = []

    if not args.artifact.exists():
        issues.append(f"missing artifact: {args.artifact}")
        artifact = {}
    else:
        artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    if artifact:
        if artifact.get("kind") != (
            "jensen_window_pf_compound_order10_localized_formula_pilot_certificate"
        ):
            issues.append(f"bad kind: {artifact.get('kind')!r}")
        if artifact.get("theorem") != (
            "z_1''(t)<=5500/t^2 for every real 1251<=t<=1251.5"
        ):
            issues.append("pilot theorem changed")
        expected_summary = {
            "rows": 4,
            "ready_rows": 4,
            "open_rows": 0,
            "accepted_blocks": 2,
            "first_summand_local_curvature_theorems": 1,
            "full_half_line_theorems": 0,
            "full_kernel_theorems": 0,
        }
        if artifact.get("summary") != expected_summary:
            issues.append(f"bad summary: {artifact.get('summary')!r}")
        parameters = artifact.get("parameters", {})
        if parameters != {
            "start_t": "1251",
            "end_t": "2503/2",
            "block_width": "1/4",
            "point_order": 7,
            "remainder_order": 8,
            "curvature_constant": 5500,
            "precision_bits": 512,
        }:
            issues.append(f"bad parameters: {parameters!r}")
        source_contract = artifact.get("source_contract", {})
        if source_contract.get("H_rows") != 1700:
            issues.append("bad H source row count")
        if source_contract.get("point_rows") != 133:
            issues.append("bad point source row count")
        if "before parsing" not in str(
            source_contract.get("precision_load_invariant", "")
        ):
            issues.append("missing precision-load invariant")
        for record in source_contract.get("files", []):
            path = target.REPO_ROOT / record.get("path", "")
            if not path.exists() or target.sha256(path) != record.get("sha256"):
                issues.append(f"source hash mismatch: {path}")

        blocks = artifact.get("blocks", [])
        expected_left = Fraction(1251)
        if len(blocks) != 2:
            issues.append(f"bad block count: {len(blocks)}")
        for index, block in enumerate(blocks):
            left = Fraction(block.get("anchor", "0"))
            right = Fraction(block.get("right", "0"))
            if left != expected_left or right - left != Fraction(1, 4):
                issues.append(f"block {index} is not contiguous: {left}..{right}")
            expected_left = right
            if block.get("passed") is not True:
                issues.append(f"block {index} did not pass")
            expected_expansion = "1251" if index == 0 else "2503/2"
            if block.get("expansion_anchor") != expected_expansion:
                issues.append(
                    f"bad expansion anchor in block {index}: "
                    f"{block.get('expansion_anchor')!r}"
                )
            try:
                scaled = Decimal(block.get("scaled_curvature_upper", "nan"))
                margin = Decimal(block.get("curvature_margin_lower", "nan"))
                point = Decimal(block.get("point_scaled_curvature", "nan"))
            except Exception as exc:
                issues.append(f"bad decimal in block {index}: {exc}")
            else:
                if not (Decimal(0) < scaled < Decimal(5500)):
                    issues.append(f"bad scaled upper in block {index}: {scaled}")
                if margin <= Decimal(0):
                    issues.append(f"bad margin in block {index}: {margin}")
                if not (Decimal(29) < point < Decimal(30)):
                    issues.append(f"bad point curvature in block {index}: {point}")
            if "-chi(W)*(W')^2" not in str(block.get("proof_formula", "")):
                issues.append(f"block {index} lost the formula boundary")
        if expected_left != Fraction(2503, 2):
            issues.append(f"blocks end at {expected_left}, not 1251.5")

        rows = artifact.get("rows", [])
        if len(rows) != 4 or len({row.get("id") for row in rows}) != 4:
            issues.append("bad theorem-ledger rows")
        if any(row.get("readiness") != "ready_to_apply" for row in rows):
            issues.append("non-ready theorem-ledger row")
        try:
            target.validate_source_manifests()
        except Exception as exc:
            issues.append(f"source manifest validation failed: {exc}")
        if args.rebuild:
            try:
                canonical = target.build_artifact()
            except Exception as exc:
                issues.append(f"canonical rebuild failed: {exc}")
            else:
                if artifact != canonical:
                    issues.append("artifact differs from full canonical rebuild")

    if not args.note.exists():
        issues.append(f"missing note: {args.note}")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in (
            "rigorous first-summand local curvature theorem",
            "z_1''(t)<=5500/t^2 for every real 1251<=t<=1251.5",
            "two exact width-`1/4` blocks",
            "not",
            "full Newman kernel",
        ):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")

    if issues:
        print(f"order-ten localized formula pilot: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-ten localized formula pilot: "
        "2 contiguous blocks, 0 issues, 1 local first-summand theorem, "
        "0 full-half-line claims, 0 full-kernel claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
