#!/usr/bin/env python3
"""Validate canonical sparse-H23 order-eleven lower-bridge segments."""

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

import jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments as source  # noqa: E402


def validate_run_contract(path: Path) -> list[str]:
    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
        canonical = source.canonical_run_contract()
    except (json.JSONDecodeError, OSError, RuntimeError, ValueError) as exc:
        return [f"run-contract validation failed: {exc}"]
    if existing != canonical:
        return ["run contract differs from canonical source and code hashes"]
    return []


def validate_block(
    block: dict,
    *,
    expected_left: Fraction,
    expected_expansion: Fraction,
    regime: str,
    label: str,
) -> list[str]:
    issues = []
    expected_right = expected_left + source.QUARTER_WIDTH
    try:
        if (
            Fraction(block["anchor"]) != expected_left
            or Fraction(block["right"]) != expected_right
            or Fraction(block["width"]) != source.QUARTER_WIDTH
            or Fraction(block["expansion_anchor"]) != expected_expansion
            or block.get("local_domain")
            != [
                str(expected_left - expected_expansion),
                str(expected_right - expected_expansion),
            ]
            or block.get("regime") != regime
            or block.get("model_degrees") != [16, 15, 14]
            or block.get("maximum_h_derivative_order") != 17
            or block.get("stable_taylor_surplus") != 4
            or block.get("passed") is not True
        ):
            issues.append(f"{label}: geometry or model contract changed")
        scaled = Decimal(block["scaled_curvature_upper"])
        margin = Decimal(block["curvature_margin_lower"])
        remainder = Decimal(block["final_uniform_remainder_upper"])
        if not scaled < Decimal(source.CURVATURE_CONSTANT):
            issues.append(f"{label}: scaled curvature upper reaches the wall")
        if not margin > 0:
            issues.append(f"{label}: curvature margin is not positive")
        if not remainder >= 0:
            issues.append(f"{label}: final remainder is negative")
        if block.get("quadrature", {}).get("shift_count") != 19:
            issues.append(f"{label}: shifted source count changed")
        stages = block.get("stage_diagnostics", [])
        if [stage.get("stage") for stage in stages] != list(range(2, 10)):
            issues.append(f"{label}: stable-stage diagnostics changed")
    except (KeyError, ValueError, TypeError) as exc:
        issues.append(f"{label}: parse failure: {exc}")
    return issues


def validate_segments(
    path: Path,
    tasks: list[tuple[int, str, Fraction, Fraction]],
) -> tuple[list[str], int, int, Decimal | None, Decimal | None]:
    issues = []
    segment_count = 0
    block_count = 0
    global_maximum = None
    global_minimum_margin = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if segment_count >= len(tasks):
                issues.append("segment cache has too many rows")
                break
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                issues.append(f"line {line_number}: invalid JSON: {exc}")
                break
            index, regime, segment_left, segment_right = tasks[segment_count]
            label = f"segment {index}"
            blocks = record.get("blocks", [])
            expected_blocks = int(
                (segment_right - segment_left) / source.QUARTER_WIDTH
            )
            if (
                record.get("kind")
                != "order11_sparse_h23_lower_bridge_segment"
                or record.get("index") != index
                or record.get("regime") != regime
                or record.get("segment_left") != str(segment_left)
                or record.get("segment_right") != str(segment_right)
                or record.get("block_count") != expected_blocks
                or len(blocks) != expected_blocks
                or record.get("passed") is not True
            ):
                issues.append(f"{label}: identity or block count changed")
                if len(issues) >= 20:
                    break
            segment_scaled = []
            segment_margins = []
            for block_index, block in enumerate(blocks):
                expected_left = (
                    segment_left + block_index * source.QUARTER_WIDTH
                )
                cell_left = (
                    segment_left
                    + (block_index // 2) * source.CELL_WIDTH
                )
                expected_expansion = (
                    cell_left if block_index % 2 == 0 else cell_left + source.CELL_WIDTH
                )
                issues.extend(
                    validate_block(
                        block,
                        expected_left=expected_left,
                        expected_expansion=expected_expansion,
                        regime=regime,
                        label=f"{label} block {block_index}",
                    )
                )
                try:
                    segment_scaled.append(Decimal(block["scaled_curvature_upper"]))
                    segment_margins.append(Decimal(block["curvature_margin_lower"]))
                except (KeyError, ValueError):
                    pass
                if len(issues) >= 20:
                    break
            if segment_scaled and (
                Decimal(record["largest_scaled_curvature_upper"])
                != max(segment_scaled)
            ):
                issues.append(f"{label}: largest scaled curvature summary changed")
            if segment_margins and (
                Decimal(record["smallest_margin_lower"])
                != min(segment_margins)
            ):
                issues.append(f"{label}: smallest margin summary changed")
            if segment_scaled:
                local_maximum = max(segment_scaled)
                global_maximum = (
                    local_maximum
                    if global_maximum is None
                    else max(global_maximum, local_maximum)
                )
            if segment_margins:
                local_minimum = min(segment_margins)
                global_minimum_margin = (
                    local_minimum
                    if global_minimum_margin is None
                    else min(global_minimum_margin, local_minimum)
                )
            block_count += len(blocks)
            segment_count += 1
            if len(issues) >= 20:
                issues.append("validation stopped after 20 issues")
                break
    return (
        issues,
        segment_count,
        block_count,
        global_maximum,
        global_minimum_margin,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_SEGMENT_CACHE)
    parser.add_argument(
        "--run-contract",
        type=Path,
        default=source.DEFAULT_RUN_CONTRACT,
    )
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    tasks = source.deterministic_segments()
    issues = validate_run_contract(args.run_contract)
    if not args.cache.exists():
        issues.append(f"missing segment cache: {args.cache}")
        segment_count = block_count = 0
        maximum = minimum_margin = None
    else:
        try:
            (
                segment_issues,
                segment_count,
                block_count,
                maximum,
                minimum_margin,
            ) = validate_segments(args.cache, tasks)
            issues.extend(segment_issues)
        except (OSError, RuntimeError, ValueError) as exc:
            issues.append(f"segment-cache validation failed: {exc}")
            segment_count = block_count = 0
            maximum = minimum_margin = None
    if args.require_complete and segment_count != len(tasks):
        issues.append(f"expected {len(tasks)} segments, found {segment_count}")

    if issues:
        print(f"order-eleven sparse H23 lower bridge: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    status = "complete" if segment_count == len(tasks) else "valid resumable prefix"
    print(
        "validated order-eleven sparse H23 lower bridge "
        f"({status}): {segment_count}/{len(tasks)} segments, "
        f"{block_count} quarter blocks, maximum scaled upper {maximum}, "
        f"minimum margin {minimum_margin}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
