#!/usr/bin/env python3
"""Independently validate compact adaptive-H23 order-eleven segments."""

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

import jensen_window_pf_compound_order11_compact_adaptive_h23_segments as source  # noqa: E402


def validate_run_contract_dict(existing: dict) -> list[str]:
    try:
        canonical = source.canonical_run_contract()
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        return [f"canonical run-contract construction failed: {exc}"]
    if existing != canonical:
        return ["run contract differs from canonical source and code hashes"]
    return []


def validate_run_contract(path: Path) -> list[str]:
    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"run-contract read failed: {exc}"]
    return validate_run_contract_dict(existing)


def expected_point_metadata(expansion: Fraction) -> dict:
    targets = [
        expansion + shift
        for shift in range(-source.STENCIL_RADIUS, source.STENCIL_RADIUS + 1)
    ]
    anchors = [source.anchor_for_target(target) for target in targets]
    return {
        "point_source_targets": [str(target) for target in targets],
        "point_source_anchors": [str(anchor) for anchor in anchors],
        "point_source_profiles": [
            source.profile_for_anchor(anchor) for anchor in anchors
        ],
        "maximum_actual_propagation_distance": str(
            max(abs(target - anchor) for target, anchor in zip(targets, anchors))
        ),
    }


def validate_block(
    block: dict,
    *,
    expected_left: Fraction,
    expected_expansion: Fraction,
    label: str,
) -> list[str]:
    issues = []
    expected_right = expected_left + source.QUARTER_WIDTH
    metadata = expected_point_metadata(expected_expansion)
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
            or block.get("regime") != "compact"
            or block.get("model_degrees") != list(source.DERIVATIVE_MODEL_DEGREES)
            or block.get("maximum_h_derivative_order") != source.MODEL_MAXIMUM_H_ORDER
            or block.get("stable_taylor_surplus") != source.STABLE_TAYLOR_SURPLUS
            or block.get("passed") is not True
        ):
            issues.append(f"{label}: geometry or model contract changed")
        for key, expected in metadata.items():
            if block.get(key) != expected:
                issues.append(f"{label}: {key} changed")
        scaled = Decimal(block["scaled_curvature_upper"])
        margin = Decimal(block["curvature_margin_lower"])
        remainder = Decimal(block["final_uniform_remainder_upper"])
        if not scaled < Decimal(source.CURVATURE_CONSTANT):
            issues.append(f"{label}: scaled curvature upper reaches the wall")
        if not margin > 0:
            issues.append(f"{label}: curvature margin is not positive")
        if not remainder >= 0:
            issues.append(f"{label}: final remainder is negative")
        quadrature = block.get("quadrature", {})
        if (
            quadrature.get("shift_count") != 19
            or len(quadrature.get("mode_brackets", [])) != 19
        ):
            issues.append(f"{label}: shifted quadrature contract changed")
        stages = block.get("stage_diagnostics", [])
        if [stage.get("stage") for stage in stages] != list(range(2, 10)):
            issues.append(f"{label}: stable-stage diagnostics changed")
    except (KeyError, ValueError, TypeError) as exc:
        issues.append(f"{label}: parse failure: {exc}")
    return issues


def validate_segment_record(
    record: dict,
    task: tuple[int, Fraction, Fraction],
) -> tuple[list[str], list[Decimal], list[Decimal]]:
    index, segment_left, segment_right = task
    label = f"segment {index}"
    issues = []
    blocks = record.get("blocks", [])
    expected_blocks = int((segment_right - segment_left) / source.QUARTER_WIDTH)
    if (
        record.get("kind") != "order11_compact_adaptive_h23_segment"
        or record.get("index") != index
        or record.get("segment_left") != str(segment_left)
        or record.get("segment_right") != str(segment_right)
        or record.get("block_count") != expected_blocks
        or len(blocks) != expected_blocks
        or record.get("passed") is not True
    ):
        issues.append(f"{label}: identity or block count changed")
    scaled_values = []
    margins = []
    expected_profiles = set()
    for block_index, block in enumerate(blocks):
        expected_left = segment_left + block_index * source.QUARTER_WIDTH
        cell_left = segment_left + (block_index // 2) * source.CELL_WIDTH
        expected_expansion = (
            cell_left if block_index % 2 == 0 else cell_left + source.CELL_WIDTH
        )
        issues.extend(
            validate_block(
                block,
                expected_left=expected_left,
                expected_expansion=expected_expansion,
                label=f"{label} block {block_index}",
            )
        )
        expected_profiles.update(
            expected_point_metadata(expected_expansion)["point_source_profiles"]
        )
        try:
            scaled_values.append(Decimal(block["scaled_curvature_upper"]))
            margins.append(Decimal(block["curvature_margin_lower"]))
        except (KeyError, ValueError):
            pass
        if len(issues) >= 20:
            break
    if record.get("source_profiles") != sorted(expected_profiles):
        issues.append(f"{label}: source profile summary changed")
    if scaled_values:
        maximum = max(scaled_values)
        maximum_block = blocks[scaled_values.index(maximum)]
        if (
            Decimal(record["largest_scaled_curvature_upper"]) != maximum
            or record.get("largest_scaled_curvature_anchor")
            != maximum_block.get("anchor")
        ):
            issues.append(f"{label}: largest scaled curvature summary changed")
    if margins:
        minimum = min(margins)
        minimum_block = blocks[margins.index(minimum)]
        if (
            Decimal(record["smallest_margin_lower"]) != minimum
            or record.get("smallest_margin_anchor") != minimum_block.get("anchor")
        ):
            issues.append(f"{label}: smallest margin summary changed")
    return issues, scaled_values, margins


def validate_segment_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
) -> tuple[list[str], int, int, Decimal | None, Decimal | None]:
    issues = []
    segment_count = 0
    block_count = 0
    global_maximum = None
    global_minimum = None
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
            segment_issues, scaled, margins = validate_segment_record(
                record,
                tasks[segment_count],
            )
            issues.extend(segment_issues)
            if scaled:
                local_maximum = max(scaled)
                global_maximum = (
                    local_maximum
                    if global_maximum is None
                    else max(global_maximum, local_maximum)
                )
            if margins:
                local_minimum = min(margins)
                global_minimum = (
                    local_minimum
                    if global_minimum is None
                    else min(global_minimum, local_minimum)
                )
            block_count += len(record.get("blocks", []))
            segment_count += 1
            if len(issues) >= 20:
                issues.append("validation stopped after 20 issues")
                break
    return issues, segment_count, block_count, global_maximum, global_minimum


def validate_pilot(path: Path) -> tuple[list[str], int, int, Decimal | None, Decimal | None]:
    try:
        artifact = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"pilot read failed: {exc}"], 0, 0, None, None
    issues = []
    contract = artifact.get("run_contract", {})
    issues.extend(validate_run_contract_dict(contract))
    if artifact.get("run_contract_sha256") != source.record_sha256(contract):
        issues.append("pilot run-contract digest changed")
    indices = artifact.get("segment_indices", [])
    records = artifact.get("segments", [])
    if (
        artifact.get("kind") != "order11_compact_adaptive_h23_segment_pilot"
        or len(indices) != len(records)
        or len(set(indices)) != len(indices)
    ):
        issues.append("pilot identity or segment index list changed")
    tasks = source.deterministic_segments()
    scaled = []
    margins = []
    for index, record in zip(indices, records):
        if not isinstance(index, int) or not 0 <= index < len(tasks):
            issues.append(f"pilot segment index is invalid: {index}")
            continue
        segment_issues, segment_scaled, segment_margins = validate_segment_record(
            record,
            tasks[index],
        )
        issues.extend(segment_issues)
        scaled.extend(segment_scaled)
        margins.extend(segment_margins)
    summary = artifact.get("summary", {})
    if (
        summary.get("segments") != len(records)
        or summary.get("quarter_blocks") != sum(
            record.get("block_count", 0) for record in records
        )
        or summary.get("failed_segments") != 0
    ):
        issues.append("pilot count summary changed")
    if scaled and Decimal(summary.get("largest_scaled_curvature_upper")) != max(scaled):
        issues.append("pilot maximum summary changed")
    if margins and Decimal(summary.get("smallest_margin_lower")) != min(margins):
        issues.append("pilot minimum margin summary changed")
    return (
        issues,
        len(records),
        sum(record.get("block_count", 0) for record in records),
        max(scaled) if scaled else None,
        min(margins) if margins else None,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_SEGMENT_CACHE)
    parser.add_argument("--run-contract", type=Path, default=source.DEFAULT_RUN_CONTRACT)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    if args.pilot is not None:
        issues, segments, blocks, maximum, minimum = validate_pilot(args.pilot)
        label = "compact adaptive H23 pilot"
    else:
        issues = validate_run_contract(args.run_contract)
        tasks = source.deterministic_segments()
        if not args.cache.exists():
            issues.append(f"missing segment cache: {args.cache}")
            segments = blocks = 0
            maximum = minimum = None
        else:
            try:
                segment_issues, segments, blocks, maximum, minimum = (
                    validate_segment_cache(args.cache, tasks)
                )
                issues.extend(segment_issues)
            except (OSError, RuntimeError, ValueError) as exc:
                issues.append(f"segment-cache validation failed: {exc}")
                segments = blocks = 0
                maximum = minimum = None
        if args.require_complete and segments != len(tasks):
            issues.append(f"expected {len(tasks)} segments, found {segments}")
        label = "compact adaptive H23 segment cache"
    if issues:
        print(f"{label}: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        f"validated {label}: {segments} segments, {blocks} quarter blocks, "
        f"maximum scaled upper {maximum}, minimum margin {minimum}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
