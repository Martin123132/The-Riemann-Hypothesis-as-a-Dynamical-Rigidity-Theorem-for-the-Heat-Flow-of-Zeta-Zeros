#!/usr/bin/env python3
"""Validate canonical sparse-H23 order-twelve lower-bridge segments or pilots."""

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

import jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments as source  # noqa: E402


def validate_preflight() -> tuple[list[str], dict]:
    issues = []
    tasks = source.deterministic_segments()
    half_cells = 0
    quarter_blocks = 0
    targets = set()
    anchors = set()
    maximum_distance = Fraction(0)
    collar_left = None
    collar_right = None
    for _, segment_left, segment_right in tasks:
        cell_left = segment_left
        while cell_left < segment_right:
            cell_right = cell_left + source.CELL_WIDTH
            half_cells += 1
            quarter_blocks += 2
            for expansion in (cell_left, cell_right):
                for shift in range(-source.STENCIL_RADIUS, source.STENCIL_RADIUS + 1):
                    target = expansion + shift
                    anchor = source.anchor_for_target(target)
                    targets.add(target)
                    anchors.add(anchor)
                    maximum_distance = max(maximum_distance, abs(target - anchor))
                    collar_left = (
                        min(target, anchor)
                        if collar_left is None
                        else min(collar_left, target, anchor)
                    )
                    collar_right = (
                        max(target, anchor)
                        if collar_right is None
                        else max(collar_right, target, anchor)
                    )
            cell_left = cell_right
    if (
        len(tasks) != source.EXPECTED_SEGMENTS
        or half_cells != source.EXPECTED_HALF_CELLS
        or quarter_blocks != source.EXPECTED_QUARTER_BLOCKS
        or tasks[0][1] != source.START_T
        or tasks[-1][2] != source.END_T
    ):
        issues.append("deterministic lower partition changed")
    if (
        min(targets) != Fraction(1493)
        or max(targets) != Fraction(5710)
        or maximum_distance != source.MAXIMUM_PROPAGATION_DISTANCE
        or source.anchor_for_target(Fraction(5710)) != Fraction(5708)
    ):
        issues.append("stencil target or sparse-anchor geometry changed")
    if (
        collar_left != Fraction(1492)
        or collar_right != Fraction(5710)
        or source.inherited.base.BROAD_H_START > collar_left
        or source.inherited.base.BROAD_H_END != Fraction(5709)
        or source.compact_h_source.DEFAULT_START_T > Fraction(5709)
        or source.compact_h_source.DEFAULT_END_T < collar_right
    ):
        issues.append("mixed H24 collar or handoff geometry changed")
    try:
        source.validate_broad_h_source()
        sparse_manifest = source.inherited.validate_sparse_manifest()
        compact_manifest = source.validate_compact_h_source()
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as exc:
        issues.append(f"source-manifest preflight failed: {exc}")
        sparse_manifest = {"cache": {}}
        compact_manifest = {"cache": {}}
    return issues, {
        "segments": len(tasks),
        "half_cells": half_cells,
        "quarter_blocks": quarter_blocks,
        "stencil_targets": len(targets),
        "target_range": [str(min(targets)), str(max(targets))],
        "sparse_anchors_used": len(anchors),
        "anchor_range": [str(min(anchors)), str(max(anchors))],
        "maximum_propagation_distance": str(maximum_distance),
        "mixed_h24_collar": [str(collar_left), str(collar_right)],
        "h24_handoff": str(source.inherited.base.BROAD_H_END),
        "sparse_rows": sparse_manifest.get("cache", {}).get("row_count"),
        "compact_h_rows": compact_manifest.get("cache", {}).get("row_count"),
    }


def validate_run_contract_dict(existing: dict) -> list[str]:
    try:
        canonical = source.canonical_run_contract()
    except (OSError, RuntimeError, ValueError) as exc:
        return [f"canonical run-contract validation failed: {exc}"]
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
        "point_source_contract_ids": [
            source.sparse_builder.ROW_CONTRACT for _ in anchors
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
            or block.get("regime") != "lower"
            or block.get("model_degrees") != list(source.DERIVATIVE_MODEL_DEGREES)
            or block.get("maximum_h_derivative_order") != source.MODEL_MAXIMUM_H_ORDER
            or block.get("stable_taylor_surplus") != source.STABLE_TAYLOR_SURPLUS
            or block.get("passed") is not True
        ):
            issues.append(f"{label}: geometry or model contract changed")
        for key, expected in metadata.items():
            if block.get(key) != expected:
                issues.append(f"{label}: {key} changed")
        kinds = block.get("h_source_kinds", [])
        if not kinds or not set(kinds) <= {"broad_tenth", "compact_unit"}:
            issues.append(f"{label}: H-source handoff metadata changed")
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
            quadrature.get("shift_count") != 21
            or len(quadrature.get("mode_brackets", [])) != 21
        ):
            issues.append(f"{label}: shifted quadrature contract changed")
        stages = block.get("stage_diagnostics", [])
        if [stage.get("stage") for stage in stages] != list(range(2, 11)):
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
        record.get("kind") != "order12_sparse_h23_lower_bridge_segment"
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
    expected_contracts = set()
    expected_h_kinds = set()
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
        expected_contracts.update(
            expected_point_metadata(expected_expansion)["point_source_contract_ids"]
        )
        expected_h_kinds.update(block.get("h_source_kinds", []))
        try:
            scaled_values.append(Decimal(block["scaled_curvature_upper"]))
            margins.append(Decimal(block["curvature_margin_lower"]))
        except (KeyError, ValueError):
            pass
        if len(issues) >= 20:
            break
    if record.get("point_source_contract_ids") != sorted(expected_contracts):
        issues.append(f"{label}: point-source contract summary changed")
    if record.get("h_source_kinds") != sorted(expected_h_kinds):
        issues.append(f"{label}: H-source summary changed")
    if scaled_values:
        maximum = max(scaled_values)
        maximum_block = blocks[scaled_values.index(maximum)]
        if (
            Decimal(record["largest_scaled_curvature_upper"]) != maximum
            or record.get("largest_scaled_curvature_anchor") != maximum_block.get("anchor")
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
    segments = 0
    blocks = 0
    maximum = None
    minimum_margin = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if segments >= len(tasks):
                issues.append("segment cache has too many rows")
                break
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                issues.append(f"line {line_number}: invalid JSON: {exc}")
                break
            record_issues, scaled, margins = validate_segment_record(
                record,
                tasks[segments],
            )
            issues.extend(record_issues)
            blocks += len(record.get("blocks", []))
            if scaled:
                local_maximum = max(scaled)
                maximum = (
                    local_maximum if maximum is None else max(maximum, local_maximum)
                )
            if margins:
                local_minimum = min(margins)
                minimum_margin = (
                    local_minimum
                    if minimum_margin is None
                    else min(minimum_margin, local_minimum)
                )
            segments += 1
            if len(issues) >= 20:
                issues.append("validation stopped after 20 issues")
                break
    return issues, segments, blocks, maximum, minimum_margin


def validate_pilot(path: Path) -> tuple[list[str], int, int, Decimal | None, Decimal | None]:
    try:
        artifact = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"pilot read failed: {exc}"], 0, 0, None, None
    issues = validate_run_contract_dict(artifact.get("run_contract", {}))
    if (
        artifact.get("kind") != "order12_sparse_h23_lower_bridge_pilot"
        or artifact.get("run_contract_sha256")
        != source.record_sha256(artifact.get("run_contract", {}))
        or artifact.get("summary", {}).get("failed_segments") != 0
    ):
        issues.append("pilot identity or contract hash changed")
    tasks = source.deterministic_segments()
    indices = artifact.get("segment_indices", [])
    records = artifact.get("segments", [])
    if len(indices) != len(records) or indices != sorted(set(indices)):
        issues.append("pilot segment index list changed")
    maximum = None
    minimum_margin = None
    block_count = 0
    for index, record in zip(indices, records):
        if not isinstance(index, int) or not 0 <= index < len(tasks):
            issues.append(f"invalid pilot segment index: {index}")
            continue
        record_issues, scaled, margins = validate_segment_record(record, tasks[index])
        issues.extend(record_issues)
        block_count += len(record.get("blocks", []))
        if scaled:
            local_maximum = max(scaled)
            maximum = local_maximum if maximum is None else max(maximum, local_maximum)
        if margins:
            local_minimum = min(margins)
            minimum_margin = (
                local_minimum
                if minimum_margin is None
                else min(minimum_margin, local_minimum)
            )
    summary = artifact.get("summary", {})
    if (
        summary.get("segments") != len(records)
        or summary.get("quarter_blocks") != block_count
        or (maximum is not None and Decimal(summary["largest_scaled_curvature_upper"]) != maximum)
        or (minimum_margin is not None and Decimal(summary["smallest_margin_lower"]) != minimum_margin)
    ):
        issues.append("pilot summary changed")
    return issues, len(records), block_count, maximum, minimum_margin


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_SEGMENT_CACHE)
    parser.add_argument("--run-contract", type=Path, default=source.DEFAULT_RUN_CONTRACT)
    parser.add_argument("--pilot", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    if args.preflight_only:
        if args.pilot is not None or args.require_complete:
            parser.error("--preflight-only cannot be combined with cache or pilot promotion")
        issues, summary = validate_preflight()
        if issues:
            print(f"order-twelve sparse H23 lower preflight: {len(issues)} issues")
            for issue in issues:
                print(f"- {issue}")
            return 1
        print(
            "validated order-twelve sparse H23 lower preflight: "
            f"{summary['segments']} segments, {summary['quarter_blocks']} quarter blocks, "
            f"{summary['stencil_targets']} stencil targets, maximum distance "
            f"{summary['maximum_propagation_distance']}, H24 handoff "
            f"{summary['h24_handoff']}, 0 issues"
        )
        return 0
    if args.pilot is not None:
        issues, segments, blocks, maximum, minimum_margin = validate_pilot(args.pilot)
        label = "pilot"
    else:
        tasks = source.deterministic_segments()
        issues = validate_run_contract(args.run_contract)
        if not args.cache.exists():
            issues.append(f"missing segment cache: {args.cache}")
            segments = blocks = 0
            maximum = minimum_margin = None
        else:
            try:
                cache_issues, segments, blocks, maximum, minimum_margin = (
                    validate_segment_cache(args.cache, tasks)
                )
                issues.extend(cache_issues)
            except (OSError, RuntimeError, ValueError) as exc:
                issues.append(f"segment-cache validation failed: {exc}")
                segments = blocks = 0
                maximum = minimum_margin = None
        if args.require_complete and segments != len(tasks):
            issues.append(f"expected {len(tasks)} segments, found {segments}")
        label = "complete" if segments == len(tasks) else "valid resumable prefix"

    if issues:
        print(f"order-twelve sparse H23 lower bridge: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        f"validated order-twelve sparse H23 lower bridge ({label}): "
        f"{segments} segments, {blocks} quarter blocks, "
        f"maximum scaled upper {maximum}, minimum margin {minimum_margin}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
