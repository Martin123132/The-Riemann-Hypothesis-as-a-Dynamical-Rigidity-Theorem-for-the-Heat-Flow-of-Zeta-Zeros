#!/usr/bin/env python3
"""Validate hashed one-cell sparse-H23 order-eleven pilot artifacts."""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments as driver  # noqa: E402


def load_record_sources(artifact: dict) -> dict[Fraction, dict]:
    paths = [REPO_ROOT / artifact["sparse_source_path"]]
    paths.extend(
        REPO_ROOT / item["path"]
        for item in artifact.get("additional_pilot_sources", [])
    )
    records = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                record = json.loads(line)
                target = Fraction(record["target_t"])
                if target in records and records[target] != record:
                    raise RuntimeError(f"conflicting source rows at t={target}")
                records[target] = record
    return records


def validate_artifact(path: Path) -> list[str]:
    issues = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if (
        artifact.get("kind") != "order11_sparse_h23_lower_bridge_pilot"
        or artifact.get("status")
        != "rigorous cell-local interval pilot over hashed inputs"
        or artifact.get("passed") is not True
    ):
        issues.append("identity or pass status changed")
    proof_boundary = artifact.get("proof_boundary", "")
    for phrase in (
        "does not certify the full lower bridge",
        "Lambda<=0",
        "RH",
    ):
        if phrase not in proof_boundary:
            issues.append(f"proof boundary misses {phrase!r}")

    try:
        cell_left = Fraction(artifact["cell_left"])
        blocks = artifact["blocks"]
        if len(blocks) != 2:
            issues.append(f"expected two quarter blocks, found {len(blocks)}")
        expected = (
            (cell_left, cell_left + Fraction(1, 4)),
            (cell_left + Fraction(1, 4), cell_left + Fraction(1, 2)),
        )
        for index, (block, (left, right)) in enumerate(zip(blocks, expected)):
            if (
                Fraction(block["anchor"]) != left
                or Fraction(block["right"]) != right
                or Fraction(block["width"]) != Fraction(1, 4)
                or block.get("model_degrees") != [16, 15, 14]
                or block.get("maximum_h_derivative_order") != 17
                or block.get("passed") is not True
            ):
                issues.append(f"block {index} geometry or model contract changed")
            scaled = Decimal(block["scaled_curvature_upper"])
            margin = Decimal(block["curvature_margin_lower"])
            if not scaled < Decimal(6000) or not margin > 0:
                issues.append(f"block {index} does not prove the curvature wall")
        if artifact.get("scaled_curvature_upper") != [
            block["scaled_curvature_upper"] for block in blocks
        ]:
            issues.append("pilot summary differs from full block records")
    except (KeyError, ValueError, TypeError) as exc:
        issues.append(f"block parse failure: {exc}")

    immutable_sources = (
        artifact.get("h_sources", [])
        + artifact.get("code_sources", [])
        + artifact.get("additional_pilot_sources", [])
    )
    if len(artifact.get("h_sources", [])) != 4 or len(
        artifact.get("code_sources", [])
    ) != 4:
        issues.append("immutable H/core source list changed")
    for source in immutable_sources:
        source_path = REPO_ROOT / source["path"]
        if not source_path.exists() or driver.sha256(source_path) != source["sha256"]:
            issues.append(f"immutable source hash changed: {source.get('path')}")

    try:
        records = load_record_sources(artifact)
        selected = artifact["selected_anchor_records"]
        selected_targets = []
        for item in selected:
            target = Fraction(item["target_t"])
            selected_targets.append(target)
            if target not in records:
                issues.append(f"selected anchor t={target} is unavailable")
            elif driver.record_sha256(records[target]) != item["sha256"]:
                issues.append(f"selected anchor hash changed at t={target}")
        required_targets = [
            expansion + shift
            for expansion in (cell_left, cell_left + Fraction(1, 2))
            for shift in range(-9, 10)
        ]
        for target in required_targets:
            if min(abs(target - anchor) for anchor in selected_targets) > 1:
                issues.append(f"selected anchors do not cover target t={target}")
                break
    except (KeyError, OSError, RuntimeError, ValueError) as exc:
        issues.append(f"anchor-source validation failed: {exc}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", type=Path, nargs="+")
    args = parser.parse_args()
    issues = []
    for path in args.artifacts:
        try:
            issues.extend(
                f"{path.name}: {issue}" for issue in validate_artifact(path)
            )
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            issues.append(f"{path.name}: artifact parse failure: {exc}")
    if issues:
        print(f"order-eleven sparse H23 lower-bridge pilots: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-eleven sparse H23 lower-bridge pilots: "
        f"{len(args.artifacts)} cells, {2 * len(args.artifacts)} quarter blocks, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
