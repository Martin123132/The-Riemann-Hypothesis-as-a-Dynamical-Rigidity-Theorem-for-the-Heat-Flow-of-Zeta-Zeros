#!/usr/bin/env python3
"""Validate the blocks 22--28 signed bridge geometry pilot."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), f"missing later-block bridge pilot: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing later-block bridge note: {gate.NOTE}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate", "later-block pilot kind drift")
    require(artifact["status"] == "all_new_later_block_exact_point_geometries_closed", "later-block pilot status drift")
    require(bool(artifact["passed"]), "later-block bridge pilot did not pass")
    rows = artifact["rows"]
    require([int(row["chain"]) for row in rows] == list(gate.WITNESS_CHAINS), "later-block pilot witness drift")
    require(len(rows) == gate.EXPECTED_GROUP_COUNT, "later-block pilot row count drift")
    require(all(22 <= int(row["block"]) <= 28 and int(row["child_length"]) == 1 for row in rows), "later-block pilot roster drift")
    require(all(bool(row["precision_overlap"]) and bool(row["signed_bridge_closed"]) for row in rows), "later-block pilot closure drift")
    require(all(bool(row["exact_correction_within_retained_majorant"]) for row in rows), "later-block pilot retained majorant drift")
    aggregate = artifact["aggregate"]
    require(aggregate["parent_lengths"] == [135, 154, 175, 200, 228, 259, 295], "later-block parent lengths drift")
    require(aggregate["blocks"] == list(range(22, 29)), "later-block block coverage drift")
    require(int(aggregate["geometry_group_count"]) == gate.EXPECTED_GROUP_COUNT, "later-block geometry count drift")
    require(int(aggregate["closed_witness_count"]) == gate.EXPECTED_GROUP_COUNT, "later-block closure count drift")
    require(int(aggregate["exact_corrections_within_retained_count"]) == gate.EXPECTED_GROUP_COUNT, "later-block retained count drift")
    require(Fraction(aggregate["maximum_formula_target_gap_upper"]) < Fraction(1, 10**12), "later-block formula gap regression")
    provenance = artifact["provenance"]
    for name, expected in (
        ("telemetry", gate.TELEMETRY),
        ("endpoints", gate.ENDPOINTS),
        ("weights", gate.WEIGHTS),
        ("bridge_evaluator", gate.bridge.BUILDER),
        ("builder", gate.BUILDER),
        ("checker", gate.CHECKER),
    ):
        require(provenance[name]["sha256"] == file_hash(expected), f"later-block pilot {name} hash drift")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("135", "295", "18", "not a proof"):
        require(token in note, f"later-block note token missing: {token}")
    print(
        "validated blocks22--28 bridge geometry pilot: "
        f"witnesses={len(rows)}, max_gap={aggregate['maximum_formula_target_gap_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
