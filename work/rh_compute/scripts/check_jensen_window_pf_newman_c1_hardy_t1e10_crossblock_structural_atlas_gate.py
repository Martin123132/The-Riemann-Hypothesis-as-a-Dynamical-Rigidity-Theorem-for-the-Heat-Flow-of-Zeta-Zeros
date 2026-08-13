#!/usr/bin/env python3
"""Validate the t=1e10 Hardy cross-block structural atlas."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing note: {gate.NOTE}")
    result = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(result["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate", "kind drift")
    require(result["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash drift")
    require(result["sources"]["checker"]["sha256"] == gate.file_hash(Path(__file__).resolve()), "checker hash drift")
    for group in ("artifacts", "dependencies"):
        for dependency in result[group].values():
            path = gate.REPO_ROOT / dependency["path"]
            require(path.is_file(), f"missing dependency: {path}")
            require(dependency["sha256"] == gate.file_hash(path), f"dependency hash drift: {path}")

    rows = result["rows"]
    blocks = result["blocks"]
    aggregate = result["aggregate"]
    require(len(rows) == 6784 and len({int(row["chain"]) for row in rows}) == 6784, "row roster drift")
    require([int(row["block"]) for row in blocks] == list(range(20, 36)), "block roster drift")
    require(all(int(row["call_count"]) == 424 for row in blocks), "calls-per-block drift")
    require(sum(int(row["recursive_call_count"]) for row in blocks) == 1414, "recursive count drift")
    require(sum(int(row["direct_call_count"]) for row in blocks) == 5370, "direct count drift")
    require(int(aggregate["call_count"]) == 6784, "aggregate call drift")
    require(int(aggregate["recursive_call_count"]) == 1414, "aggregate recursive drift")
    require(int(aggregate["direct_call_count"]) == 5370, "aggregate direct drift")
    require(int(aggregate["last_recursive_block"]) == 28, "last recursive block drift")
    require(int(aggregate["first_all_direct_block"]) == 29, "all-direct transition drift")
    require(Fraction(aggregate["minimum_recursive_a1_selector_margin"]) > 0, "a1 margin drift")
    require(Fraction(aggregate["minimum_recursive_dual_fractional_margin"]) > 0, "dual margin drift")
    require(Fraction(aggregate["minimum_recursive_curvature_floor"]) > 0, "curvature margin drift")
    require(Fraction(aggregate["minimum_recursive_child_radical_one_plus_z"]) > 0, "radical margin drift")
    require(bool(aggregate["checkpoint_state_equal_ignoring_run_id"]), "checkpoint equivalence drift")
    require(bool(aggregate["displayed_hardy_values_exact"]), "displayed equivalence drift")
    require(result["status"] == "rigorous_t1e10_blocks_20_35_binary_structural_atlas_with_route_transition", "status drift")
    require("not a proof" in gate.NOTE.read_text(encoding="utf-8").lower(), "note proof boundary missing")
    print("validated t=1e10 cross-block structural atlas: 6784 calls, 1414 recursive, blocks 29-35 all direct")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
