#!/usr/bin/env python3
"""Validate the t=1e10 cross-block retained endpoint gate."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing note: {gate.NOTE}")
    result = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(result["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate", "kind drift")
    require(result["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash drift")
    require(result["sources"]["checker"]["sha256"] == gate.file_hash(Path(__file__).resolve()), "checker hash drift")
    for group in ("dependencies", "sources"):
        for name, dependency in result[group].items():
            if group == "sources" and name in ("builder", "checker"):
                continue
            path = gate.REPO_ROOT / dependency["path"]
            require(path.is_file(), f"missing dependency: {path}")
            require(dependency["sha256"] == gate.file_hash(path), f"dependency hash drift: {path}")
    aggregate = result["aggregate"]
    rows = result["rows"]
    blocks = result["blocks"]
    require(len(rows) == 1414 and len({int(row["chain"]) for row in rows}) == 1414, "row roster drift")
    require([int(row["block"]) for row in blocks] == list(range(20, 29)), "block roster drift")
    require(sum(int(row["recursive_call_count"]) for row in blocks) == 1414, "block call drift")
    require(int(aggregate["recursive_call_count"]) == 1414, "aggregate call drift")
    require(int(aggregate["newly_evaluated_call_count"]) == 1040, "new call drift")
    require(int(aggregate["reused_block20_call_count"]) == 374, "block20 reuse drift")
    require(int(aggregate["precision_overlap_count"]) == 1414, "precision overlap drift")
    require(int(aggregate["w3_count"]) + int(aggregate["w4_count"]) == 1414, "selector count drift")
    require(Fraction(aggregate["maximum_complete_majorant_upper"]) > 0, "maximum majorant drift")
    require(all(bool(row["precision_overlap"]) for row in rows), "row precision flag drift")
    require(result["status"] == "rigorous_retained_decay_endpoint_majorants_for_all_1414_recursive_calls_in_blocks_20_28", "status drift")
    require("not a proof" in gate.NOTE.read_text(encoding="utf-8").lower(), "note proof boundary missing")
    print(
        "validated t=1e10 cross-block retained endpoint gate: "
        f"1414 calls, blocks 20-28, max {aggregate['maximum_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
