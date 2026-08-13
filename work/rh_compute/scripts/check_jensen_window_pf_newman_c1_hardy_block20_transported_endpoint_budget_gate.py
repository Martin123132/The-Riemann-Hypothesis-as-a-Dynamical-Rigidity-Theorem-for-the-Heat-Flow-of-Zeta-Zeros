#!/usr/bin/env python3
"""Validate the transported block-20 endpoint budget gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as gate


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate",
        "kind drift",
    )
    require(len(artifact["rows"]) == 15, "output-row drift")
    require(artifact["aggregate"]["total_call_count"] == 424, "call-count drift")
    require(artifact["aggregate"]["recursive_call_count"] == 374, "recursive-count drift")
    require(artifact["aggregate"]["direct_call_count"] == 50, "direct-count drift")
    require(
        artifact["aggregate"]["outputs_with_endpoint_majorant_within_requested_scale"] == 0,
        "absolute-budget classification drift",
    )
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    require(file_hash(gate.BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(Path(__file__)) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    rows, aggregate = gate.compute()
    require(rows == artifact["rows"], "independent transported rows drift")
    require(aggregate == artifact["aggregate"], "independent aggregate drift")
    require(
        all(Fraction(row["transported_endpoint_majorant_upper"]) > gate.REQUESTED_ERROR_SCALE for row in rows),
        "an output unexpectedly closes under absolute propagation",
    )
    require(
        all(
            Fraction(row["transported_direct_majorant_upper"])
            < Fraction(row["transported_endpoint_majorant_upper"]) * Fraction(1, 10**20)
            for row in rows
        ),
        "direct-kernel budget is no longer negligible",
    )
    print(
        "validated transported endpoint budget: 424 calls, 15 outputs, "
        "absolute majorant exceeds 0.005 on every output"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
