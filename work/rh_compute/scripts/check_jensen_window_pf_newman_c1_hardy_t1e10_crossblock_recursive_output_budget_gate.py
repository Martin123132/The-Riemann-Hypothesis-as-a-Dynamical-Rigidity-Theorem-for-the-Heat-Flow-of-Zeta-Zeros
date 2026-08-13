#!/usr/bin/env python3
"""Validate the t=1e10 cross-block recursive output budget gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate",
        "kind drift",
    )
    require(len(artifact["block_rows"]) == 9, "block-row count drift")
    require(len(artifact["output_rows"]) == 15, "output-row count drift")
    aggregate = artifact["aggregate"]
    require(aggregate["recursive_call_count"] == 1414, "recursive-call count drift")
    require(aggregate["block20_output_anchor_exact"], "block-20 anchor was not exact")
    require(sum(row["recursive_call_count"] for row in artifact["block_rows"]) == 1414, "block call sum drift")
    for row in artifact["output_rows"]:
        block20 = Fraction(row["block20_majorant_upper"])
        later = Fraction(row["later_blocks_majorant_upper"])
        total = Fraction(row["all_recursive_blocks_majorant_upper"])
        require(
            abs(block20 + later - total) <= Fraction(1, 10**50),
            f"output decomposition drift at {row['output_index']}",
        )
        require(row["within_requested_scale"] == (total < Fraction(5, 1000)), "scale flag drift")
    for record in (*artifact["dependencies"].values(), *artifact["sources"].values()):
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing artifact: {path}")
        require(file_hash(path) == record["sha256"], f"artifact hash drift: {path}")
    require(NOTE.is_file(), "missing markdown note")
    print(
        "validated t=1e10 cross-block recursive output budget: "
        f"1414 calls, max {aggregate['maximum_all_recursive_blocks_majorant_upper']}, "
        f"{aggregate['outputs_within_requested_scale']}/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
