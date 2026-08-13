#!/usr/bin/env python3
"""Validate the t=1e10 cross-block direct-kernel output gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate",
        "kind drift",
    )
    require(len(artifact["rows"]) == 5370, "direct row-count drift")
    require(len(artifact["block_rows"]) == 16, "block row-count drift")
    require(len(artifact["output_rows"]) == 15, "output row-count drift")
    require(sum(row["direct_call_count"] for row in artifact["block_rows"]) == 5370, "block call sum drift")
    require(artifact["local_aggregate"]["precision_overlap_count"] == 5370, "precision overlap drift")
    aggregate = artifact["aggregate"]
    require(aggregate["direct_call_count"] == 5370, "direct aggregate drift")
    require(aggregate["recursive_call_count"] == 1414, "recursive aggregate drift")
    require(aggregate["total_call_count"] == 6784, "total call aggregate drift")
    require(aggregate["recursive_output_anchor_exact"], "recursive output anchor drift")
    require(aggregate["block20_direct_output_anchor_exact"], "block-20 direct output anchor drift")
    for row in artifact["output_rows"]:
        recursive = Fraction(row["recursive_endpoint_majorant_upper"])
        direct = Fraction(row["direct_kernel_majorant_upper"])
        combined = Fraction(row["combined_exact_point_majorant_upper"])
        require(
            abs(recursive + direct - combined) <= Fraction(1, 10**50),
            f"combined output decomposition drift at {row['output_index']}",
        )
        require(row["within_requested_scale"] == (combined < Fraction(5, 1000)), "scale flag drift")
    for record in (*artifact["dependencies"].values(), *artifact["sources"].values()):
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing artifact: {path}")
        require(file_hash(path) == record["sha256"], f"artifact hash drift: {path}")
    require(NOTE.is_file(), "missing markdown note")
    print(
        "validated t=1e10 cross-block direct-kernel output gate: "
        f"5370 direct, combined max {aggregate['maximum_combined_exact_point_majorant_upper']}, "
        f"{aggregate['outputs_within_requested_scale']}/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
