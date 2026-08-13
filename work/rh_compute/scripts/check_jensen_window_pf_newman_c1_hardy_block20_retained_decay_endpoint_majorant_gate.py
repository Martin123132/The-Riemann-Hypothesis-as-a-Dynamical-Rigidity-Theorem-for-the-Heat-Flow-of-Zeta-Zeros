#!/usr/bin/env python3
"""Validate retained-decay endpoint majorants and transported closure."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate",
        "kind drift",
    )
    require(len(artifact["rows"]) == 374, "call roster drift")
    require(len(artifact["transported_outputs"]) == 15, "output roster drift")
    require(artifact["aggregate"]["precision_overlap_count"] == 374, "precision overlap drift")
    require(artifact["aggregate"]["w3_count"] == 165, "W3 count drift")
    require(artifact["aggregate"]["w4_count"] == 209, "W4 count drift")
    require(artifact["aggregate"]["outputs_within_requested_scale"] == 15, "transport closure drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(
        all(
            Fraction(row["complete_majorant_lower"]) > 0
            and Fraction(row["complete_improvement_factor_lower"]) > 1
            and bool(row["precision_overlap"])
            for row in artifact["rows"]
        ),
        "local retained-decay row invariant failed",
    )
    require(
        all(
            Fraction(row["retained_decay_transported_majorant_upper"]) < Fraction(5, 1000)
            and Fraction(row["transport_improvement_factor"]) > 1
            and bool(row["within_requested_scale"])
            for row in artifact["transported_outputs"]
        ),
        "transported retained-decay invariant failed",
    )
    print("validated retained-decay endpoint majorants: 374 calls, 15/15 transported outputs below 0.005")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
