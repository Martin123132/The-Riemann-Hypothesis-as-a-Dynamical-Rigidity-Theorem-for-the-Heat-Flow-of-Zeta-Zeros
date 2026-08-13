#!/usr/bin/env python3
"""Validate the all-call endpoint-linearization decomposition certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["status"] == "rigorous_all_374_endpoint_linearization_decomposition_validated", "full decomposition status drift")
    require(artifact["precisions_decimal_digits"] == [70, 110], "full decomposition precision drift")
    require(len(artifact["rows"]) == 374, "full decomposition row count drift")
    require(len({int(row["chain"]) for row in artifact["rows"]}) == 374, "full decomposition duplicate chain")
    aggregate = artifact["aggregate"]
    require(aggregate["call_count"] == 374, "full decomposition call count drift")
    require(aggregate["generic_w5_identity_count"] == 374, "full generic/W5 identity count drift")
    require(aggregate["correction_decomposition_identity_count"] == 374, "full correction identity count drift")
    require(float(aggregate["maximum_generic_w5_gap_upper"]) < 1e-40, "full generic/W5 gap drift")
    require(float(aggregate["maximum_correction_decomposition_gap_upper"]) < 1e-10, "full correction gap drift")
    require(float(aggregate["maximum_component_triangle_sum_upper"]) > 0, "full triangle sum drift")
    require(float(aggregate["maximum_component_triangle_to_correction_ratio_upper"]) >= 1, "full triangle ratio drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing full decomposition dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("all 374", "endpoint_half+a_endpoint", "triangle/correction ratio", "not a fitted"):
        require(token in note, f"full decomposition note token missing: {token}")
    print("validated endpoint-linearization decomposition full gate: 374/374 calls, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
