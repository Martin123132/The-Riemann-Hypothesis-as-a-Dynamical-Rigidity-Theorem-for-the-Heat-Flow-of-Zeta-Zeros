#!/usr/bin/env python3
"""Validate the four-branch endpoint-linearization decomposition pilot."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "rigorous_four_branch_endpoint_linearization_decomposition_validated",
        "endpoint-linearization pilot status drift",
    )
    require(artifact["precisions_decimal_digits"] == [70, 110], "endpoint-linearization precision drift")
    require(artifact["witness_chains"] == [1, 2, 3, 36], "endpoint-linearization witness drift")
    require(len(artifact["rows"]) == 4, "endpoint-linearization row count drift")
    combos = {(row["phi1_sign"], row["phi3_sign"]) for row in artifact["rows"]}
    require(len(combos) == 4, "endpoint-linearization sign coverage drift")
    for row in artifact["rows"]:
        require(float(row["generic_w5_gap_upper"]) < 1e-40, f"chain {row['chain']} generic/W5 gap drift")
        require(
            float(row["correction_decomposition_gap_upper"]) < 1e-10,
            f"chain {row['chain']} correction gap drift",
        )
    require(float(artifact["aggregate"]["maximum_generic_w5_gap_upper"]) < 1e-40, "aggregate W5 gap drift")
    require(
        float(artifact["aggregate"]["maximum_correction_decomposition_gap_upper"]) < 1e-10,
        "aggregate correction gap drift",
    )
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing endpoint-linearization dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("endpoint-linear model", "Q_exact-Q_paper", "W2 and W3", "not an all-height"):
        require(token in note, f"endpoint-linearization note token missing: {token}")
    print("validated endpoint-linearization decomposition pilot: 4/4 branches, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
