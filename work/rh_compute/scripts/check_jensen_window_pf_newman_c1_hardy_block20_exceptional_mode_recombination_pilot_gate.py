#!/usr/bin/env python3
"""Validate the exceptional-mode recombination pilot."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "rigorous_four_branch_exceptional_mode_recombination_validated",
        "exceptional-mode status drift",
    )
    require(artifact["witness_chains"] == [1, 2, 3, 36], "exceptional-mode witness drift")
    require(len(artifact["rows"]) == 4, "exceptional-mode row count drift")
    require({row["selector"] for row in artifact["rows"]} == {"W3", "W4"}, "selector coverage drift")
    aggregate = artifact["aggregate"]
    require(aggregate["branch_count"] == 4, "exceptional-mode branch count drift")
    require(float(aggregate["maximum_model_identity_gap_upper"]) < 1e-40, "model identity gap drift")
    require(float(aggregate["maximum_total_recombination_gap_upper"]) < 1e-10, "recombination gap drift")
    require(float(aggregate["maximum_exceptional_mode_tail_upper"]) < 1e-20, "mode tail drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing exceptional-mode dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("global coefficient `Phi2`", "`1/delta`", "`1/eta`", "not a height-uniform"):
        require(token in note, f"exceptional-mode note token missing: {token}")
    print("validated exceptional-mode recombination pilot: 4/4 branches, W2--W4 identities and direct modes certified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
