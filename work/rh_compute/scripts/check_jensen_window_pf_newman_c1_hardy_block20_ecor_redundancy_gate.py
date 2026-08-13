#!/usr/bin/env python3
"""Validate the exact cubic ecor redundancy artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing ecor gate output")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "exact_displayed_ecor_coefficient_already_present_in_source_gc_phase",
        "ecor status drift",
    )
    require(artifact["identity"]["quartic_match_identity"] == "0", "quartic identity drift")
    require(
        artifact["identity"]["phase_quartic_coefficient"]
        == artifact["identity"]["commented_ecor_coefficient"],
        "displayed ecor mismatch",
    )
    require(artifact["source_audit"]["executable_ecor_occurrence_count"] == 2, "ecor use count drift")
    require(not artifact["source_audit"]["ecor_used_in_phase_expression"], "ecor unexpectedly active")

    for name in ("source",):
        record = artifact[name]
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing {name} source")
        require(file_hash(path) == record["sha256"], f"{name} hash drift")
    for name in ("builder", "checker"):
        record = artifact["sources"][name]
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing {name}")
        require(file_hash(path) == record["sha256"], f"{name} hash drift")

    note = " ".join(NOTE.read_text(encoding="utf-8").split())
    for token in (
        "g(c_src)",
        "9*pi/(32*phi2^5)",
        "double the coefficient",
        "stale comment",
        "does not prove",
    ):
        require(token in note, f"ecor note token missing: {token}")
    print("validated Hardy block-20 ecor redundancy gate: exact quartic match, no executable phase use")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
