#!/usr/bin/env python3
"""Validate the full block-20 Hardy chain telemetry fixture."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full"
RESULT = ROOT / "fixture_result.json"
RUNNER = REPO_ROOT / "work/rh_compute/scripts/run_hardy_block20_full_telemetry_fixture.py"
CHECKER = Path(__file__).resolve()
ACCEPTED_SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, RUNNER, CHECKER):
        require(path.is_file(), f"missing full-block fixture artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "hardy_block20_full_chain_telemetry_fixture", "fixture kind drift")
    require(
        artifact["status"] == "exact_424_chain_block20_fixture_with_evaluator_equivalence_validated",
        "fixture status drift",
    )
    require(artifact["max_chains"] == 424, "fixture cap drift")
    require(artifact["block"] == 20 and artifact["sum_count"] == 212, "fixture roster drift")
    require(artifact["branch_count"] == 424, "fixture branch count drift")
    require(artifact["sum_index_range"] == [1, 212], "fixture sum-index drift")
    require(artifact["branch_roster"] == [1, 2], "fixture branch roster drift")
    require(sum(artifact["mit_histogram"].values()) == 424, "fixture MIT histogram drift")
    require(artifact["telemetry_validation"]["chain_count"] == 424, "telemetry validation count drift")
    require(artifact["equivalence"]["checkpoint_state_equal_ignoring_run_id"] is True, "checkpoint equivalence drift")
    require(artifact["equivalence"]["displayed_hardy_values_exact"] is True, "displayed-value equivalence drift")
    require(artifact["equivalence"]["displayed_value_count"] == 15, "displayed-value count drift")
    require(artifact["equivalence"]["timing_text_excluded"] is True, "timing-text boundary drift")
    require(artifact["sources"]["accepted_source_sha256"] == ACCEPTED_SOURCE_SHA256, "accepted source drift")
    require(file_hash(RUNNER) == artifact["sources"]["runner_sha256"], "runner hash drift")

    for key in ("telemetry", "checkpoint", "output", "input"):
        record = artifact["artifacts"][key]
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing full-block {key}")
        require(file_hash(path) == record["sha256"], f"full-block {key} hash drift")

    telemetry_path = REPO_ROOT / artifact["artifacts"]["telemetry"]["path"]
    chains = []
    for line in telemetry_path.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        if record.get("type") == "chain":
            chains.append(record)
    require(len(chains) == 424, "full-block chain record count drift")
    require([int(record["chain"]) for record in chains] == list(range(1, 425)), "chain sequence drift")
    expected = [(index, branch) for index in range(1, 213) for branch in (1, 2)]
    actual = [(int(record["sum_index"]), int(record["branch"])) for record in chains]
    require(actual == expected, "full-block physical roster gap")
    require(all(int(record["block"]) == 20 for record in chains), "chain escaped block 20")

    require("Finite exact-bit source telemetry" in artifact["proof_boundary"], "fixture boundary drift")
    print(
        "validated Hardy block-20 full telemetry fixture: "
        f"424 chains, {artifact['telemetry_validation']['recurrence_count']} recurrences, "
        "exact state/displayed values"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
