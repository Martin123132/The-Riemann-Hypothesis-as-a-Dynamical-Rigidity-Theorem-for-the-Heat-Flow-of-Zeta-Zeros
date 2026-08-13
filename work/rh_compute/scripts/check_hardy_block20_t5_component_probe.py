#!/usr/bin/env python3
"""Validate the source-derived block-20 t5 component probe fixture."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/t5_component_probe/t1e10_block20/probe_result.json"
)
RUNNER = REPO_ROOT / "work/rh_compute/scripts/run_hardy_block20_t5_component_probe.py"
CHECKER = Path(__file__).resolve()
HEX128 = re.compile(r"[0-9A-F]{32}\Z")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, RUNNER, CHECKER):
        require(path.is_file(), f"missing t5 component probe artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "hardy_block20_t5_component_probe", "t5-probe kind drift")
    require(
        artifact["status"] == "source_derived_374_call_probe_with_exact_t5_replay",
        "t5-probe status drift",
    )
    require(artifact["recursive_chain_count"] == 374, "t5-probe chain count drift")
    validation = artifact["validation"]
    require(validation["row_count"] == 374, "t5-probe row count drift")
    require(validation["saved_t5_bit_match_count"] == 374, "t5-probe saved-bit replay drift")
    require(validation["upper_erfc_call_count"] == 374, "t5-probe upper-call count drift")
    require(
        validation["intrinsic_erfc_call_count"]
        == validation["lower_erfc_call_count"] + validation["upper_erfc_call_count"],
        "t5-probe erfc total drift",
    )
    require(
        validation["two_call_row_count"] + validation["three_call_row_count"] == 374,
        "t5-probe slot-row partition drift",
    )
    require(
        validation["intrinsic_erfc_call_count"]
        == 2 * validation["two_call_row_count"] + 3 * validation["three_call_row_count"],
        "t5-probe slot count drift",
    )
    require(all(HEX128.fullmatch(value) for value in validation["constants_hex"].values()), "t5-probe constant hex drift")
    require(artifact["build"]["cpu_cap"] == 1 and artifact["run"]["cpu_cap"] == 1, "t5-probe CPU cap drift")
    require(artifact["build"]["nice"] == 10 and artifact["run"]["nice"] == 10, "t5-probe priority drift")

    for item in artifact["artifacts"].values():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing t5-probe output: {path}")
        require(file_hash(path) == item["sha256"], f"t5-probe output hash drift: {path}")
    for key in ("accepted_source", "telemetry"):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"t5-probe {key} hash drift")
    require(file_hash(RUNNER) == artifact["sources"]["runner"]["sha256"], "t5-probe runner hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "t5-probe checker hash drift")
    require("does not itself prove erfc accuracy" in artifact["proof_boundary"], "t5-probe boundary drift")
    print(
        "validated Hardy block-20 t5 component probe: "
        f"374 calls, {validation['intrinsic_erfc_call_count']} intrinsic erfc values, "
        "374 exact t5 bit replays"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
