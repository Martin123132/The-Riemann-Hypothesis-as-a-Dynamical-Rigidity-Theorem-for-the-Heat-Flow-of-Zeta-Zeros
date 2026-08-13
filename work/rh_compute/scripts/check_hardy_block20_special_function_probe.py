#!/usr/bin/env python3
"""Validate the source-derived block-20 PSI/ERF probe fixture."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/special_function_probe/t1e10_block20/probe_result.json"
)
RUNNER = REPO_ROOT / "work/rh_compute/scripts/run_hardy_block20_special_function_probe.py"
CHECKER = Path(__file__).resolve()
HEX128 = re.compile(r"[0-9A-F]{32}\Z")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, RUNNER, CHECKER):
        require(path.is_file(), f"missing special-function probe artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "hardy_block20_special_function_probe", "probe kind drift")
    require(
        artifact["status"] == "source_derived_374_call_probe_with_exact_t1_t2_t4_replay",
        "probe status drift",
    )
    require(artifact["recursive_chain_count"] == 374, "probe chain count drift")
    validation = artifact["validation"]
    require(validation["row_count"] == 374, "probe row count drift")
    require(validation["psi_call_count"] == 2244, "probe PSI count drift")
    require(validation["erf_call_count"] == 748, "probe ERF count drift")
    require(validation["saved_component_bit_match_count"] == 1122, "probe saved-bit replay drift")
    require(all(HEX128.fullmatch(value) for value in validation["constants_hex"].values()), "probe constant hex drift")
    require(artifact["build"]["cpu_cap"] == 1 and artifact["run"]["cpu_cap"] == 1, "probe CPU cap drift")
    require(artifact["build"]["nice"] == 10 and artifact["run"]["nice"] == 10, "probe priority drift")

    for item in artifact["artifacts"].values():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing probe output: {path}")
        require(file_hash(path) == item["sha256"], f"probe output hash drift: {path}")
    for key in ("accepted_source", "telemetry"):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"probe {key} hash drift")
    require(file_hash(RUNNER) == artifact["sources"]["runner"]["sha256"], "probe runner hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "probe checker hash drift")
    require("does not itself prove special-function accuracy" in artifact["proof_boundary"], "probe boundary drift")
    print(
        "validated Hardy block-20 special-function probe: "
        "374 calls, 2992 source values, 1122 exact t1/t2/t4 bit replays"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
