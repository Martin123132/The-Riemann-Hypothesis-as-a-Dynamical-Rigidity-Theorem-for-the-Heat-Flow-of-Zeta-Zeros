#!/usr/bin/env python3
"""Validate the pinned binary128 rounding-mode probe."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/binary128_rounding_mode_probe/probe_result.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file(), "missing binary128 rounding-mode result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "hardy_binary128_rounding_mode_probe", "rounding probe kind drift")
    require(artifact["status"] == "pinned_source_environment_starts_binary128_in_round_to_nearest_mode", "rounding probe status drift")
    require(bool(artifact["passed"]), "rounding probe failed")
    values = artifact["probe_values"]
    expected = {
        "radix": "2",
        "digits": "113",
        "min_exponent": "-16381",
        "max_exponent": "16384",
        "ieee_datatype": "T",
        "ieee_nearest_supported": "T",
        "rounding_mode": "nearest",
        "one_hex": "3FFF0000000000000000000000000000",
        "epsilon_hex": "3F8F0000000000000000000000000000",
        "tiny_hex": "00010000000000000000000000000000",
    }
    for key, value in expected.items():
        require(values.get(key) == value, f"rounding probe {key} drift")
    require(all(int(count) == 0 for count in artifact["source_rounding_mode_mutator_counts"].values()), "accepted source changes rounding mode")
    require(artifact["build"]["cpu_cap"] == 1 and artifact["build"]["container_nice"] == 10, "rounding probe resource drift")
    for record in artifact["artifacts"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing rounding probe artifact {path}")
        require(file_hash(path) == record["sha256"], f"rounding probe hash drift {path}")
    print("validated binary128 rounding-mode probe: radix=2, digits=113, mode=nearest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
