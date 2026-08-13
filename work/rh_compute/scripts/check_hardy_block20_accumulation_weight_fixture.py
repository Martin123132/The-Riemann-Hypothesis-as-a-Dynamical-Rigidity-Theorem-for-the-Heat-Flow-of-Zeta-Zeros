#!/usr/bin/env python3
"""Validate the block-20 accumulation-weight observer fixture."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import run_hardy_chain_telemetry_fixture as base


EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/accumulation_weights/t1e10_block20"
REFERENCE_ROOT = EXTERNAL_ROOT / "fixtures/chain_telemetry/t1e10_block20_full/enabled"
RESULT = FIXTURE_ROOT / "fixture_result.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "hardy_block20_accumulation_weight_fixture", "kind drift")
    require(artifact["block"] == 20, "block drift")
    require(artifact["weights"]["row_count"] == 3180, "row-count drift")
    require(artifact["weights"]["binary128_payload_count"] == 25440, "payload-count drift")
    require(artifact["weights"]["roster_exact"], "roster was not exact")

    for record in artifact["artifacts"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing artifact: {path}")
        require(file_hash(path) == record["sha256"], f"artifact hash drift: {path}")

    weights_path = REPO_ROOT / artifact["artifacts"]["weights"]["path"]
    lines = weights_path.read_text(encoding="ascii").splitlines()
    pattern = re.compile(r"([0-9]+) ([0-9]+)(?: ([0-9A-F]{32})){8} ?")
    roster: list[tuple[int, int]] = []
    for line_number, line in enumerate(lines, 1):
        require(pattern.fullmatch(line) is not None, f"invalid row {line_number}")
        fields = line.split()
        roster.append((int(fields[0]), int(fields[1])))
    expected = [(sum_index, output_index) for sum_index in range(1, 213) for output_index in range(1, 16)]
    require(roster == expected, "weight roster mismatch")

    run_root = FIXTURE_ROOT / "run"
    new_journal = base.load_jsonl(run_root / "checkpoint.jsonl", decimal=True)
    reference_journal = base.load_jsonl(REFERENCE_ROOT / "checkpoint.jsonl", decimal=True)
    require(
        base.without_run_id(new_journal) == base.without_run_id(reference_journal),
        "checkpoint equivalence drift",
    )
    require(
        base.parse_grand_totals(run_root / "run.output.txt")
        == base.parse_grand_totals(REFERENCE_ROOT / "run.output.txt"),
        "displayed Hardy equivalence drift",
    )
    print("validated block-20 accumulation weights: 3180 rows, exact evaluator equivalence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
