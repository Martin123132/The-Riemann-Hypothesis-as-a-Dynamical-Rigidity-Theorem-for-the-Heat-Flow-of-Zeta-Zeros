#!/usr/bin/env python3
"""Validate the Hardy blocks 20--35 accumulation-weight observer fixture."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import run_hardy_chain_telemetry_fixture as chain


EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/accumulation_weights/t1e10_crossblock"
REFERENCE_ROOT = EXTERNAL_ROOT / "fixtures/chain_telemetry/t1e10_crossblock_probe/enabled"
RESULT = FIXTURE_ROOT / "fixture_result.json"
FIRST_BLOCK = 20
LAST_BLOCK = 35
SUM_COUNT = 212
OUTPUT_COUNT = 15
EXPECTED_ROWS = (LAST_BLOCK - FIRST_BLOCK + 1) * SUM_COUNT * OUTPUT_COUNT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "hardy_crossblock_accumulation_weight_fixture", "kind drift")
    weights = artifact["weights"]
    require(weights["row_count"] == EXPECTED_ROWS, "row-count drift")
    require(weights["block_count"] == 16, "block-count drift")
    require(weights["binary128_payload_count"] == EXPECTED_ROWS * 8, "payload-count drift")
    require(weights["roster_exact"], "roster was not exact")

    for record in artifact["artifacts"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing artifact: {path}")
        require(file_hash(path) == record["sha256"], f"artifact hash drift: {path}")

    weights_path = REPO_ROOT / artifact["artifacts"]["weights"]["path"]
    lines = weights_path.read_text(encoding="ascii").splitlines()
    pattern = re.compile(r"([0-9]+) ([0-9]+) ([0-9]+)(?: ([0-9A-F]{32})){8} ?")
    roster: list[tuple[int, int, int]] = []
    for line_number, line in enumerate(lines, 1):
        require(pattern.fullmatch(line) is not None, f"invalid row {line_number}")
        fields = line.split()
        roster.append((int(fields[0]), int(fields[1]), int(fields[2])))
    expected = [
        (block, sum_index, output_index)
        for block in range(FIRST_BLOCK, LAST_BLOCK + 1)
        for sum_index in range(1, SUM_COUNT + 1)
        for output_index in range(1, OUTPUT_COUNT + 1)
    ]
    require(roster == expected, "weight roster mismatch")

    run_root = FIXTURE_ROOT / "run"
    new_journal = chain.load_jsonl(run_root / "checkpoint.jsonl", decimal=True)
    reference_journal = chain.load_jsonl(REFERENCE_ROOT / "checkpoint.jsonl", decimal=True)
    require(
        chain.without_run_id(new_journal) == chain.without_run_id(reference_journal),
        "checkpoint equivalence drift",
    )
    require(
        chain.parse_grand_totals(run_root / "run.output.txt")
        == chain.parse_grand_totals(REFERENCE_ROOT / "run.output.txt"),
        "displayed Hardy equivalence drift",
    )
    print(f"validated cross-block accumulation weights: {EXPECTED_ROWS} rows, exact evaluator equivalence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
