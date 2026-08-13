#!/usr/bin/env python3
"""Build the telemetry derivative and capture every cubic chain in block 20."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import run_hardy_chain_telemetry_fixture as base


FIXTURE_ROOT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full"
)
RUN_ROOT = FIXTURE_ROOT / "enabled"
RESULT = FIXTURE_ROOT / "fixture_result.json"
REFERENCE_ROOT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/enabled"
)
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\telemetry_block20_full_20260806"
)
MAX_CHAINS = 424


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    args = parser.parse_args()

    source_check = subprocess.run(
        [sys.executable, str(base.SOURCE_CHECKER)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    require(source_check.returncode == 0, f"telemetry source contract failed: {source_check.stderr}")
    image_id = base.inspect_image_id()
    binary, build = base.build(args.build_root)

    base.FIXTURE_ROOT = FIXTURE_ROOT
    base.reset_run_root(RUN_ROOT)
    run = base.invoke(binary, RUN_ROOT, image_id, telemetry=True, max_chains=MAX_CHAINS)
    telemetry = base.validate_telemetry(RUN_ROOT / "chain_telemetry.jsonl", MAX_CHAINS)
    require(telemetry["chain_count"] == MAX_CHAINS, "full block telemetry stopped early")

    records = base.load_jsonl(RUN_ROOT / "chain_telemetry.jsonl")
    chains = [record for record in records if record.get("type") == "chain"]
    expected_roster = [
        (sum_index, branch)
        for sum_index in range(1, 213)
        for branch in (1, 2)
    ]
    actual_roster = [(int(record["sum_index"]), int(record["branch"])) for record in chains]
    require(actual_roster == expected_roster, "block-20 chain roster gap")
    require(all(int(record["block"]) == 20 for record in chains), "telemetry escaped block 20")
    require(all(int(record["base_degree"]) == 3 for record in chains), "noncubic block-20 chain")

    new_journal = base.load_jsonl(RUN_ROOT / "checkpoint.jsonl", decimal=True)
    reference_journal = base.load_jsonl(REFERENCE_ROOT / "checkpoint.jsonl", decimal=True)
    require(
        base.without_run_id(new_journal) == base.without_run_id(reference_journal),
        "full telemetry changed evaluator checkpoint state",
    )
    displayed_values = base.parse_grand_totals(RUN_ROOT / "run.output.txt")
    reference_values = base.parse_grand_totals(REFERENCE_ROOT / "run.output.txt")
    require(displayed_values == reference_values, "full telemetry changed displayed Hardy values")

    chain_lengths: dict[str, int] = {}
    for record in chains:
        key = str(int(record["mit"]))
        chain_lengths[key] = chain_lengths.get(key, 0) + 1

    artifact = {
        "kind": "hardy_block20_full_chain_telemetry_fixture",
        "status": "exact_424_chain_block20_fixture_with_evaluator_equivalence_validated",
        "max_chains": MAX_CHAINS,
        "block": 20,
        "sum_count": 212,
        "branch_count": 424,
        "sum_index_range": [1, 212],
        "branch_roster": [1, 2],
        "mit_histogram": chain_lengths,
        "source_contract_output": source_check.stdout.strip(),
        "image_id": image_id,
        "build": build,
        "run": run,
        "telemetry_validation": telemetry,
        "equivalence": {
            "checkpoint_records": len(new_journal),
            "checkpoint_state_equal_ignoring_run_id": True,
            "displayed_hardy_values_exact": True,
            "displayed_value_count": len(displayed_values),
            "timing_text_excluded": True,
            "reference_journal_sha256": file_hash(REFERENCE_ROOT / "checkpoint.jsonl"),
            "reference_output_sha256": file_hash(REFERENCE_ROOT / "run.output.txt"),
        },
        "artifacts": {
            "telemetry": {
                "path": relative(RUN_ROOT / "chain_telemetry.jsonl"),
                "sha256": file_hash(RUN_ROOT / "chain_telemetry.jsonl"),
            },
            "checkpoint": {
                "path": relative(RUN_ROOT / "checkpoint.jsonl"),
                "sha256": file_hash(RUN_ROOT / "checkpoint.jsonl"),
            },
            "output": {
                "path": relative(RUN_ROOT / "run.output.txt"),
                "sha256": file_hash(RUN_ROOT / "run.output.txt"),
            },
            "input": {
                "path": relative(RUN_ROOT / "inputs3.nml"),
                "sha256": file_hash(RUN_ROOT / "inputs3.nml"),
            },
        },
        "sources": {
            "accepted_source_sha256": file_hash(base.RESUMABLE_ROOT / "zeta14cubicmult_resumable.f90"),
            "checkpoint_module_sha256": file_hash(base.RESUMABLE_ROOT / "rh_hardy_checkpoint.f90"),
            "telemetry_module_sha256": file_hash(base.TELEMETRY_ROOT / "rh_hardy_chain_telemetry.f90"),
            "telemetry_source_sha256": file_hash(base.TELEMETRY_ROOT / "zeta14cubicmult_telemetry.f90"),
            "runner_sha256": file_hash(Path(__file__).resolve()),
            "base_runner_sha256": file_hash(Path(base.__file__).resolve()),
        },
        "proof_boundary": (
            "Finite exact-bit source telemetry for one low-height block only. It proves observer "
            "equivalence for this run, not special-function accuracy, a uniform recurrence theorem, "
            "physical-height control, RH, or a prize-level result."
        ),
    }
    FIXTURE_ROOT.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 full telemetry fixture: "
        f"{MAX_CHAINS} chains, {telemetry['recurrence_count']} recurrences, exact state/displayed values"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
