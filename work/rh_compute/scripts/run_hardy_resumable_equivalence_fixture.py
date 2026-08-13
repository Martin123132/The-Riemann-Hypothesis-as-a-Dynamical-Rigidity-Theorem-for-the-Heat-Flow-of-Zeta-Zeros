#!/usr/bin/env python3
"""Build and exercise the one-rank resumable Hardy diagnostic evaluator."""

from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
RESUMABLE_ROOT = EXTERNAL_ROOT / "resumable"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/resume_equivalence"
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\resumable_20260805"
)
IMAGE = "rh-hardy-fastcode:bookworm"
BINARY_NAME = "zeta14cubicmult.resumable"
PARKED_EXIT_CODES = {75}


CASES = {
    "t1e10": {
        "input": EXTERNAL_ROOT / "fixtures/t1e10/inputs3_et005.nml",
        "accepted": EXTERNAL_ROOT
        / "fixtures/t1e10/zeta14cubicmult_et005_output.txt",
        "stop_block": 5,
        "exercise_all_controls": True,
    },
    "t1e12": {
        "input": EXTERNAL_ROOT / "fixtures/t1e12/inputs3.nml",
        "accepted": EXTERNAL_ROOT
        / "fixtures/t1e12/zeta14cubicmult_et005_output.txt",
        "stop_block": 23,
        "exercise_all_controls": False,
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def docker_mount(path: Path) -> str:
    return str(path.resolve())


def run_command(command: list[str], *, expected: set[int] = {0}) -> subprocess.CompletedProcess[str]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    completed.elapsed_seconds = time.perf_counter() - started  # type: ignore[attr-defined]
    if completed.returncode not in expected:
        joined = " ".join(command)
        raise RuntimeError(
            f"command exited {completed.returncode}, expected {sorted(expected)}:\n"
            f"{joined}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def build(build_root: Path) -> Path:
    build_root.mkdir(parents=True, exist_ok=True)
    binary = build_root / BINARY_NAME
    command = [
        "docker",
        "run",
        "--rm",
        "--cpus=1",
        "-v",
        f"{docker_mount(REPO_ROOT)}:/rh:ro",
        "-v",
        f"{docker_mount(build_root)}:/build",
        "-w",
        "/build",
        IMAGE,
        "bash",
        "-lc",
        (
            "nice -n 10 mpif90 -O3 -z noexecstack "
            "-fallow-argument-mismatch -ffree-line-length-none "
            "/rh/work/rh_compute/external/hardy_fastcode/resumable/"
            "rh_hardy_checkpoint.f90 "
            "/rh/work/rh_compute/external/hardy_fastcode/resumable/"
            "zeta14cubicmult_resumable.f90 -lpari "
            f"-o /build/{BINARY_NAME}"
        ),
    ]
    completed = run_command(command)
    (build_root / "build.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (build_root / "build.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    require(binary.is_file(), "resumable binary was not produced")
    return binary


def inspect_image_id() -> str:
    completed = run_command(
        ["docker", "image", "inspect", "--format={{.Id}}", IMAGE]
    )
    image_id = completed.stdout.strip()
    require(re.fullmatch(r"sha256:[0-9a-f]{64}", image_id) is not None, "Docker image-id format drift")
    return image_id


def remove_fixture_tree(path: Path) -> None:
    resolved = path.resolve()
    fixture_root = FIXTURE_ROOT.resolve()
    require(resolved.is_relative_to(fixture_root), f"refusing to remove outside fixture root: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


def prepare_run(run_root: Path, input_path: Path, *, reset: bool) -> None:
    if reset and run_root.exists():
        remove_fixture_tree(run_root)
    run_root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, run_root / "inputs3.nml")
    (run_root / "zeta14v6resa").write_text("", encoding="ascii")


def provenance_id(binary: Path, run_root: Path, image_id: str) -> str:
    material = {
        "binary_sha256": file_hash(binary),
        "module_sha256": file_hash(RESUMABLE_ROOT / "rh_hardy_checkpoint.f90"),
        "source_sha256": file_hash(
            RESUMABLE_ROOT / "zeta14cubicmult_resumable.f90"
        ),
        "input_sha256": file_hash(run_root / "inputs3.nml"),
        "image_id": image_id,
        "mpi_ranks": 1,
        "rs_chunk_size": 257,
    }
    encoded = json.dumps(material, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def invoke(
    binary: Path,
    run_root: Path,
    image_id: str,
    label: str,
    *,
    resume: bool,
    extra_environment: dict[str, str] | None = None,
    parked: bool = False,
    expected_return_codes: set[int] | None = None,
) -> dict:
    (run_root / "zeta14v6resa").write_text("", encoding="ascii")
    environment = {
        "RH_HARDY_CHECKPOINT": "/run/checkpoint.dat",
        "RH_HARDY_JOURNAL": "/run/checkpoint.jsonl",
        "RH_HARDY_STOP_FILE": "/run/stop.request",
        "RH_HARDY_RESUME": "1" if resume else "0",
        "RH_HARDY_RS_CHUNK_SIZE": "257",
        "RH_HARDY_RUN_ID": provenance_id(binary, run_root, image_id),
    }
    if extra_environment:
        environment.update(extra_environment)

    command = ["docker", "run", "--rm", "--cpus=1"]
    for key, value in environment.items():
        command.extend(["-e", f"{key}={value}"])
    command.extend(
        [
            "-e",
            "OMPI_ALLOW_RUN_AS_ROOT=1",
            "-e",
            "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
            "-v",
            f"{docker_mount(binary.parent)}:/build:ro",
            "-v",
            f"{docker_mount(run_root)}:/run",
            "-w",
            "/run",
            IMAGE,
            "bash",
            "-lc",
            f"nice -n 10 mpirun --allow-run-as-root -np 1 /build/{binary.name}",
        ]
    )
    expected = expected_return_codes or ({0} | (PARKED_EXIT_CODES if parked else set()))
    completed = run_command(command, expected=expected)
    if parked:
        require(completed.returncode in PARKED_EXIT_CODES, f"{label} did not park")
    elif expected_return_codes is None:
        require(completed.returncode == 0, f"{label} did not complete")
    (run_root / f"{label}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_root / f"{label}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if not parked:
        shutil.copy2(run_root / "zeta14v6resa", run_root / f"{label}.output.txt")
    return {
        "label": label,
        "resume": resume,
        "returncode": completed.returncode,
        "elapsed_seconds": completed.elapsed_seconds,  # type: ignore[attr-defined]
        "environment": environment,
    }


def load_journal(path: Path) -> list[dict]:
    records = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        require(line.strip(), f"blank JSONL record at line {line_number}")
        records.append(json.loads(line, parse_float=Decimal))
    require(records, f"empty journal: {path}")
    return records


def validate_journal(records: list[dict]) -> dict:
    units = [(int(record["stage"]), int(record["unit"])) for record in records]
    require(len(units) == len(set(units)), "duplicate stage/unit checkpoint")
    require(units[0] == (0, 0), "journal must begin with block zero")
    require(units[-1] == (2, 0), "journal must end with completed state")
    run_ids = {record["run_id"] for record in records}
    require(len(run_ids) == 1, "journal provenance hash drift")
    require(re.fullmatch(r"[0-9a-f]{64}", next(iter(run_ids))) is not None, "run-id format drift")
    for record in records:
        require(record["format"] == "RH_HARDY_CHECKPOINT_V1", "journal format drift")
        require(record["comm_size"] == 1, "fixture communicator size drift")
        require(record["nchalf"] == 15, "fixture shift count drift")
        require(record["rs_chunk_size"] == 257, "fixture RS chunk-size drift")
        round_trip_digits = math.ceil(
            1 + record["real_digits"] * math.log10(record["real_radix"])
        )
        require(
            record["serialized_significand_digits"] >= round_trip_digits,
            "checkpoint precision is insufficient for exact round trip",
        )
        for field in ("zsum", "rszsumtot", "rszsum"):
            values = record[field]
            require(len(values) == 15, f"{field} channel count drift")
            require(
                all(isinstance(value, (int, Decimal)) and Decimal(value).is_finite() for value in values),
                f"{field} type or finiteness drift",
            )
    zp_units = [unit for stage, unit in units if stage == 0]
    require(zp_units == list(range(records[0]["totblock"] + 1)), "Gaussian block sequence gap")
    rs_units = [unit for stage, unit in units if stage == 1]
    require(rs_units == list(range(1, len(rs_units) + 1)), "RS unit sequence gap")
    require(len(rs_units) == records[-2]["nmax"] + 1, "RS unit count drift")
    return {
        "record_count": len(records),
        "zp_record_count": len(zp_units),
        "rs_record_count": len(rs_units),
        "complete_record_count": sum(stage == 2 for stage, _ in units),
        "totblock": records[0]["totblock"],
        "nmax": records[-2]["nmax"],
        "run_id": next(iter(run_ids)),
    }


def parse_grand_totals(path: Path) -> list[Decimal]:
    values = re.findall(
        r"Grand total of Hardy function Z\(t[^=]*=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == 15, f"expected 15 final Hardy values in {path}")
    return [Decimal(value) for value in values]


def run_case(binary: Path, image_id: str, case_name: str, case: dict) -> dict:
    case_root = FIXTURE_ROOT / case_name
    uninterrupted_root = case_root / "uninterrupted"
    resumed_root = case_root / "resumed"

    prepare_run(uninterrupted_root, case["input"], reset=True)
    uninterrupted_runs = [
        invoke(binary, uninterrupted_root, image_id, "uninterrupted", resume=False)
    ]

    prepare_run(resumed_root, case["input"], reset=True)
    resumed_runs: list[dict] = []
    if case["exercise_all_controls"]:
        resumed_runs.append(
            invoke(
                binary,
                resumed_root,
                image_id,
                "park_runtime",
                resume=False,
                extra_environment={"RH_HARDY_MAX_SECONDS": "0.000000001"},
                parked=True,
            )
        )
        resumed_runs.append(
            invoke(
                binary,
                resumed_root,
                image_id,
                "park_unit",
                resume=True,
                extra_environment={
                    "RH_HARDY_STOP_AFTER_STAGE": "0",
                    "RH_HARDY_STOP_AFTER_UNIT": str(case["stop_block"]),
                },
                parked=True,
            )
        )
        (resumed_root / "stop.request").write_text("stop after next unit\n", encoding="ascii")
        resumed_runs.append(
            invoke(binary, resumed_root, image_id, "park_stop_file", resume=True, parked=True)
        )
        (resumed_root / "stop.request").unlink()
    else:
        resumed_runs.append(
            invoke(
                binary,
                resumed_root,
                image_id,
                "park_unit",
                resume=False,
                extra_environment={
                    "RH_HARDY_STOP_AFTER_STAGE": "0",
                    "RH_HARDY_STOP_AFTER_UNIT": str(case["stop_block"]),
                },
                parked=True,
            )
        )

    resumed_runs.append(
        invoke(
            binary,
            resumed_root,
            image_id,
            "park_rs",
            resume=True,
            extra_environment={
                "RH_HARDY_STOP_AFTER_STAGE": "1",
                "RH_HARDY_STOP_AFTER_UNIT": "1",
            },
            parked=True,
        )
    )
    resumed_runs.append(
        invoke(binary, resumed_root, image_id, "resumed_complete", resume=True)
    )

    uninterrupted_journal = load_journal(uninterrupted_root / "checkpoint.jsonl")
    resumed_journal = load_journal(resumed_root / "checkpoint.jsonl")
    uninterrupted_summary = validate_journal(uninterrupted_journal)
    resumed_summary = validate_journal(resumed_journal)
    require(uninterrupted_journal == resumed_journal, "full-precision journal state drift")

    uninterrupted_values = parse_grand_totals(
        uninterrupted_root / "uninterrupted.output.txt"
    )
    resumed_values = parse_grand_totals(resumed_root / "resumed_complete.output.txt")
    accepted_values = parse_grand_totals(case["accepted"])
    require(uninterrupted_values == resumed_values, "resumed final output drift")
    accepted_max_error = max(
        abs(left - right) for left, right in zip(uninterrupted_values, accepted_values)
    )

    negative_runs = []
    if case["exercise_all_controls"]:
        for label, overrides in (
            ("reject_run_id", {"RH_HARDY_RUN_ID": "0" * 64}),
            ("reject_chunk_size", {"RH_HARDY_RS_CHUNK_SIZE": "258"}),
        ):
            negative_root = case_root / label
            if negative_root.exists():
                remove_fixture_tree(negative_root)
            shutil.copytree(resumed_root, negative_root)
            before_snapshot = file_hash(negative_root / "checkpoint.dat")
            before_journal = file_hash(negative_root / "checkpoint.jsonl")
            negative_runs.append(
                invoke(
                    binary,
                    negative_root,
                    image_id,
                    label,
                    resume=True,
                    extra_environment=overrides,
                    expected_return_codes={2},
                )
            )
            require(
                file_hash(negative_root / "checkpoint.dat") == before_snapshot,
                f"{label} changed the authoritative snapshot",
            )
            require(
                file_hash(negative_root / "checkpoint.jsonl") == before_journal,
                f"{label} changed the append-only journal",
            )

    return {
        "case": case_name,
        "uninterrupted_runs": uninterrupted_runs,
        "resumed_runs": resumed_runs,
        "uninterrupted_journal": uninterrupted_summary,
        "resumed_journal": resumed_summary,
        "journal_exact_match": True,
        "final_output_exact_match": True,
        "accepted_output_max_abs_error": str(accepted_max_error),
        "negative_resume_tests": negative_runs,
        "hashes": {
            "input": file_hash(case["input"]),
            "accepted_output": file_hash(case["accepted"]),
            "uninterrupted_journal": file_hash(uninterrupted_root / "checkpoint.jsonl"),
            "resumed_journal": file_hash(resumed_root / "checkpoint.jsonl"),
            "uninterrupted_output": file_hash(
                uninterrupted_root / "uninterrupted.output.txt"
            ),
            "resumed_output": file_hash(resumed_root / "resumed_complete.output.txt"),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    parser.add_argument("--case", choices=["all", *CASES], default="all")
    args = parser.parse_args()

    image_id = inspect_image_id()
    binary = build(args.build_root)
    selected = CASES if args.case == "all" else {args.case: CASES[args.case]}
    case_results = [
        run_case(binary, image_id, name, case) for name, case in selected.items()
    ]
    result = {
        "kind": "hardy_resumable_equivalence_fixture",
        "image": IMAGE,
        "image_id": image_id,
        "mpi_ranks": 1,
        "cpu_limit": 1,
        "binary_path": str(binary),
        "binary_sha256": file_hash(binary),
        "module_sha256": file_hash(RESUMABLE_ROOT / "rh_hardy_checkpoint.f90"),
        "source_sha256": file_hash(
            RESUMABLE_ROOT / "zeta14cubicmult_resumable.f90"
        ),
        "cases": case_results,
    }
    FIXTURE_ROOT.mkdir(parents=True, exist_ok=True)
    result_path = FIXTURE_ROOT / "fixture_result.json"
    temporary = result_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, result_path)
    print(
        "validated resumable Hardy fixtures: "
        f"{len(case_results)} cases, exact journals and final outputs"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
