#!/usr/bin/env python3
"""Build and exercise the bounded one-rank Hardy chain telemetry scout."""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
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
TELEMETRY_ROOT = EXTERNAL_ROOT / "telemetry"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/chain_telemetry/t1e10"
INPUT = EXTERNAL_ROOT / "fixtures/t1e10/inputs3_et005.nml"
UPSTREAM_OUTPUT = EXTERNAL_ROOT / "fixtures/t1e10/zeta14cubicmult_et005_output.txt"
REFERENCE_ROOT = EXTERNAL_ROOT / "fixtures/resume_equivalence/t1e10/uninterrupted"
REFERENCE_JOURNAL = REFERENCE_ROOT / "checkpoint.jsonl"
REFERENCE_OUTPUT = REFERENCE_ROOT / "uninterrupted.output.txt"
SOURCE_CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_hardy_chain_telemetry_source_contract.py"
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\telemetry_20260806"
)
IMAGE = "rh-hardy-fastcode:bookworm"
BINARY_NAME = "zeta14cubicmult.telemetry"
ACCEPTED_SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def docker_mount(path: Path) -> str:
    return str(path.resolve())


def run_command(
    command: list[str], *, expected: set[int] = {0}
) -> subprocess.CompletedProcess[str]:
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
        raise RuntimeError(
            f"command exited {completed.returncode}, expected {sorted(expected)}:\n"
            f"{' '.join(command)}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def inspect_image_id() -> str:
    completed = run_command(
        ["docker", "image", "inspect", "--format={{.Id}}", IMAGE]
    )
    image_id = completed.stdout.strip()
    require(
        re.fullmatch(r"sha256:[0-9a-f]{64}", image_id) is not None,
        "Docker image-id format drift",
    )
    return image_id


def build(build_root: Path) -> tuple[Path, dict]:
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
            "/rh/work/rh_compute/external/hardy_fastcode/telemetry/"
            "rh_hardy_chain_telemetry.f90 "
            "/rh/work/rh_compute/external/hardy_fastcode/telemetry/"
            "zeta14cubicmult_telemetry.f90 -lpari "
            f"-o /build/{BINARY_NAME}"
        ),
    ]
    completed = run_command(command)
    (build_root / "build.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (build_root / "build.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    require(binary.is_file(), "telemetry binary was not produced")
    return binary, {
        "elapsed_seconds": completed.elapsed_seconds,  # type: ignore[attr-defined]
        "binary_sha256": file_hash(binary),
        "stdout_sha256": file_hash(build_root / "build.stdout.txt"),
        "stderr_sha256": file_hash(build_root / "build.stderr.txt"),
    }


def reset_run_root(path: Path) -> None:
    resolved = path.resolve()
    fixture_root = FIXTURE_ROOT.resolve()
    require(
        resolved.is_relative_to(fixture_root),
        f"refusing to reset outside telemetry fixture root: {resolved}",
    )
    if resolved.exists():
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    shutil.copy2(INPUT, resolved / "inputs3.nml")
    (resolved / "zeta14v6resa").write_text("", encoding="ascii")


def provenance_id(binary: Path, image_id: str) -> str:
    material = {
        "binary_sha256": file_hash(binary),
        "checkpoint_module_sha256": file_hash(RESUMABLE_ROOT / "rh_hardy_checkpoint.f90"),
        "telemetry_module_sha256": file_hash(TELEMETRY_ROOT / "rh_hardy_chain_telemetry.f90"),
        "telemetry_source_sha256": file_hash(TELEMETRY_ROOT / "zeta14cubicmult_telemetry.f90"),
        "accepted_source_sha256": file_hash(RESUMABLE_ROOT / "zeta14cubicmult_resumable.f90"),
        "input_sha256": file_hash(INPUT),
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
    *,
    telemetry: bool,
    max_chains: int,
) -> dict:
    environment = {
        "RH_HARDY_CHECKPOINT": "/run/checkpoint.dat",
        "RH_HARDY_JOURNAL": "/run/checkpoint.jsonl",
        "RH_HARDY_STOP_FILE": "/run/stop.request",
        "RH_HARDY_RESUME": "0",
        "RH_HARDY_RS_CHUNK_SIZE": "257",
        "RH_HARDY_RUN_ID": provenance_id(binary, image_id),
    }
    if telemetry:
        environment.update(
            {
                "RH_HARDY_CHAIN_TELEMETRY": "/run/chain_telemetry.jsonl",
                "RH_HARDY_TELEMETRY_MAX_CHAINS": str(max_chains),
                "RH_HARDY_TELEMETRY_MAX_DIRECT_TERMS": "100000",
            }
        )

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
    completed = run_command(command)
    (run_root / "run.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_root / "run.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    shutil.copy2(run_root / "zeta14v6resa", run_root / "run.output.txt")
    return {
        "telemetry_enabled": telemetry,
        "elapsed_seconds": completed.elapsed_seconds,  # type: ignore[attr-defined]
        "returncode": completed.returncode,
        "environment": environment,
        "journal_sha256": file_hash(run_root / "checkpoint.jsonl"),
        "output_sha256": file_hash(run_root / "run.output.txt"),
    }


def load_jsonl(path: Path, *, decimal: bool = False) -> list[dict]:
    records: list[dict] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        require(line.strip(), f"blank JSONL record at {path}:{line_number}")
        if decimal:
            records.append(json.loads(line, parse_float=Decimal))
        else:
            records.append(json.loads(line))
    require(records, f"empty JSONL file: {path}")
    return records


def without_run_id(records: list[dict]) -> list[dict]:
    return [{key: value for key, value in record.items() if key != "run_id"} for record in records]


def parse_grand_totals(path: Path) -> list[Decimal]:
    values = re.findall(
        r"Grand total of Hardy function Z\(t[^=]*=\s*([-+0-9.Ee]+)",
        path.read_text(encoding="utf-8"),
    )
    require(len(values) == 15, f"expected 15 final Hardy values in {path}")
    return [Decimal(value) for value in values]


def complex_decimal(value: list[str]) -> tuple[Decimal, Decimal]:
    require(isinstance(value, list) and len(value) == 2, "invalid telemetry complex value")
    result = (Decimal(value[0]), Decimal(value[1]))
    require(all(part.is_finite() for part in result), "non-finite telemetry complex value")
    return result


def cadd(*values: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return (sum(value[0] for value in values), sum(value[1] for value in values))


def cconj(value: tuple[Decimal, Decimal]) -> tuple[Decimal, Decimal]:
    return (value[0], -value[1])


def validate_telemetry(path: Path, max_chains: int) -> dict:
    records = load_jsonl(path)
    require(records[0]["type"] == "header", "telemetry header missing")
    require(records[0]["schema"] == "rh_hardy_chain_telemetry_v1", "telemetry schema drift")
    require(records[-1]["type"] == "summary", "telemetry summary missing")

    chains = [record for record in records if record["type"] == "chain"]
    levels = [record for record in records if record["type"] == "level"]
    q_terms = [record for record in records if record["type"] == "q_terms"]
    recurrences = [record for record in records if record["type"] == "recurrence"]
    ends = [record for record in records if record["type"] == "chain_end"]
    chain_ids = [record["chain"] for record in chains]

    require(0 < len(chains) <= max_chains, "telemetry chain bound failure")
    require(chain_ids == list(range(1, len(chains) + 1)), "telemetry chain sequence gap")
    require(len(ends) == len(chains), "telemetry chain-end count drift")
    require(records[-1]["chains_recorded"] == len(chains), "telemetry summary count drift")
    require(all(record["direct_ok"] for record in recurrences), "bounded direct sum was skipped")

    expected_levels = sum(int(record["mit"]) for record in chains)
    expected_steps = sum(int(record["mit"]) - 1 for record in chains)
    require(len(levels) == expected_levels, "telemetry level count drift")
    require(len(q_terms) == expected_steps, "telemetry q-term count drift")
    require(len(recurrences) == expected_steps, "telemetry recurrence count drift")

    hex128 = re.compile(r"[0-9A-F]{32}")

    def require_hex(value: object, label: str) -> None:
        require(
            isinstance(value, str) and hex128.fullmatch(value) is not None,
            f"invalid binary128 payload: {label}",
        )

    def require_complex_hex(value: object, label: str) -> None:
        require(isinstance(value, list) and len(value) == 2, f"invalid complex payload: {label}")
        require_hex(value[0], f"{label}.real")
        require_hex(value[1], f"{label}.imag")

    for record in chains:
        require_hex(record.get("tpm_hex"), f"chain {record['chain']} tpm")
        initial_hex = record.get("initial_coefficients_hex")
        require(
            isinstance(initial_hex, list)
            and len(initial_hex) == int(record["base_degree"]),
            "initial coefficient hex/degree mismatch",
        )
        for index, value in enumerate(initial_hex, 1):
            require_hex(value, f"chain {record['chain']} initial coefficient {index}")
    for record in levels:
        require_hex(record.get("phi0_hex"), f"level {record['chain']}/{record['level']} phi0")
        require_hex(record.get("xr_hex"), f"level {record['chain']}/{record['level']} xr")
        require_hex(
            record.get("frac_length_hex"),
            f"level {record['chain']}/{record['level']} frac_length",
        )
        coefficient_hex = record.get("coefficients_hex")
        require(
            isinstance(coefficient_hex, list)
            and len(coefficient_hex) == int(record["degree"]),
            "coefficient hex/degree mismatch",
        )
        for index, value in enumerate(coefficient_hex, 1):
            require_hex(value, f"level {record['chain']}/{record['level']} coefficient {index}")
    for record in q_terms:
        for field in ("t1", "t2", "t3", "t4", "t5", "endpoint", "qq"):
            require_complex_hex(record.get(f"{field}_hex"), f"q {record['chain']}/{record['nit']} {field}")
    for record in recurrences:
        for field in (
            "raw_parent",
            "raw_child",
            "adapted_parent",
            "adapted_child",
            "state_before",
            "multiplier",
            "qq",
            "state_pre_transform",
            "state_after",
            "model_pre_transform",
            "model_after",
            "local_defect",
            "child_state_defect",
            "parent_state_defect",
        ):
            require_complex_hex(
                record.get(f"{field}_hex"),
                f"recurrence {record['chain']}/{record['nit']} {field}",
            )

    q_by_key = {(record["chain"], record["nit"]): record for record in q_terms}
    max_q_identity_error = Decimal(0)
    max_q_identity_scaled_error = Decimal(0)
    max_logged_qq_disagreement = Decimal(0)
    local_defect_norms: list[Decimal] = []
    child_state_defect_norms: list[Decimal] = []
    parent_state_defect_norms: list[Decimal] = []

    with localcontext() as context:
        context.prec = 70
        for recurrence in recurrences:
            key = (recurrence["chain"], recurrence["nit"])
            require(key in q_by_key, f"missing q record for recurrence {key}")
            q_record = q_by_key[key]
            t1 = complex_decimal(q_record["t1"])
            t2 = complex_decimal(q_record["t2"])
            t3 = complex_decimal(q_record["t3"])
            t4 = complex_decimal(q_record["t4"])
            t5 = complex_decimal(q_record["t5"])
            qq = complex_decimal(q_record["qq"])
            reconstructed = cadd(cconj(cadd(t1, t2, t4)), t3, t5)
            q_identity_error = max(
                abs(qq[0] - reconstructed[0]), abs(qq[1] - reconstructed[1])
            )
            q_scale = max(
                abs(part)
                for value in (t1, t2, t3, t4, t5, qq)
                for part in value
            )
            max_q_identity_error = max(max_q_identity_error, q_identity_error)
            max_q_identity_scaled_error = max(
                max_q_identity_scaled_error, q_identity_error / (Decimal(1) + q_scale)
            )

            recurrence_qq = complex_decimal(recurrence["qq"])
            max_logged_qq_disagreement = max(
                max_logged_qq_disagreement,
                abs(qq[0] - recurrence_qq[0]),
                abs(qq[1] - recurrence_qq[1]),
            )
            for field, destination in (
                ("local_defect", local_defect_norms),
                ("child_state_defect", child_state_defect_norms),
                ("parent_state_defect", parent_state_defect_norms),
            ):
                value = complex_decimal(recurrence[field])
                destination.append((value[0] * value[0] + value[1] * value[1]).sqrt())

    # The operands and result are independently rounded binary128 values.  A
    # scale-normalized guard tests their operation-rounding consistency without
    # pretending the decimal strings obey exact real arithmetic.
    require(
        max_q_identity_scaled_error <= Decimal("1e-32"),
        "logged q decomposition identity drift",
    )
    require(max_logged_qq_disagreement == 0, "q/recurrence qq record mismatch")

    return {
        "record_count": len(records),
        "chain_count": len(chains),
        "level_count": len(levels),
        "recurrence_count": len(recurrences),
        "q_term_count": len(q_terms),
        "exact_binary128_payloads_validated": True,
        "max_q_identity_abs_error": str(max_q_identity_error),
        "max_q_identity_scaled_error": str(max_q_identity_scaled_error),
        "max_local_defect_norm": str(max(local_defect_norms, default=Decimal(0))),
        "median_local_defect_norm": str(
            sorted(local_defect_norms)[len(local_defect_norms) // 2]
            if local_defect_norms
            else Decimal(0)
        ),
        "max_child_state_defect_norm": str(max(child_state_defect_norms, default=Decimal(0))),
        "max_parent_state_defect_norm": str(max(parent_state_defect_norms, default=Decimal(0))),
        "sha256": file_hash(path),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    parser.add_argument("--max-chains", type=int, default=64)
    args = parser.parse_args()
    require(1 <= args.max_chains <= 256, "max-chains must be in [1, 256]")

    require(file_hash(RESUMABLE_ROOT / "zeta14cubicmult_resumable.f90") == ACCEPTED_SOURCE_SHA256, "accepted source hash drift")
    run_command([sys.executable, str(SOURCE_CHECKER)])
    image_id = inspect_image_id()
    binary, build_result = build(args.build_root)

    disabled_root = FIXTURE_ROOT / "disabled"
    enabled_root = FIXTURE_ROOT / "enabled"
    reset_run_root(disabled_root)
    reset_run_root(enabled_root)
    disabled_run = invoke(binary, disabled_root, image_id, telemetry=False, max_chains=args.max_chains)
    enabled_run = invoke(binary, enabled_root, image_id, telemetry=True, max_chains=args.max_chains)

    disabled_journal = load_jsonl(disabled_root / "checkpoint.jsonl", decimal=True)
    enabled_journal = load_jsonl(enabled_root / "checkpoint.jsonl", decimal=True)
    reference_journal = load_jsonl(REFERENCE_JOURNAL, decimal=True)
    require(disabled_journal == enabled_journal, "enabled telemetry changed full-precision checkpoint state")
    require(
        without_run_id(disabled_journal) == without_run_id(reference_journal),
        "telemetry derivative changed accepted checkpoint state beyond provenance",
    )

    disabled_values = parse_grand_totals(disabled_root / "run.output.txt")
    enabled_values = parse_grand_totals(enabled_root / "run.output.txt")
    reference_values = parse_grand_totals(REFERENCE_OUTPUT)
    upstream_values = parse_grand_totals(UPSTREAM_OUTPUT)
    require(disabled_values == enabled_values, "enabled telemetry changed final displayed values")
    require(disabled_values == reference_values, "telemetry derivative changed accepted displayed values")
    upstream_max_error = max(abs(left - right) for left, right in zip(disabled_values, upstream_values))

    telemetry_result = validate_telemetry(
        enabled_root / "chain_telemetry.jsonl", args.max_chains
    )
    result = {
        "kind": "hardy_chain_telemetry_fixture",
        "image": IMAGE,
        "image_id": image_id,
        "mpi_ranks": 1,
        "cpu_limit": 1,
        "accepted_source_sha256": ACCEPTED_SOURCE_SHA256,
        "binary_path": str(binary),
        "build": build_result,
        "disabled_run": disabled_run,
        "enabled_run": enabled_run,
        "enabled_vs_disabled_journal_exact": True,
        "accepted_reference_journal_exact_except_run_id": True,
        "enabled_vs_disabled_output_exact": True,
        "accepted_reference_output_exact": True,
        "upstream_output_max_abs_error": str(upstream_max_error),
        "telemetry": telemetry_result,
        "hashes": {
            "checkpoint_module": file_hash(RESUMABLE_ROOT / "rh_hardy_checkpoint.f90"),
            "telemetry_module": file_hash(TELEMETRY_ROOT / "rh_hardy_chain_telemetry.f90"),
            "telemetry_source": file_hash(TELEMETRY_ROOT / "zeta14cubicmult_telemetry.f90"),
            "input": file_hash(INPUT),
            "reference_journal": file_hash(REFERENCE_JOURNAL),
            "disabled_journal": file_hash(disabled_root / "checkpoint.jsonl"),
            "enabled_journal": file_hash(enabled_root / "checkpoint.jsonl"),
        },
    }
    FIXTURE_ROOT.mkdir(parents=True, exist_ok=True)
    result_path = FIXTURE_ROOT / "fixture_result.json"
    temporary = result_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, result_path)
    print(
        "validated Hardy chain telemetry fixture: "
        f"{telemetry_result['chain_count']} chains, "
        f"{telemetry_result['recurrence_count']} recurrence rows, exact evaluator state"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
