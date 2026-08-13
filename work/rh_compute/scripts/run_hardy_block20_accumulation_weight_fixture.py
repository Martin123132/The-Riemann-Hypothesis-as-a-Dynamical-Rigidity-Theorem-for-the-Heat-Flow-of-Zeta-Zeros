#!/usr/bin/env python3
"""Emit the exact binary128 block-20 outer accumulation weights."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import run_hardy_chain_telemetry_fixture as base


EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
ACCEPTED_SOURCE = EXTERNAL_ROOT / "resumable/zeta14cubicmult_resumable.f90"
CHECKPOINT_MODULE = EXTERNAL_ROOT / "resumable/rh_hardy_checkpoint.f90"
REFERENCE_ROOT = EXTERNAL_ROOT / "fixtures/chain_telemetry/t1e10_block20_full/enabled"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/accumulation_weights/t1e10_block20"
RUN_ROOT = FIXTURE_ROOT / "run"
GENERATED_SOURCE = FIXTURE_ROOT / "generated/zeta14cubicmult_accumulation_weights.f90"
WEIGHTS = RUN_ROOT / "block20_accumulation_weights.txt"
RESULT = FIXTURE_ROOT / "fixture_result.json"
CHECKER = SCRIPT_ROOT / "check_hardy_block20_accumulation_weight_fixture.py"
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\accumulation_weights_20260807"
)
BINARY_NAME = "zeta14cubicmult.accumulation_weights"
IMAGE = base.IMAGE


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def replace_once(source: str, old: str, new: str, label: str) -> str:
    require(source.count(old) == 1, f"{label} source anchor count drift")
    return source.replace(old, new, 1)


def generate_source() -> dict[str, str]:
    original = ACCEPTED_SOURCE.read_text(encoding="utf-8")
    generated = original
    generated = replace_once(
        generated,
        "call alphasum(MTM,rae,pcsum)",
        "call alphasum(iblock,jsum,MTM,rae,pcsum)",
        "alphasum call",
    )
    generated = replace_once(
        generated,
        "subroutine alphasum(MTM,rae,pcsum)",
        "subroutine alphasum(iblock,jsum,MTM,rae,pcsum)",
        "alphasum signature",
    )
    generated = replace_once(
        generated,
        "  integer, parameter :: dp = selected_real_kind(33)  \n"
        "  integer, parameter :: dp1 = selected_int_kind(16) \n\n"
        "  integer       :: K,ip,mc,i\n"
        "  integer(dp1)  :: MTM,numbercalc,nchalf,jj",
        "  integer, parameter :: dp = selected_real_kind(33)  \n"
        "  integer, parameter :: dp1 = selected_int_kind(16) \n"
        "  integer, parameter :: hex_kind = selected_int_kind(32)\n\n"
        "  integer       :: K,ip,mc,i,iblock,rh_weight_unit\n"
        "  integer(dp1)  :: jsum,MTM,numbercalc,nchalf,jj",
        "alphasum declarations",
    )
    anchor = (
        "  do i=1,nchalf\n"
        "   c1=(0.0,1.0)*wp(i)\n"
        "   c1=exp(c1)\n"
        "   c2=(0.0,1.0)*wm(i)\n"
        "   c2=exp(c2)\n"
        "   pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))\n"
        "  enddo \n"
        "  return"
    )
    observer = (
        "  do i=1,nchalf\n"
        "   c1=(0.0,1.0)*wp(i)\n"
        "   c1=exp(c1)\n"
        "   c2=(0.0,1.0)*wm(i)\n"
        "   c2=exp(c2)\n"
        "   pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))\n"
        "  enddo \n\n"
        "  ! RH_ACCUMULATION_WEIGHT_OBSERVER_BEGIN\n"
        "  if (iblock.eq.20) then\n"
        "     open(newunit=rh_weight_unit,file='block20_accumulation_weights.txt', &\n"
        "          status='unknown',position='append',action='write')\n"
        "     do i=1,nchalf\n"
        "        c1=exp((0.0,1.0)*wp(i))\n"
        "        c2=exp((0.0,1.0)*wm(i))\n"
        "        write(rh_weight_unit,'(I0,1X,I0,1X,8(Z32.32,1X))') jsum,i, &\n"
        "             transfer(rae,0_hex_kind),transfer(amp(i),0_hex_kind), &\n"
        "             transfer(wp(i),0_hex_kind),transfer(wm(i),0_hex_kind), &\n"
        "             transfer(real(c1,kind=dp),0_hex_kind), &\n"
        "             transfer(aimag(c1),0_hex_kind), &\n"
        "             transfer(real(c2,kind=dp),0_hex_kind), &\n"
        "             transfer(aimag(c2),0_hex_kind)\n"
        "     enddo\n"
        "     flush(rh_weight_unit)\n"
        "     close(rh_weight_unit)\n"
        "  endif\n"
        "  ! RH_ACCUMULATION_WEIGHT_OBSERVER_END\n"
        "  return"
    )
    generated = replace_once(generated, anchor, observer, "accumulation observer")
    GENERATED_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    GENERATED_SOURCE.write_text(generated, encoding="utf-8")
    require(ACCEPTED_SOURCE.read_text(encoding="utf-8") == original, "accepted source changed")
    return {
        "accepted_source_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "generated_source_sha256": hashlib.sha256(generated.encode()).hexdigest(),
    }


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
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
    require(
        completed.returncode == 0,
        f"command failed ({completed.returncode}): {' '.join(command)}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
    )
    return completed


def build(build_root: Path) -> tuple[Path, dict[str, object]]:
    build_root.mkdir(parents=True, exist_ok=True)
    binary = build_root / BINARY_NAME
    command = [
        "docker",
        "run",
        "--rm",
        "--cpus=1",
        "-v",
        f"{base.docker_mount(REPO_ROOT)}:/rh:ro",
        "-v",
        f"{base.docker_mount(build_root)}:/build",
        "-w",
        "/build",
        IMAGE,
        "bash",
        "-lc",
        (
            "nice -n 10 mpif90 -O3 -z noexecstack -fallow-argument-mismatch "
            "-ffree-line-length-none /rh/"
            f"{relative(CHECKPOINT_MODULE)} /rh/{relative(GENERATED_SOURCE)} "
            f"-lpari -o /build/{BINARY_NAME}"
        ),
    ]
    completed = run_command(command)
    (build_root / "build.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (build_root / "build.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    require(binary.is_file(), "accumulation-weight binary was not produced")
    return binary, {
        "elapsed_seconds": completed.elapsed_seconds,  # type: ignore[attr-defined]
        "binary_sha256": file_hash(binary),
    }


def reset_run_root() -> None:
    resolved = RUN_ROOT.resolve()
    require(resolved.is_relative_to(FIXTURE_ROOT.resolve()), "unsafe run-root reset")
    if resolved.exists():
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    shutil.copy2(REFERENCE_ROOT / "inputs3.nml", resolved / "inputs3.nml")
    (resolved / "zeta14v6resa").write_text("", encoding="ascii")


def invoke(binary: Path) -> dict[str, object]:
    run_id = hashlib.sha256(
        (file_hash(binary) + file_hash(GENERATED_SOURCE)).encode("ascii")
    ).hexdigest()
    environment = {
        "RH_HARDY_CHECKPOINT": "/run/checkpoint.dat",
        "RH_HARDY_JOURNAL": "/run/checkpoint.jsonl",
        "RH_HARDY_STOP_FILE": "/run/stop.request",
        "RH_HARDY_RESUME": "0",
        "RH_HARDY_RS_CHUNK_SIZE": "257",
        "RH_HARDY_RUN_ID": run_id,
        "OMPI_ALLOW_RUN_AS_ROOT": "1",
        "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
    }
    command = ["docker", "run", "--rm", "--cpus=1"]
    for key, value in environment.items():
        command.extend(["-e", f"{key}={value}"])
    command.extend(
        [
            "-v",
            f"{base.docker_mount(binary.parent)}:/build:ro",
            "-v",
            f"{base.docker_mount(RUN_ROOT)}:/run",
            "-w",
            "/run",
            IMAGE,
            "bash",
            "-lc",
            f"nice -n 10 mpirun --allow-run-as-root -np 1 /build/{binary.name}",
        ]
    )
    completed = run_command(command)
    (RUN_ROOT / "run.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (RUN_ROOT / "run.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    shutil.copy2(RUN_ROOT / "zeta14v6resa", RUN_ROOT / "run.output.txt")
    return {
        "elapsed_seconds": completed.elapsed_seconds,  # type: ignore[attr-defined]
        "run_id": run_id,
    }


def validate_weights() -> dict[str, object]:
    lines = WEIGHTS.read_text(encoding="ascii").splitlines()
    require(len(lines) == 212 * 15, "accumulation-weight row-count drift")
    pattern = re.compile(r"([0-9]+) ([0-9]+)(?: ([0-9A-F]{32})){8} ?")
    roster: list[tuple[int, int]] = []
    for line_number, line in enumerate(lines, 1):
        require(pattern.fullmatch(line) is not None, f"invalid weight row {line_number}")
        fields = line.split()
        roster.append((int(fields[0]), int(fields[1])))
    expected = [(sum_index, output_index) for sum_index in range(1, 213) for output_index in range(1, 16)]
    require(roster == expected, "accumulation-weight roster drift")
    return {
        "row_count": len(lines),
        "sum_count": 212,
        "output_count": 15,
        "roster_exact": True,
        "binary128_payload_count": len(lines) * 8,
        "sha256": file_hash(WEIGHTS),
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    args = parser.parse_args()

    source = generate_source()
    image_id = base.inspect_image_id()
    binary, build_record = build(args.build_root)
    reset_run_root()
    run_record = invoke(binary)
    weights = validate_weights()

    new_journal = base.load_jsonl(RUN_ROOT / "checkpoint.jsonl", decimal=True)
    reference_journal = base.load_jsonl(REFERENCE_ROOT / "checkpoint.jsonl", decimal=True)
    require(
        base.without_run_id(new_journal) == base.without_run_id(reference_journal),
        "weight observer changed evaluator checkpoint state",
    )
    new_values = base.parse_grand_totals(RUN_ROOT / "run.output.txt")
    reference_values = base.parse_grand_totals(REFERENCE_ROOT / "run.output.txt")
    require(new_values == reference_values, "weight observer changed displayed Hardy values")

    artifact = {
        "kind": "hardy_block20_accumulation_weight_fixture",
        "status": "exact_binary128_212_by_15_weight_fixture_with_evaluator_equivalence_validated",
        "block": 20,
        "source": source,
        "image_id": image_id,
        "build": build_record,
        "run": run_record,
        "weights": weights,
        "equivalence": {
            "checkpoint_records": len(new_journal),
            "checkpoint_state_equal_ignoring_run_id": True,
            "displayed_hardy_values_exact": True,
            "displayed_value_count": len(new_values),
            "reference_checkpoint_sha256": file_hash(REFERENCE_ROOT / "checkpoint.jsonl"),
            "reference_output_sha256": file_hash(REFERENCE_ROOT / "run.output.txt"),
        },
        "artifacts": {
            "generated_source": {"path": relative(GENERATED_SOURCE), "sha256": file_hash(GENERATED_SOURCE)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "checkpoint": {"path": relative(RUN_ROOT / "checkpoint.jsonl"), "sha256": file_hash(RUN_ROOT / "checkpoint.jsonl")},
            "output": {"path": relative(RUN_ROOT / "run.output.txt"), "sha256": file_hash(RUN_ROOT / "run.output.txt")},
            "runner": {"path": relative(Path(__file__)), "sha256": file_hash(Path(__file__))},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Finite source-equivalent observer fixture for the exact binary128 outer amplitude and phase "
            "factors of block 20. It proves no accuracy estimate, height-uniform transport theorem, outer "
            "Hardy remainder, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print("built block-20 accumulation-weight fixture: 3180 rows, evaluator state unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
