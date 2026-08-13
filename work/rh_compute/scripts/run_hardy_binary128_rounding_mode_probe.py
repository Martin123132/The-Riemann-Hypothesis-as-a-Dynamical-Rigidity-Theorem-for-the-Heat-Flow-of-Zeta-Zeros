#!/usr/bin/env python3
"""Compile and run the pinned binary128 rounding-mode probe."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
FIXTURE_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/binary128_rounding_mode_probe"
PROBE_SOURCE = FIXTURE_ROOT / "probe.f90"
PROBE_STDOUT = FIXTURE_ROOT / "probe.stdout.txt"
PROBE_STDERR = FIXTURE_ROOT / "probe.stderr.txt"
RESULT = FIXTURE_ROOT / "probe_result.json"
ACCEPTED_SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_hardy_binary128_rounding_mode_probe.py"
DEFAULT_BUILD_ROOT = Path(r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\binary128_rounding_mode_probe_20260809")
IMAGE = "rh-hardy-fastcode:bookworm"
FLAGS = ("-O3", "-ffree-line-length-none")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def run(command: list[str]) -> tuple[subprocess.CompletedProcess[str], float]:
    started = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace")
    elapsed = time.perf_counter() - started
    require(
        completed.returncode == 0,
        f"command exited {completed.returncode}: {' '.join(command)}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
    )
    return completed, elapsed


def parse_probe(stdout: str) -> tuple[str, dict[str, str]]:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    require(lines and "GNU Fortran" in lines[0], "missing compiler banner")
    values: dict[str, str] = {}
    for line in lines[1:]:
        require("=" in line, f"malformed probe line: {line}")
        key, value = line.split("=", 1)
        require(key not in values, f"duplicate probe key {key}")
        values[key] = value
    return lines[0], values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    args = parser.parse_args()
    for path in (PROBE_SOURCE, ACCEPTED_SOURCE, CHECKER):
        require(path.is_file(), f"missing rounding-mode probe dependency: {path}")
    priority = set_low_priority()
    args.build_root.mkdir(parents=True, exist_ok=True)

    image_result, _ = run(["docker", "image", "inspect", IMAGE, "--format", "{{.Id}}"])
    image_id = image_result.stdout.strip()
    command = [
        "docker",
        "run",
        "--rm",
        "--cpus=1",
        "--memory=512m",
        "--pids-limit=128",
        "-v",
        f"{PROBE_SOURCE.parent.resolve()}:/probe:ro",
        "-v",
        f"{args.build_root.resolve()}:/build",
        "-w",
        "/build",
        IMAGE,
        "bash",
        "-lc",
        "gfortran --version | head -n 1; "
        "gfortran -O3 -ffree-line-length-none /probe/probe.f90 -o rounding_mode_probe; "
        "nice -n 10 ./rounding_mode_probe",
    ]
    completed, elapsed = run(command)
    PROBE_STDOUT.write_text(completed.stdout, encoding="utf-8")
    PROBE_STDERR.write_text(completed.stderr, encoding="utf-8")
    compiler, values = parse_probe(completed.stdout)

    accepted_text = ACCEPTED_SOURCE.read_text(encoding="utf-8").lower()
    mode_set_tokens = ("ieee_set_rounding_mode", "fesetround")
    mode_set_counts = {token: accepted_text.count(token) for token in mode_set_tokens}
    passed = (
        values.get("radix") == "2"
        and values.get("digits") == "113"
        and values.get("min_exponent") == "-16381"
        and values.get("max_exponent") == "16384"
        and values.get("ieee_datatype") == "T"
        and values.get("ieee_nearest_supported") == "T"
        and values.get("rounding_mode") == "nearest"
        and all(count == 0 for count in mode_set_counts.values())
    )
    artifact = {
        "kind": "hardy_binary128_rounding_mode_probe",
        "status": "pinned_source_environment_starts_binary128_in_round_to_nearest_mode" if passed else "binary128_rounding_mode_probe_failed",
        "passed": passed,
        "probe_values": values,
        "source_rounding_mode_mutator_counts": mode_set_counts,
        "build": {
            "image": IMAGE,
            "image_id": image_id,
            "compiler": compiler,
            "compiler_flags": list(FLAGS),
            "cpu_cap": 1,
            "memory_cap": "512m",
            "container_nice": 10,
            "controller_priority": priority,
            "elapsed_seconds": elapsed,
        },
        "artifacts": {
            "probe_source": {"path": relative(PROBE_SOURCE), "sha256": file_hash(PROBE_SOURCE)},
            "stdout": {"path": relative(PROBE_STDOUT), "sha256": file_hash(PROBE_STDOUT)},
            "stderr": {"path": relative(PROBE_STDERR), "sha256": file_hash(PROBE_STDERR)},
            "accepted_source": {"path": relative(ACCEPTED_SOURCE), "sha256": file_hash(ACCEPTED_SOURCE)},
            "runner": {"path": relative(Path(__file__).resolve()), "sha256": file_hash(Path(__file__).resolve())},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": "Pins the process-start rounding mode and binary format in the accepted Docker/compiler environment. It does not by itself certify every source operation, intrinsic special-function error, other build environment, or RH.",
    }
    require(passed, "binary128 rounding-mode probe did not pass")
    tmp = RESULT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(RESULT)
    print(f"built binary128 rounding-mode probe: {values['rounding_mode']}, digits={values['digits']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
