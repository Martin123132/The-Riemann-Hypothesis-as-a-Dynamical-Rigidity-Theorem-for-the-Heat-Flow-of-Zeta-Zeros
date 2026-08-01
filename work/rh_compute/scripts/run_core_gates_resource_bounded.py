#!/usr/bin/env python3
"""Run core proof gates serially with resumable progress and CPU parking."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import psutil

import check_core_proof_programme_gates as core


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PROGRESS = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "core_resource_bounded_progress.jsonl"
)
DEFAULT_SUMMARY = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "core_resource_bounded_summary.json"
)
POLL_SECONDS = 0.25
CPU_SAMPLE_SECONDS = 5.0


def set_below_normal() -> None:
    try:
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def append_jsonl(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def read_records(path: Path) -> list[dict]:
    if not path.exists():
        return []
    records: list[dict] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"invalid progress JSONL line {line_number}: {exc}"
            ) from exc
        if not isinstance(record, dict):
            raise RuntimeError(
                f"non-object progress JSONL line {line_number}"
            )
        records.append(record)
    return records


def output_tail(text: str, lines: int = 12) -> str:
    return "\n".join(text.splitlines()[-lines:])


def baseline_samples(count: int = 6) -> list[float]:
    return [
        round(psutil.cpu_percent(interval=1.0), 1)
        for _ in range(count)
    ]


def command_for(spec: core.GateSpec) -> list[str]:
    return [sys.executable, *spec.command]


def run_one(
    index: int,
    spec: core.GateSpec,
    timeout_seconds: int,
    cpu_threshold: float,
) -> tuple[dict, bool]:
    started = time.perf_counter()
    creationflags = getattr(
        subprocess,
        "BELOW_NORMAL_PRIORITY_CLASS",
        0,
    )
    with tempfile.TemporaryFile(
        mode="w+",
        encoding="utf-8",
    ) as stdout_handle, tempfile.TemporaryFile(
        mode="w+",
        encoding="utf-8",
    ) as stderr_handle:
        process = subprocess.Popen(
            command_for(spec),
            cwd=REPO_ROOT,
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
            creationflags=creationflags,
        )
        cpu_samples: list[float] = []
        consecutive_high = 0
        park_requested = False
        last_cpu_sample = time.perf_counter()
        psutil.cpu_percent(interval=None)

        while process.poll() is None:
            elapsed = time.perf_counter() - started
            if elapsed > timeout_seconds:
                process.terminate()
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=15)
                stdout_handle.seek(0)
                stderr_handle.seek(0)
                stdout = stdout_handle.read()
                stderr = stderr_handle.read()
                return (
                    {
                        "type": "gate",
                        "index": index,
                        "name": spec.name,
                        "category": spec.category,
                        "command": " ".join(command_for(spec)),
                        "ok": False,
                        "timed_out": True,
                        "elapsed_seconds": round(
                            time.perf_counter() - started,
                            3,
                        ),
                        "cpu_samples": cpu_samples,
                        "stdout_tail": output_tail(stdout),
                        "stderr_tail": output_tail(stderr),
                    },
                    True,
                )
            now = time.perf_counter()
            if now - last_cpu_sample >= CPU_SAMPLE_SECONDS:
                cpu = round(psutil.cpu_percent(interval=None), 1)
                cpu_samples.append(cpu)
                consecutive_high = (
                    consecutive_high + 1
                    if cpu > cpu_threshold
                    else 0
                )
                if consecutive_high >= 2:
                    park_requested = True
                last_cpu_sample = now
            time.sleep(POLL_SECONDS)

        process.wait()
        stdout_handle.seek(0)
        stderr_handle.seek(0)
        stdout = stdout_handle.read()
        stderr = stderr_handle.read()
    combined = stdout + "\n" + stderr
    missing = [
        needle for needle in spec.expected if needle not in combined
    ]
    ok = process.returncode == 0 and not missing
    record = {
        "type": "gate",
        "index": index,
        "name": spec.name,
        "category": spec.category,
        "command": " ".join(command_for(spec)),
        "returncode": process.returncode,
        "ok": ok,
        "timed_out": False,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "cpu_samples": cpu_samples,
        "missing_expected": missing,
        "stdout_tail": output_tail(stdout),
        "stderr_tail": output_tail(stderr),
    }
    return record, park_requested


def build_summary(
    mode: str,
    specs: list[core.GateSpec],
    records: list[dict],
    baseline: list[float],
    status: str,
    progress_path: Path,
) -> dict:
    gate_records = [
        record for record in records if record.get("type") == "gate"
    ]
    latest_by_index = {
        int(record["index"]): record for record in gate_records
    }
    completed = [
        record
        for record in latest_by_index.values()
        if record.get("ok")
    ]
    failed = [
        record
        for record in latest_by_index.values()
        if not record.get("ok")
    ]
    return {
        "kind": "rh_core_resource_bounded_summary",
        "status": status,
        "mode": mode,
        "core_runner_sha256": digest(Path(core.__file__)),
        "progress_sha256": digest(progress_path),
        "total_selected_gates": len(specs),
        "completed_gates": len(completed),
        "failed_gates": len(failed),
        "remaining_gates": len(specs) - len(completed),
        "baseline_cpu_samples": baseline,
        "maximum_observed_cpu": max(
            (
                sample
                for record in gate_records
                for sample in record.get("cpu_samples", [])
            ),
            default=None,
        ),
        "failed_names": [record["name"] for record in failed],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("day", "night"),
        default="day",
    )
    parser.add_argument("--skip-slow", action="store_true")
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--max-new-gates", type=int)
    parser.add_argument("--progress", type=Path, default=DEFAULT_PROGRESS)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def main() -> int:
    set_below_normal()
    args = parse_args()
    specs = [
        spec
        for spec in core.GATES
        if not (args.skip_slow and spec.slow)
    ]
    core_hash = digest(Path(core.__file__))
    records = read_records(args.progress)
    headers = [
        record for record in records if record.get("type") == "session"
    ]
    if headers:
        header = headers[0]
        if header.get("core_runner_sha256") != core_hash:
            raise RuntimeError(
                "progress belongs to a different core registry hash"
            )
        if header.get("selected_gate_count") != len(specs):
            raise RuntimeError(
                "progress belongs to a different selected gate set"
            )
    else:
        append_jsonl(
            args.progress,
            {
                "type": "session",
                "core_runner_sha256": core_hash,
                "selected_gate_count": len(specs),
                "skip_slow": args.skip_slow,
            },
        )
        records = read_records(args.progress)

    successful = {
        int(record["index"])
        for record in records
        if record.get("type") == "gate" and record.get("ok")
    }
    baseline = baseline_samples()
    average = sum(baseline) / len(baseline)
    maximum = max(baseline)
    if args.mode == "day":
        threshold = 75.0
        if average > 60.0:
            summary = build_summary(
                args.mode,
                specs,
                records,
                baseline,
                "deferred_busy_baseline",
                args.progress,
            )
            write_json(args.summary, summary)
            print(
                "deferred RH core replay: daytime baseline "
                f"average={average:.1f}, max={maximum:.1f}"
            )
            return 75
    else:
        threshold = 85.0
        if average + 20.0 > 75.0:
            summary = build_summary(
                args.mode,
                specs,
                records,
                baseline,
                "deferred_busy_baseline",
                args.progress,
            )
            write_json(args.summary, summary)
            print(
                "deferred RH core replay: night baseline plus one-worker "
                f"allowance exceeds 75, average={average:.1f}, "
                f"max={maximum:.1f}"
            )
            return 75

    launched = 0
    for index, spec in enumerate(specs, start=1):
        if index in successful:
            continue
        if (
            args.max_new_gates is not None
            and launched >= args.max_new_gates
        ):
            break
        record, park_requested = run_one(
            index,
            spec,
            args.timeout,
            threshold,
        )
        append_jsonl(args.progress, record)
        records.append(record)
        launched += 1
        status = "OK" if record["ok"] else "FAIL"
        print(
            f"{status} resource-bounded core gate "
            f"{index}/{len(specs)}: {spec.name} "
            f"({record['elapsed_seconds']}s)",
            flush=True,
        )
        if not record["ok"]:
            summary = build_summary(
                args.mode,
                specs,
                records,
                baseline,
                "failed",
                args.progress,
            )
            write_json(args.summary, summary)
            return 1
        if park_requested:
            summary = build_summary(
                args.mode,
                specs,
                records,
                baseline,
                "resource_parked",
                args.progress,
            )
            write_json(args.summary, summary)
            print(
                "parked RH core replay after the completed atomic gate "
                f"{index}/{len(specs)}",
                flush=True,
            )
            return 75

    records = read_records(args.progress)
    completed = {
        int(record["index"])
        for record in records
        if record.get("type") == "gate" and record.get("ok")
    }
    status = (
        "complete"
        if len(completed) == len(specs)
        else "bounded_batch_complete"
    )
    summary = build_summary(
        args.mode,
        specs,
        records,
        baseline,
        status,
        args.progress,
    )
    write_json(args.summary, summary)
    print(
        f"resource-bounded RH core replay status={status}: "
        f"{len(completed)}/{len(specs)} gates",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
