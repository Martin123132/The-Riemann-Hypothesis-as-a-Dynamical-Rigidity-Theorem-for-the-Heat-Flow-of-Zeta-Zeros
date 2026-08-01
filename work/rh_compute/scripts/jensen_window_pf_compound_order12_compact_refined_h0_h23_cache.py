#!/usr/bin/env python3
"""Build window-16 H0-H23 jets on the inherited step-four compact lattice."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import ctypes
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache as base  # noqa: E402


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_compact_refined_h0_h23_step4.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_compact_refined_h0_h23_cache.json"
)
DEFAULT_PILOT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_compact_refined_h0_h23_pilot.json"
)

START_T = base.START_T
END_T = base.END_T
STEP_T = base.STEP_T
PROFILES = base.PROFILES
MAX_MOMENT = 23
WINDOW_Y = 16
TAYLOR_ORDER = 34
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def record_sha256(record: dict) -> str:
    payload = json.dumps(
        record,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def profile_for_target(target: Fraction) -> str:
    return base.profile_for_target(target)


def deterministic_tasks() -> list[tuple[int, Fraction, str]]:
    tasks = base.deterministic_tasks()
    if len(tasks) != 8085 or tasks[0][1] != START_T or tasks[-1][1] != END_T:
        raise RuntimeError("order-twelve refined task geometry changed")
    return tasks


def source_contract() -> dict:
    payload = {
        "geometry": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_t": str(STEP_T),
            "role": "window-16 replacement on the inherited step-four lattice",
            "lower_last_t": str(base.LOWER_LAST_T),
            "upper_first_t": str(base.UPPER_FIRST_T),
        },
        "profiles": PROFILES,
        "derivative_orders": [0, MAX_MOMENT],
        "window_y": WINDOW_Y,
        "taylor_order": TAYLOR_ORDER,
        "inherited_geometry_contract": base.SOURCE_CONTRACT,
    }
    identifier = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return {"id": identifier, **payload}


SOURCE_CONTRACT = source_contract()


def _set_below_normal_priority() -> None:
    if os.name == "nt":
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
        kernel32.SetPriorityClass.restype = ctypes.c_int
        handle = kernel32.GetCurrentProcess()
        if not kernel32.SetPriorityClass(handle, 0x00004000):
            raise ctypes.WinError(ctypes.get_last_error())


def initialize_worker() -> None:
    _set_below_normal_priority()
    base.source.MAX_MOMENT = MAX_MOMENT
    base.source.WINDOW_Y = WINDOW_Y
    base.source.TAYLOR_ORDER = TAYLOR_ORDER


def exact_task(task: tuple[int, Fraction, str]) -> dict:
    index, target, profile_name = task
    base.source.MAX_MOMENT = MAX_MOMENT
    base.source.WINDOW_Y = WINDOW_Y
    base.source.TAYLOR_ORDER = TAYLOR_ORDER
    base.source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    row = base.source.exact_task((index, target, profile_name))
    return {
        **row,
        "kind": "order12_compact_refined_h0_h23_jet",
        "source_contract_id": SOURCE_CONTRACT["id"],
    }


def validate_record(record: dict, task: tuple[int, Fraction, str]) -> None:
    index, target, profile_name = task
    base.source.MAX_MOMENT = MAX_MOMENT
    base.source.WINDOW_Y = WINDOW_Y
    base.source.TAYLOR_ORDER = TAYLOR_ORDER
    base.source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    if (
        record.get("kind") != "order12_compact_refined_h0_h23_jet"
        or record.get("source_contract_id") != SOURCE_CONTRACT["id"]
        or record.get("index") != index
        or record.get("target_t") != str(target)
        or record.get("profile") != profile_name
        or record.get("contract_id") != base.source.row_contract(profile_name)
        or record.get("passed") is not True
        or set(record.get("h_derivatives", {}))
        != {str(order) for order in range(MAX_MOMENT + 1)}
    ):
        raise RuntimeError(f"invalid order-twelve refined row {index} at t={target}")


def load_cache(
    path: Path,
    expected: list[tuple[int, Fraction, str]],
) -> list[dict]:
    if not path.exists():
        return []
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid refined JSONL row {line_number}"
                ) from exc
            if len(records) >= len(expected):
                raise RuntimeError("refined cache has too many rows")
            validate_record(record, expected[len(records)])
            records.append(record)
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction, str]],
    *,
    workers: int,
    runtime_seconds: float,
    max_points: int | None,
) -> tuple[list[dict], int, float]:
    if runtime_seconds <= 0:
        raise ValueError("runtime_seconds must be positive")
    records = load_cache(path, expected)
    stop = len(expected) if max_points is None else min(len(expected), max_points)
    if len(records) >= stop:
        return records, 0, 0.0
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    initial = len(records)
    executor = None
    try:
        if workers > 1:
            executor = ProcessPoolExecutor(
                max_workers=workers,
                initializer=initialize_worker,
            )
        else:
            initialize_worker()
        with path.open("a", encoding="utf-8", newline="\n") as handle:
            while len(records) < stop:
                if perf_counter() - started >= runtime_seconds:
                    break
                batch = expected[len(records) : min(stop, len(records) + workers)]
                completed = (
                    [exact_task(task) for task in batch]
                    if executor is None
                    else list(executor.map(exact_task, batch, chunksize=1))
                )
                for record, task in zip(completed, batch):
                    validate_record(record, task)
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    records.append(record)
                handle.flush()
                os.fsync(handle.fileno())
                if len(records) % (workers * 5) == 0 or len(records) == stop:
                    elapsed = perf_counter() - started
                    rate = (len(records) - initial) / elapsed if elapsed else 0.0
                    print(
                        "order-twelve refined H0-H23 rows: "
                        f"{len(records)}/{stop} ({elapsed:.1f}s; {rate:.3f} rows/s)",
                        flush=True,
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records, len(records) - initial, perf_counter() - started


def build_pilot(anchors: list[Fraction], *, workers: int) -> dict:
    expected_by_target = {target: task for task in deterministic_tasks() for target in (task[1],)}
    unique = sorted(set(anchors))
    if not unique or len(unique) != len(anchors):
        raise ValueError("pilot anchors must be nonempty and distinct")
    try:
        tasks = [expected_by_target[target] for target in unique]
    except KeyError as exc:
        raise ValueError(
            f"pilot target is not on the inherited step-four lattice: {exc.args[0]}"
        ) from exc
    if workers > 1:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(tasks)),
            initializer=initialize_worker,
        ) as executor:
            rows = list(executor.map(exact_task, tasks, chunksize=1))
    else:
        initialize_worker()
        rows = [exact_task(task) for task in tasks]
    for row, task in zip(rows, tasks):
        validate_record(row, task)
    return {
        "kind": "order12_compact_refined_h0_h23_pilot",
        "date": "2026-07-22",
        "status": "rigorous named window-16 H0-H23 source rows",
        "proof_boundary": (
            "These rows test the proposed refined step-four source. They do not "
            "certify contiguous compact curvature coverage or prove RH."
        ),
        "source_contract": SOURCE_CONTRACT,
        "anchors": [str(target) for target in unique],
        "rows": rows,
        "row_sha256": [record_sha256(row) for row in rows],
        "summary": {"rows": len(rows), "failed_rows": 0},
        "generator": relative(Path(__file__)),
    }


def write_manifest(
    path: Path,
    cache_path: Path,
    records: list[dict],
    expected: list[tuple[int, Fraction, str]],
    *,
    added: int,
    elapsed_seconds: float,
) -> dict:
    complete = len(records) == len(expected)
    next_task = expected[len(records)] if not complete else None
    lower_rows = sum(
        record["profile"] == base.LOWER_PROFILE_NAME for record in records
    )
    upper_rows = sum(
        record["profile"] == base.UPPER_PROFILE_NAME for record in records
    )
    manifest = {
        "kind": "jensen_window_pf_compound_order12_compact_refined_h0_h23_cache",
        "date": "2026-07-22",
        "status": (
            "rigorous complete refined H0-H23 compact source"
            if complete
            else "rigorous resumable refined H0-H23 compact source prefix"
        ),
        "proof_boundary": (
            "This cache supplies tighter exact jets for rigorous order-24 Taylor "
            "transport on the inherited step-four lattice."
        ),
        "source_contract": SOURCE_CONTRACT,
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_t": str(STEP_T),
            "max_moment": MAX_MOMENT,
            "window_y": WINDOW_Y,
            "taylor_order": TAYLOR_ORDER,
            "required_rows": len(expected),
            "profiles": PROFILES,
        },
        "cache": {
            "path": relative(cache_path),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "required_row_count": len(expected),
            "all_rows_passed": all(row.get("passed") is True for row in records),
            "complete": complete,
            "lower_profile_rows": lower_rows,
            "upper_profile_rows": upper_rows,
        },
        "latest_run": {
            "rows_added": added,
            "elapsed_seconds": f"{elapsed_seconds:.1f}",
            "next_index": next_task[0] if next_task else None,
            "next_target_t": str(next_task[1]) if next_task else None,
            "next_profile": next_task[2] if next_task else None,
        },
        "generator": relative(Path(__file__)),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order12_compact_refined_h0_h23_cache.py"
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--runtime-seconds", type=float, default=3600.0)
    parser.add_argument("--max-points", type=int)
    parser.add_argument("--pilot-anchor", action="append", type=Fraction, default=[])
    parser.add_argument("--pilot-output", type=Path, default=DEFAULT_PILOT)
    args = parser.parse_args()
    workers = max(1, min(4, args.workers))
    if args.max_points is not None and args.max_points < 0:
        parser.error("--max-points cannot be negative")
    if args.pilot_anchor:
        artifact = build_pilot(args.pilot_anchor, workers=workers)
        args.pilot_output.parent.mkdir(parents=True, exist_ok=True)
        args.pilot_output.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote order-twelve refined pilot: {len(artifact['rows'])} rows")
        return 0
    expected = deterministic_tasks()
    records, added, elapsed = build_cache(
        args.cache,
        expected,
        workers=workers,
        runtime_seconds=args.runtime_seconds,
        max_points=args.max_points,
    )
    manifest = write_manifest(
        args.manifest,
        args.cache,
        records,
        expected,
        added=added,
        elapsed_seconds=elapsed,
    )
    cache = manifest["cache"]
    latest = manifest["latest_run"]
    print(
        "order-twelve refined compact source: "
        f"{cache['row_count']}/{cache['required_row_count']} rows, "
        f"added={latest['rows_added']}, elapsed={latest['elapsed_seconds']}s, "
        f"next={latest['next_target_t']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
