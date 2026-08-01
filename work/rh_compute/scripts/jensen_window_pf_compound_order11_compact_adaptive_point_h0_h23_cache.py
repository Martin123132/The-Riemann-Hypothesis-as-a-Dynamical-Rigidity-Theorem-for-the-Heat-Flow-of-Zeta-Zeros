#!/usr/bin/env python3
"""Build the adaptive exact H0-H23 source for the full compact bridge."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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

import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_step4.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.json"
)
START_T = Fraction(5692)
END_T = Fraction(38028)
STEP_T = Fraction(4)
LOWER_LAST_T = Fraction(25016)
UPPER_FIRST_T = Fraction(25020)
LOWER_PROFILE_NAME = "order11_compact_lower_p896_n96"
UPPER_PROFILE_NAME = "order11_compact_upper_p896_n108"
PROFILES = {
    LOWER_PROFILE_NAME: {
        "precision_bits": 896,
        "mode_bisections": 200,
        "panels": 96,
        "row_contract": "order11-compact-adaptive-h0-h23-p896-b200-n96-w15-t30-v1",
    },
    UPPER_PROFILE_NAME: {
        "precision_bits": 896,
        "mode_bisections": 200,
        "panels": 108,
        "row_contract": "order11-compact-adaptive-h0-h23-p896-b200-n108-w15-t30-v1",
    },
}
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def profile_for_target(target: Fraction) -> str:
    if START_T <= target <= LOWER_LAST_T:
        return LOWER_PROFILE_NAME
    if UPPER_FIRST_T <= target <= END_T:
        return UPPER_PROFILE_NAME
    raise ValueError(f"compact target is outside the profile partition: {target}")


def deterministic_tasks() -> list[tuple[int, Fraction, str]]:
    quotient = (END_T - START_T) / STEP_T
    if quotient.denominator != 1:
        raise RuntimeError("compact exact source domain does not align")
    tasks = []
    for index in range(quotient.numerator + 1):
        target = START_T + index * STEP_T
        tasks.append((index, target, profile_for_target(target)))
    if tasks[-1][1] != END_T:
        raise RuntimeError("compact exact source task partition does not close")
    return tasks


def source_contract() -> dict:
    payload = {
        "geometry": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_t": str(STEP_T),
            "lower_last_t": str(LOWER_LAST_T),
            "upper_first_t": str(UPPER_FIRST_T),
        },
        "profiles": PROFILES,
        "derivative_orders": [0, source.MAX_MOMENT],
        "window_y": source.WINDOW_Y,
        "taylor_order": source.TAYLOR_ORDER,
    }
    contract_id = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return {"id": contract_id, **payload}


SOURCE_CONTRACT = source_contract()


def upstream_row_contract(profile_name: str) -> str:
    source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    return source.row_contract(profile_name)


def contract_provenance() -> dict:
    return {
        "record_contract_id_semantics": (
            "Historical identifier emitted by the reused exact H0-H23 quadrature "
            "kernel; its order10 and step8 tokens name implementation provenance, "
            "not the adaptive source geometry."
        ),
        "authoritative_outer_binding": (
            "Every row source_contract_id equals source_contract.id; that outer "
            "contract fixes the step-four geometry and adaptive profile partition."
        ),
        "profiles": {
            profile_name: {
                "adaptive_profile_contract_id": profile["row_contract"],
                "record_contract_id": upstream_row_contract(profile_name),
            }
            for profile_name, profile in PROFILES.items()
        },
    }


def exact_task(task: tuple[int, Fraction, str]) -> dict:
    index, target, profile_name = task
    source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    row = source.exact_task((index, target, profile_name))
    return {
        **row,
        "kind": "order11_compact_adaptive_point_h0_h23_jet",
        "source_contract_id": SOURCE_CONTRACT["id"],
    }


def validate_record(record: dict, task: tuple[int, Fraction, str]) -> None:
    index, target, profile_name = task
    source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    if (
        record.get("kind") != "order11_compact_adaptive_point_h0_h23_jet"
        or record.get("source_contract_id") != SOURCE_CONTRACT["id"]
        or record.get("index") != index
        or record.get("target_t") != str(target)
        or record.get("profile") != profile_name
        or record.get("contract_id") != upstream_row_contract(profile_name)
        or record.get("passed") is not True
        or set(record.get("h_derivatives", {}))
        != {str(order) for order in range(source.MAX_MOMENT + 1)}
    ):
        raise RuntimeError(f"invalid adaptive compact source row {index}")


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
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid adaptive compact JSONL row {line_number}"
                ) from exc
    if len(records) > len(expected):
        raise RuntimeError("adaptive compact source cache has too many rows")
    for record, task in zip(records, expected):
        validate_record(record, task)
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
    with ProcessPoolExecutor(max_workers=workers) as pool, path.open(
        "a", encoding="utf-8", newline="\n"
    ) as handle:
        while len(records) < stop:
            if perf_counter() - started >= runtime_seconds:
                break
            batch = expected[len(records) : min(stop, len(records) + workers)]
            batch_rows = list(pool.map(exact_task, batch))
            for record, task in zip(batch_rows, batch):
                validate_record(record, task)
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
            handle.flush()
            os.fsync(handle.fileno())
    return records, len(records) - initial, perf_counter() - started


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
        record["profile"] == LOWER_PROFILE_NAME for record in records
    )
    upper_rows = sum(
        record["profile"] == UPPER_PROFILE_NAME for record in records
    )
    manifest = {
        "kind": "jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache",
        "date": "2026-07-18",
        "status": (
            "rigorous complete adaptive exact H0-H23 compact source"
            if complete
            else "rigorous resumable adaptive exact H0-H23 compact source prefix"
        ),
        "proof_boundary": (
            "This cache supplies exact first-summand H0-H23 jets. It does not "
            "itself certify any compact curvature block."
        ),
        "contract_provenance": contract_provenance(),
        "source_contract": SOURCE_CONTRACT,
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_t": str(STEP_T),
            "required_rows": len(expected),
            "profiles": PROFILES,
        },
        "cache": {
            "path": relative(cache_path),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "required_row_count": len(expected),
            "all_rows_passed": all(record.get("passed") is True for record in records),
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
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py"
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
    args = parser.parse_args()
    expected = deterministic_tasks()
    records, added, elapsed = build_cache(
        args.cache,
        expected,
        workers=max(1, args.workers),
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
        "adaptive compact exact source: "
        f"{cache['row_count']}/{cache['required_row_count']} rows, "
        f"added={latest['rows_added']}, elapsed={latest['elapsed_seconds']}s, "
        f"next={latest['next_target_t']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
