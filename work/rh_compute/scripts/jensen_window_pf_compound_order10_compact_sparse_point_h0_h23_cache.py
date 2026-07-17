#!/usr/bin/env python3
"""Build sparse exact H0-H23 jets for the order-ten compact bridge."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order10_compact_point_h0_h7_cache as profiles  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_jet_taylor_quadrature,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_step8.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache.json"
)
START_T = Fraction(5692)
END_T = Fraction(38020)
STEP_T = Fraction(8)
MEDIUM_LAST_T = Fraction(29996)
FAR_FIRST_T = Fraction(30004)
MAX_MOMENT = 23
WINDOW_Y = profiles.WINDOW_Y
TAYLOR_ORDER = profiles.TAYLOR_ORDER
DEFAULT_WORKERS = max(1, min(4, os.cpu_count() or 1))
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache.py"
)


def profile_for_target(target: Fraction) -> str:
    if START_T <= target <= MEDIUM_LAST_T:
        return "medium"
    if FAR_FIRST_T <= target <= END_T:
        return "far"
    raise ValueError(f"sparse target is outside the profile partition: {target}")


def row_contract(profile_name: str) -> str:
    profile = profiles.PROFILE_SPECS[profile_name]
    return (
        "order10-sparse-point-h0-h23-step8-"
        f"p{profile['precision_bits']}-b{profile['mode_bisections']}-"
        f"n{profile['panels']}-w{WINDOW_Y}-t{TAYLOR_ORDER}-v1"
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple[int, Fraction, str]]:
    quotient = (END_T - START_T) / STEP_T
    if quotient.denominator != 1:
        raise RuntimeError("sparse exact domain does not align with step eight")
    return [
        (index, START_T + index * STEP_T, profile_for_target(START_T + index * STEP_T))
        for index in range(quotient.numerator + 1)
    ]


def initialize_worker() -> None:
    flint.ctx.prec = profiles.PROFILE_SPECS["far"]["precision_bits"]


def exact_task(task: tuple[int, Fraction, str]) -> dict:
    index, target, profile_name = task
    profile = profiles.PROFILE_SPECS[profile_name]
    flint.ctx.prec = profile["precision_bits"]
    mode_left, mode_right = profiles.mode_bracket(target, profile_name)
    result = integrate_h_jet_taylor_quadrature(
        mode_left,
        mode_right,
        profile["panels"],
        window_y=WINDOW_Y,
        taylor_order=TAYLOR_ORDER,
        max_moment=MAX_MOMENT,
    )
    if not result.get("passed"):
        return {
            "contract_id": row_contract(profile_name),
            "index": index,
            "target_t": str(target),
            "profile": profile_name,
            "passed": False,
            "failure": result.get("failure"),
        }
    mode_ball = (
        (arb_rational(mode_left) + arb_rational(mode_right)) / 2
        + flint.arb(
            0,
            (arb_rational(mode_right) - arb_rational(mode_left)) / 2,
        )
    )
    target_ball = potential_jet_arb(mode_ball, 1)[1]
    if not bool(target_ball.contains(arb_rational(target))):
        return {
            "contract_id": row_contract(profile_name),
            "index": index,
            "target_t": str(target),
            "profile": profile_name,
            "passed": False,
            "failure": "target-not-in-mode-bracket",
        }
    return {
        "kind": "order10_compact_sparse_point_h0_h23_jet",
        "contract_id": row_contract(profile_name),
        "index": index,
        "target_t": str(target),
        "profile": profile_name,
        "mode_left": str(mode_left),
        "mode_right": str(mode_right),
        "target_t_ball": target_ball.str(50).replace("e", "E"),
        "passed": True,
        "h_derivatives": {
            str(order): compact.interval_text(result["h_derivatives"][order])
            for order in range(MAX_MOMENT + 1)
        },
        "maximum_panel_error_upper": result["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
    }


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
                    f"invalid sparse exact JSONL row {line_number}"
                ) from exc
    if len(records) > len(expected):
        raise RuntimeError("sparse exact cache has too many rows")
    derivative_keys = {str(order) for order in range(MAX_MOMENT + 1)}
    for record, task in zip(records, expected):
        index, target, profile_name = task
        if (
            record.get("kind") != "order10_compact_sparse_point_h0_h23_jet"
            or record.get("contract_id") != row_contract(profile_name)
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("profile") != profile_name
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != derivative_keys
        ):
            raise RuntimeError(f"sparse exact cache mismatch at row {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction, str]],
    *,
    workers: int,
    overwrite: bool,
    max_points: int | None,
    runtime_seconds: float | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, expected)
    stop = len(expected) if max_points is None else min(len(expected), max_points)
    remaining = expected[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker()
        iterator = map(exact_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(exact_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"sparse exact row {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 5 == 0:
                    handle.flush()
                if completed % 20 == 0:
                    elapsed = perf_counter() - started
                    rate = completed / elapsed if elapsed else 0.0
                    eta = (len(remaining) - completed) / rate if rate else math.inf
                    print(
                        "order-ten sparse exact H0-H23 rows: "
                        f"{len(records)}/{stop} ({elapsed:.1f}s; ETA {eta / 3600:.2f}h)",
                        flush=True,
                    )
                elapsed = perf_counter() - started
                if runtime_seconds is not None and elapsed >= runtime_seconds:
                    handle.flush()
                    print(
                        "order-ten sparse exact H0-H23 runtime checkpoint: "
                        f"{len(records)}/{stop} rows after {elapsed:.1f}s; "
                        "finishing in-flight workers and parking the cache",
                        flush=True,
                    )
                    break
            handle.flush()
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def write_manifest(path: Path, cache_path: Path, records: list[dict]) -> dict:
    profile_counts = {
        profile_name: sum(record["profile"] == profile_name for record in records)
        for profile_name in ("medium", "far")
    }
    manifest = {
        "kind": "jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache",
        "date": "2026-07-16",
        "status": (
            "rigorous sparse exact H0-H23 jet cache; computational input "
            "for order-24 Taylor propagation"
        ),
        "proof_boundary": (
            "This cache encloses exact first-summand H jets on an integer "
            "step-eight lattice. It does not itself certify the propagated "
            "integer grid or a curvature theorem."
        ),
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_t": str(STEP_T),
            "max_moment": MAX_MOMENT,
            "window_y": WINDOW_Y,
            "taylor_order": TAYLOR_ORDER,
            "profile_partition": {
                "medium": [str(START_T), str(MEDIUM_LAST_T)],
                "far": [str(FAR_FIRST_T), str(END_T)],
            },
            "profiles": {
                profile_name: {
                    **profiles.PROFILE_SPECS[profile_name],
                    "row_contract": row_contract(profile_name),
                }
                for profile_name in ("medium", "far")
            },
        },
        "cache": {
            "path": cache_path.resolve().relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "profile_row_counts": profile_counts,
            "all_rows_passed": True,
            "h_derivative_orders": [0, MAX_MOMENT],
        },
        "diagnostics": {
            "maximum_panel_error_upper": max(
                (record["maximum_panel_error_upper"] for record in records),
                key=float,
            ),
            "maximum_tail_moment_upper": max(
                (record["maximum_tail_moment_upper"] for record in records),
                key=float,
            ),
        },
        "generator": GENERATOR_PATH,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-points", type=int)
    parser.add_argument(
        "--runtime-seconds",
        type=float,
        help=(
            "Gracefully checkpoint after this many seconds: finish the current "
            "row, flush the JSONL cache, cancel queued futures, and wait for "
            "in-flight workers. Use 12600 for the workspace 3.5-hour soft stop."
        ),
    )
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()
    if args.runtime_seconds is not None and args.runtime_seconds <= 0:
        parser.error("--runtime-seconds must be positive")

    expected = deterministic_tasks()
    records = build_cache(
        args.cache,
        expected,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_points=args.max_points,
        runtime_seconds=args.runtime_seconds,
    )
    print(f"order-ten sparse exact H0-H23 cache rows: {len(records)}/{len(expected)}")
    if args.cache_only:
        return 0
    if len(records) != len(expected):
        if args.runtime_seconds is not None and args.max_points is None:
            print(
                "runtime-limited cache checkpoint is valid and resumable; "
                "rerun the same command to continue"
            )
            return 0
        raise RuntimeError("complete the sparse exact cache before writing its manifest")
    manifest = write_manifest(args.manifest, args.cache, records)
    print(
        "wrote order-ten sparse exact H0-H23 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
