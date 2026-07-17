#!/usr/bin/env python3
"""Build resumable exact H0-H23 anchors for the order-eleven lower bridge."""

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
    "jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_step2.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_step2_cache.json"
)
DEFAULT_START_T = Fraction(1244)
DEFAULT_END_T = Fraction(5708)
DEFAULT_STEP_T = Fraction(2)
PRECISION_BITS = 896
MODE_BISECTIONS = 180
PANELS = 80
WINDOW_Y = 15
TAYLOR_ORDER = 30
MAX_MOMENT = 23
DEFAULT_WORKERS = max(1, min(4, os.cpu_count() or 1))
ROW_CONTRACT = "order11-lower-sparse-point-h0-h23-step2-p896-b180-n80-w15-t30-v1"
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache.py"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mode_bracket(target: Fraction) -> tuple[Fraction, Fraction]:
    left = Fraction(9, 10)
    right = Fraction(2)
    target_arb = arb_rational(target)
    for _ in range(MODE_BISECTIONS):
        midpoint = (left + right) / 2
        value = potential_jet_arb(arb_rational(midpoint), 1)[1]
        if bool(value < target_arb):
            left = midpoint
        elif bool(value > target_arb):
            right = midpoint
        else:
            raise RuntimeError(f"inconclusive lower-sparse mode bracket at {target}")
    if not bool(
        potential_jet_arb(arb_rational(left), 1)[1]
        < target_arb
        < potential_jet_arb(arb_rational(right), 1)[1]
    ):
        raise RuntimeError(f"invalid lower-sparse mode bracket at {target}")
    return left, right


def deterministic_tasks(
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
) -> list[tuple[int, Fraction]]:
    if not 0 < start_t <= end_t or step_t <= 0:
        raise ValueError("invalid lower-sparse point range")
    quotient = (end_t - start_t) / step_t
    if quotient.denominator != 1:
        raise ValueError("lower-sparse point range does not align with its step")
    return [
        (index, start_t + index * step_t)
        for index in range(quotient.numerator + 1)
    ]


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def exact_task(task: tuple[int, Fraction]) -> dict:
    index, target = task
    flint.ctx.prec = PRECISION_BITS
    mode_left, mode_right = mode_bracket(target)
    result = integrate_h_jet_taylor_quadrature(
        mode_left,
        mode_right,
        PANELS,
        window_y=WINDOW_Y,
        taylor_order=TAYLOR_ORDER,
        max_moment=MAX_MOMENT,
    )
    if not result.get("passed"):
        return {
            "contract_id": ROW_CONTRACT,
            "index": index,
            "target_t": str(target),
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
            "contract_id": ROW_CONTRACT,
            "index": index,
            "target_t": str(target),
            "passed": False,
            "failure": "target-not-in-mode-bracket",
        }
    return {
        "kind": "order11_lower_sparse_point_h0_h23_jet",
        "contract_id": ROW_CONTRACT,
        "index": index,
        "target_t": str(target),
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


def load_cache(path: Path, expected: list[tuple[int, Fraction]]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(expected):
        raise RuntimeError("lower-sparse point cache has too many rows")
    derivative_keys = {str(order) for order in range(MAX_MOMENT + 1)}
    for record, task in zip(records, expected):
        index, target = task
        if (
            record.get("kind") != "order11_lower_sparse_point_h0_h23_jet"
            or record.get("contract_id") != ROW_CONTRACT
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != derivative_keys
        ):
            raise RuntimeError(f"lower-sparse point cache mismatch at row {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction]],
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
    executor = None
    if workers > 1:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
    else:
        initialize_worker()
    try:
        with path.open("a", encoding="utf-8") as handle:
            completed = 0
            cursor = 0
            while cursor < len(remaining):
                batch = remaining[cursor : cursor + workers]
                if executor is None:
                    batch_records = [exact_task(task) for task in batch]
                else:
                    batch_records = list(
                        executor.map(exact_task, batch, chunksize=1)
                    )
                for record in batch_records:
                    if record.get("passed") is not True:
                        raise RuntimeError(
                            f"lower-sparse row {record.get('index')} failed: "
                            f"{record.get('failure')}"
                        )
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    records.append(record)
                    completed += 1
                    if completed % 5 == 0:
                        handle.flush()
                    if completed % 20 == 0:
                        elapsed = perf_counter() - started
                        rate = completed / elapsed if elapsed else 0.0
                        eta = (
                            (len(remaining) - completed) / rate
                            if rate
                            else math.inf
                        )
                        print(
                            "order-eleven lower sparse H0-H23 rows: "
                            f"{len(records)}/{stop} "
                            f"({elapsed:.1f}s; ETA {eta / 3600:.2f}h)",
                            flush=True,
                        )
                cursor += len(batch)
                if (
                    runtime_seconds is not None
                    and perf_counter() - started >= runtime_seconds
                ):
                    handle.flush()
                    print(
                        "order-eleven lower sparse runtime checkpoint: "
                        f"{len(records)}/{stop} rows; "
                        "parking after the completed worker batch",
                        flush=True,
                    )
                    break
            handle.flush()
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def write_manifest(
    path: Path,
    cache_path: Path,
    records: list[dict],
    *,
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
) -> dict:
    manifest = {
        "kind": "jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache",
        "date": "2026-07-17",
        "status": "rigorous sparse exact H0-H23 cache; computational input only",
        "proof_boundary": (
            "These rows enclose exact H jets on a step-two lattice. They do not "
            "alone certify Taylor propagation or order-eleven curvature."
        ),
        "parameters": {
            "start_t": str(start_t),
            "end_t": str(end_t),
            "step_t": str(step_t),
            "precision_bits": PRECISION_BITS,
            "mode_bisections": MODE_BISECTIONS,
            "panels": PANELS,
            "window_y": WINDOW_Y,
            "taylor_order": TAYLOR_ORDER,
            "max_moment": MAX_MOMENT,
            "row_contract": ROW_CONTRACT,
        },
        "cache": {
            "path": cache_path.resolve().relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [0, MAX_MOMENT],
        },
        "diagnostics": {
            "maximum_panel_error_upper": max(
                (record["maximum_panel_error_upper"] for record in records), key=float
            ),
            "maximum_tail_moment_upper": max(
                (record["maximum_tail_moment_upper"] for record in records), key=float
            ),
        },
        "generator": GENERATOR_PATH,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--start-t", type=Fraction, default=DEFAULT_START_T)
    parser.add_argument("--end-t", type=Fraction, default=DEFAULT_END_T)
    parser.add_argument("--step-t", type=Fraction, default=DEFAULT_STEP_T)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-points", type=int)
    parser.add_argument("--runtime-seconds", type=float)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()
    if args.runtime_seconds is not None and args.runtime_seconds <= 0:
        parser.error("--runtime-seconds must be positive")
    if args.max_points is not None and args.max_points < 0:
        parser.error("--max-points cannot be negative")
    expected = deterministic_tasks(args.start_t, args.end_t, args.step_t)
    records = build_cache(
        args.cache,
        expected,
        workers=max(1, min(4, args.workers)),
        overwrite=args.overwrite,
        max_points=args.max_points,
        runtime_seconds=args.runtime_seconds,
    )
    print(f"order-eleven lower sparse point rows: {len(records)}/{len(expected)}")
    if args.cache_only:
        return 0
    if len(records) != len(expected):
        if args.runtime_seconds is not None and args.max_points is None:
            print("runtime-limited lower-sparse checkpoint is valid and resumable")
            return 0
        raise RuntimeError("complete the lower-sparse cache before writing its manifest")
    manifest = write_manifest(
        args.manifest,
        args.cache,
        records,
        start_t=args.start_t,
        end_t=args.end_t,
        step_t=args.step_t,
    )
    print(f"wrote lower-sparse manifest: {manifest['cache']['row_count']} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
