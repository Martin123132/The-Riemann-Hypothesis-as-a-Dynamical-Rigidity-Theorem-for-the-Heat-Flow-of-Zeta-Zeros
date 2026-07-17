#!/usr/bin/env python3
"""Build reusable exact-point H0-H8 jets for order-nine continuation."""

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
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.json"
)
DEFAULT_START_T = Fraction(1243)
DEFAULT_END_T = Fraction(5707)
DEFAULT_STEP_T = Fraction(1, 2)
PRECISION_BITS = 512
MODE_BISECTIONS = 120
PANELS = 64
WINDOW_Y = 15
TAYLOR_ORDER = 30
MAX_MOMENT = 8
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


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
            raise RuntimeError(f"inconclusive point mode bracket at t={target}")
    if not bool(
        potential_jet_arb(arb_rational(left), 1)[1]
        < target_arb
        < potential_jet_arb(arb_rational(right), 1)[1]
    ):
        raise RuntimeError(f"invalid point mode bracket at t={target}")
    return left, right


def deterministic_tasks(
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
) -> list[tuple]:
    flint.ctx.prec = PRECISION_BITS
    if not 0 < start_t <= end_t or step_t <= 0:
        raise ValueError("invalid point H cache range")
    quotient = (end_t - start_t) / step_t
    if quotient.denominator != 1:
        raise ValueError("point H cache range must align with its step")
    targets = [
        start_t + index * step_t
        for index in range(quotient.numerator + 1)
    ]
    brackets = [mode_bracket(target) for target in targets]
    return [
        (index, target, brackets[index][0], brackets[index][1])
        for index, target in enumerate(targets)
    ]


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def point_task(task: tuple) -> dict:
    index, target, mode_left, mode_right = task
    flint.ctx.prec = PRECISION_BITS
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
            "index": index,
            "target_t": str(target),
            "mode_left": str(mode_left),
            "mode_right": str(mode_right),
            "passed": False,
            "failure": result.get("failure"),
        }
    target_ball = potential_jet_arb(
        (arb_rational(mode_left) + arb_rational(mode_right)) / 2
        + flint.arb(
            0,
            (arb_rational(mode_right) - arb_rational(mode_left)) / 2,
        ),
        1,
    )[1]
    if not bool(target_ball.contains(arb_rational(target))):
        return {
            "index": index,
            "target_t": str(target),
            "passed": False,
            "failure": "target-not-in-mode-bracket",
        }
    return {
        "kind": "order9_shifted_point_h0_h8_jet",
        "index": index,
        "target_t": str(target),
        "mode_left": str(mode_left),
        "mode_right": str(mode_right),
        "target_t_ball": target_ball.str(40).replace("e", "E"),
        "passed": True,
        "h_derivatives": {
            str(order): compact.interval_text(result["h_derivatives"][order])
            for order in range(MAX_MOMENT + 1)
        },
        "maximum_panel_error_upper": result["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
    }


def load_cache(path: Path, expected: list[tuple]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(expected):
        raise RuntimeError("point H cache has too many rows")
    for record, task in zip(records, expected):
        index, target, mode_left, mode_right = task
        if (
            record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("mode_left") != str(mode_left)
            or record.get("mode_right") != str(mode_right)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {}))
            != {str(order) for order in range(MAX_MOMENT + 1)}
        ):
            raise RuntimeError(f"point H cache mismatch at row {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple],
    *,
    workers: int,
    overwrite: bool,
    max_points: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, expected)
    stop = len(expected) if max_points is None else min(len(expected), max_points)
    remaining = expected[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    if workers == 1:
        initialize_worker()
        iterator = map(point_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(point_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"point H row {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 50 == 0:
                    handle.flush()
                    print(
                        "order-nine shifted point H0-H8 rows: "
                        f"{len(records)}/{stop} ({perf_counter() - start:.1f}s)"
                    )
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
    maximum_panel_error = max(
        (record["maximum_panel_error_upper"] for record in records),
        key=float,
    )
    maximum_tail_moment = max(
        (record["maximum_tail_moment_upper"] for record in records),
        key=float,
    )
    manifest = {
        "kind": "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache",
        "date": "2026-07-14",
        "status": "rigorous reusable exact-point H0-H8 jet cache; computational input only",
        "proof_boundary": (
            "This cache encloses H Taylor jets at exact half-grid targets. "
            "It does not itself prove any stable-coordinate or curvature bound."
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
        },
        "cache": {
            "path": cache_path.resolve().relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [0, MAX_MOMENT],
        },
        "diagnostics": {
            "maximum_panel_error_upper": maximum_panel_error,
            "maximum_tail_moment_upper": maximum_tail_moment,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.py"
        ),
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
    parser.add_argument("--start-t", type=Fraction, default=DEFAULT_START_T)
    parser.add_argument("--end-t", type=Fraction, default=DEFAULT_END_T)
    parser.add_argument("--step-t", type=Fraction, default=DEFAULT_STEP_T)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-points", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks(args.start_t, args.end_t, args.step_t)
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_points=args.max_points,
    )
    print(f"order-nine point H0-H8 cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(records) != len(tasks):
        raise RuntimeError("complete the point H cache before writing its manifest")
    manifest = write_manifest(
        args.manifest,
        args.cache,
        records,
        start_t=args.start_t,
        end_t=args.end_t,
        step_t=args.step_t,
    )
    print(
        "wrote order-nine point H0-H8 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
