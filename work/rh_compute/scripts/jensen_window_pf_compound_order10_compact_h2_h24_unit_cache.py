#!/usr/bin/env python3
"""Build unit-tile H2-H24 enclosures for the order-ten compact bridge."""

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
    integrate_h_derivatives,
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
    "jensen_window_pf_compound_order10_compact_h2_h24_unit_tiles.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_h2_h24_unit_cache.json"
)
DEFAULT_START_T = Fraction(5692)
DEFAULT_END_T = Fraction(38028)
DEFAULT_TILE_WIDTH_T = Fraction(1)
MODE_LEFT = Fraction(3, 2)
MODE_RIGHT = Fraction(201, 100)
PRECISION_BITS = 256
MODE_BISECTIONS = 120
PANELS = 200
WINDOW_Y = 15
MAX_MOMENT = 24
DEFAULT_WORKERS = max(1, min(6, os.cpu_count() or 1))
ROW_CONTRACT = "order10-h2-h24-p256-b120-n200-w15-unit-v1"
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order10_compact_h2_h24_unit_cache.py"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mode_bracket(target: Fraction) -> tuple[Fraction, Fraction]:
    """Enclose the unique saddle mode whose potential derivative is target."""
    left = MODE_LEFT
    right = MODE_RIGHT
    target_arb = arb_rational(target)
    left_value = potential_jet_arb(arb_rational(left), 1)[1]
    right_value = potential_jet_arb(arb_rational(right), 1)[1]
    if not bool(left_value < target_arb < right_value):
        raise RuntimeError(f"initial tile mode bracket misses t={target}")
    for _ in range(MODE_BISECTIONS):
        midpoint = (left + right) / 2
        value = potential_jet_arb(arb_rational(midpoint), 1)[1]
        if bool(value < target_arb):
            left = midpoint
        elif bool(value > target_arb):
            right = midpoint
        else:
            raise RuntimeError(f"inconclusive tile mode bracket at t={target}")
    if not bool(
        potential_jet_arb(arb_rational(left), 1)[1]
        < target_arb
        < potential_jet_arb(arb_rational(right), 1)[1]
    ):
        raise RuntimeError(f"invalid tile mode bracket at t={target}")
    return left, right


def deterministic_tasks(
    start_t: Fraction,
    end_t: Fraction,
    tile_width_t: Fraction,
) -> list[tuple[int, Fraction, Fraction]]:
    if not 0 < start_t < end_t or tile_width_t <= 0:
        raise ValueError("invalid compact H cache range")
    quotient = (end_t - start_t) / tile_width_t
    if quotient.denominator != 1:
        raise ValueError("compact H cache range must align with its tile width")
    return [
        (
            index,
            start_t + index * tile_width_t,
            start_t + (index + 1) * tile_width_t,
        )
        for index in range(quotient.numerator)
    ]


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def tile_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, t_left, t_right = task
    flint.ctx.prec = PRECISION_BITS
    mode_left = mode_bracket(t_left)[0]
    mode_right = mode_bracket(t_right)[1]
    result = integrate_h_derivatives(
        mode_left,
        mode_right,
        PANELS,
        window_y=WINDOW_Y,
        eighth_envelope_bound=compact.EIGHTH_ENVELOPE,
        max_moment=MAX_MOMENT,
    )
    if not result.get("passed"):
        return {
            "contract_id": ROW_CONTRACT,
            "index": index,
            "target_t_left": str(t_left),
            "target_t_right": str(t_right),
            "mode_left": str(mode_left),
            "mode_right": str(mode_right),
            "passed": False,
            "failure": result.get("failure"),
        }
    actual_t_left = potential_jet_arb(arb_rational(mode_left), 1)[1]
    actual_t_right = potential_jet_arb(arb_rational(mode_right), 1)[1]
    if not bool(actual_t_left < arb_rational(t_left)) or not bool(
        actual_t_right > arb_rational(t_right)
    ):
        return {
            "contract_id": ROW_CONTRACT,
            "index": index,
            "target_t_left": str(t_left),
            "target_t_right": str(t_right),
            "mode_left": str(mode_left),
            "mode_right": str(mode_right),
            "passed": False,
            "failure": "mode-bracket-coverage",
        }
    return {
        "kind": "order10_compact_h2_h24_unit_tile",
        "contract_id": ROW_CONTRACT,
        "index": index,
        "target_t_left": str(t_left),
        "target_t_right": str(t_right),
        "mode_left": str(mode_left),
        "mode_right": str(mode_right),
        "actual_t_left": compact.interval_text(actual_t_left),
        "actual_t_right": compact.interval_text(actual_t_right),
        "passed": True,
        "h_derivatives": {
            str(order): compact.interval_text(result["h_derivatives"][order])
            for order in range(2, MAX_MOMENT + 1)
        },
        "normalizer_lower": result["normalizer_lower"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
    }


def load_cache(
    path: Path,
    expected: list[tuple[int, Fraction, Fraction]],
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
                    f"invalid compact H JSONL row {line_number}"
                ) from exc
    if len(records) > len(expected):
        raise RuntimeError("compact H cache has too many rows")
    derivative_keys = {str(order) for order in range(2, MAX_MOMENT + 1)}
    for record, task in zip(records, expected):
        index, t_left, t_right = task
        if (
            record.get("kind") != "order10_compact_h2_h24_unit_tile"
            or record.get("contract_id") != ROW_CONTRACT
            or record.get("index") != index
            or record.get("target_t_left") != str(t_left)
            or record.get("target_t_right") != str(t_right)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != derivative_keys
        ):
            raise RuntimeError(f"compact H cache mismatch at tile {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    overwrite: bool,
    max_tiles: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, expected)
    stop = len(expected) if max_tiles is None else min(len(expected), max_tiles)
    remaining = expected[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker()
        iterator = map(tile_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(tile_task, remaining, chunksize=4)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"compact H tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 50 == 0:
                    handle.flush()
                if completed % 250 == 0:
                    elapsed = perf_counter() - started
                    rate = completed / elapsed if elapsed else 0.0
                    remaining_seconds = (
                        (len(remaining) - completed) / rate if rate else math.inf
                    )
                    print(
                        "order-ten compact H2-H24 unit tiles: "
                        f"{len(records)}/{stop} ({elapsed:.1f}s; "
                        f"ETA {remaining_seconds / 60:.1f}m)",
                        flush=True,
                    )
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
    tile_width_t: Fraction,
) -> dict:
    manifest = {
        "kind": "jensen_window_pf_compound_order10_compact_h2_h24_unit_cache",
        "date": "2026-07-16",
        "status": (
            "rigorous t-uniform H2-H24 interval cache; computational input, "
            "not a curvature theorem"
        ),
        "proof_boundary": (
            f"This cache encloses first-summand H derivatives on "
            f"{start_t}<=t<={end_t}. It does not itself prove a stable-gap "
            "or curvature inequality."
        ),
        "parameters": {
            "start_t": str(start_t),
            "end_t": str(end_t),
            "tile_width_t": str(tile_width_t),
            "initial_mode_bracket": [str(MODE_LEFT), str(MODE_RIGHT)],
            "mode_bisections": MODE_BISECTIONS,
            "max_moment": MAX_MOMENT,
            "precision_bits": PRECISION_BITS,
            "panels": PANELS,
            "window_y": WINDOW_Y,
            "eighth_envelope": str(compact.EIGHTH_ENVELOPE),
            "row_contract": ROW_CONTRACT,
        },
        "cache": {
            "path": cache_path.resolve().relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [2, MAX_MOMENT],
        },
        "diagnostics": {
            "maximum_tail_moment_upper": max(
                (record["maximum_tail_moment_upper"] for record in records),
                key=float,
            ),
            "maximum_simpson_error_upper": max(
                (record["maximum_simpson_error_upper"] for record in records),
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
    parser.add_argument("--start-t", type=Fraction, default=DEFAULT_START_T)
    parser.add_argument("--end-t", type=Fraction, default=DEFAULT_END_T)
    parser.add_argument(
        "--tile-width-t",
        type=Fraction,
        default=DEFAULT_TILE_WIDTH_T,
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks(args.start_t, args.end_t, args.tile_width_t)
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_tiles=args.max_tiles,
    )
    print(f"order-ten compact H2-H24 cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(records) != len(tasks):
        raise RuntimeError("complete the compact H cache before writing its manifest")
    manifest = write_manifest(
        args.manifest,
        args.cache,
        records,
        start_t=args.start_t,
        end_t=args.end_t,
        tile_width_t=args.tile_width_t,
    )
    print(
        "wrote order-ten compact H2-H24 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
