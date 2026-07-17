#!/usr/bin/env python3
"""Build an exact-t H-derivative cache for localized order-nine curvature."""

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
    "jensen_window_pf_compound_order9_shifted_jet_compact_h2_h18_half_tiles.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_jet_compact_h2_h18_cache.json"
)
DEFAULT_START_T = Fraction(1243)
DEFAULT_END_T = Fraction(10007)
DEFAULT_TILE_WIDTH_T = Fraction(1, 2)
MODE_BISECTIONS = 60
MAX_MOMENT = 18
PRECISION_BITS = 256
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
            raise RuntimeError(f"inconclusive mode bracket at t={target}")
    if not bool(
        potential_jet_arb(arb_rational(left), 1)[1]
        < target_arb
        < potential_jet_arb(arb_rational(right), 1)[1]
    ):
        raise RuntimeError(f"invalid mode bracket at t={target}")
    return left, right


def deterministic_tasks(
    start_t: Fraction,
    end_t: Fraction,
    tile_width_t: Fraction,
    max_moment: int = MAX_MOMENT,
) -> list[tuple]:
    flint.ctx.prec = PRECISION_BITS
    if not 0 < start_t < end_t or tile_width_t <= 0:
        raise ValueError("invalid exact-t cache range")
    quotient = (end_t - start_t) / tile_width_t
    if quotient.denominator != 1:
        raise ValueError("cache range must contain an integral number of tiles")
    targets = [
        start_t + index * tile_width_t
        for index in range(quotient.numerator + 1)
    ]
    brackets = [mode_bracket(target) for target in targets]
    return [
        (
            index,
            targets[index],
            targets[index + 1],
            brackets[index][0],
            brackets[index + 1][1],
            max_moment,
        )
        for index in range(quotient.numerator)
    ]


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def tile_task(task: tuple) -> dict:
    index, t_left, t_right, mode_left, mode_right, max_moment = task
    flint.ctx.prec = PRECISION_BITS
    result = integrate_h_derivatives(
        mode_left,
        mode_right,
        compact.PANELS,
        window_y=compact.WINDOW_Y,
        eighth_envelope_bound=compact.EIGHTH_ENVELOPE,
        max_moment=max_moment,
    )
    if not result.get("passed"):
        return {
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
            "index": index,
            "passed": False,
            "failure": "mode-bracket-coverage",
        }
    return {
        "kind": "order9_shifted_jet_compact_h2_h18_exact_t_tile",
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
            for order in range(2, max_moment + 1)
        },
        "normalizer_lower": result["normalizer_lower"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
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
        raise RuntimeError("order-nine exact-t H cache has too many rows")
    for record, task in zip(records, expected):
        index, t_left, t_right, mode_left, mode_right, max_moment = task
        if (
            record.get("index") != index
            or record.get("target_t_left") != str(t_left)
            or record.get("target_t_right") != str(t_right)
            or record.get("mode_left") != str(mode_left)
            or record.get("mode_right") != str(mode_right)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {}))
            != {str(order) for order in range(2, max_moment + 1)}
        ):
            raise RuntimeError(f"order-nine exact-t H cache mismatch at {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple],
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
    start = perf_counter()
    if workers == 1:
        iterator = map(tile_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(tile_task, remaining, chunksize=8)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"exact-t H tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 250 == 0:
                    handle.flush()
                    print(
                        "order-nine exact-t H2-H18 tiles: "
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
    tile_width_t: Fraction,
    max_moment: int,
) -> dict:
    manifest = {
        "kind": "jensen_window_pf_compound_order9_shifted_jet_compact_h_derivative_cache",
        "date": "2026-07-14",
        "status": (
            f"rigorous exact-t H2-H{max_moment} interval cache; "
            "computational input only"
        ),
        "proof_boundary": (
            "This cache encloses first-summand H derivatives on exact rational "
            "t tiles. It does not itself prove a stable-gap or curvature bound."
        ),
        "parameters": {
            "start_t": str(start_t),
            "end_t": str(end_t),
            "tile_width_t": str(tile_width_t),
            "mode_bisections": MODE_BISECTIONS,
            "max_moment": max_moment,
            "precision_bits": PRECISION_BITS,
            "panels": compact.PANELS,
            "window_y": compact.WINDOW_Y,
            "eighth_envelope": str(compact.EIGHTH_ENVELOPE),
        },
        "cache": {
            "path": cache_path.resolve().relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [2, max_moment],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_shifted_jet_compact_h2_h18_cache.py"
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
    parser.add_argument(
        "--tile-width-t", type=Fraction, default=DEFAULT_TILE_WIDTH_T
    )
    parser.add_argument("--max-moment", type=int, default=MAX_MOMENT)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    if args.max_moment < 2:
        raise ValueError("max moment must be at least two")
    tasks = deterministic_tasks(
        args.start_t,
        args.end_t,
        args.tile_width_t,
        args.max_moment,
    )
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_tiles=args.max_tiles,
    )
    print(
        f"order-nine exact-t H2-H{args.max_moment} cache rows: "
        f"{len(records)}/{len(tasks)}"
    )
    if args.cache_only:
        return 0
    if len(records) != len(tasks):
        raise RuntimeError("complete the exact-t H cache before writing its manifest")
    manifest = write_manifest(
        args.manifest,
        args.cache,
        records,
        start_t=args.start_t,
        end_t=args.end_t,
        tile_width_t=args.tile_width_t,
        max_moment=args.max_moment,
    )
    print(
        f"wrote order-nine exact-t H2-H{args.max_moment} cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
