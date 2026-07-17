#!/usr/bin/env python3
"""Build a t-uniform H2-H22 cache for the order-seven compact bridge."""

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
    "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_tenth_tiles.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache.json"
)
START_T = Fraction(315)
END_T = Fraction(1005)
TILE_WIDTH_T = Fraction(1, 10)
MODE_BISECTIONS = 60
MAX_MOMENT = 22
PRECISION_BITS = 256
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))
CACHE_TILE_LABEL = "tenth"
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache.py"
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
            raise RuntimeError(f"inconclusive mode bracket at t={target}")
    if not bool(
        potential_jet_arb(arb_rational(left), 1)[1]
        < target_arb
        < potential_jet_arb(arb_rational(right), 1)[1]
    ):
        raise RuntimeError(f"invalid mode bracket at t={target}")
    return left, right


def deterministic_tasks() -> list[tuple]:
    flint.ctx.prec = PRECISION_BITS
    tile_count = int((END_T - START_T) / TILE_WIDTH_T)
    targets = [START_T + index * TILE_WIDTH_T for index in range(tile_count + 1)]
    brackets = [mode_bracket(target) for target in targets]
    return [
        (
            index,
            targets[index],
            targets[index + 1],
            brackets[index][0],
            brackets[index + 1][1],
        )
        for index in range(tile_count)
    ]


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def tile_task(task: tuple) -> dict:
    index, t_left, t_right, mode_left, mode_right = task
    flint.ctx.prec = PRECISION_BITS
    result = integrate_h_derivatives(
        mode_left,
        mode_right,
        compact.PANELS,
        window_y=compact.WINDOW_Y,
        eighth_envelope_bound=compact.EIGHTH_ENVELOPE,
        max_moment=MAX_MOMENT,
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
        "kind": f"order7_shifted_jet_compact_h2_h22_{CACHE_TILE_LABEL}_tile",
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


def load_cache(path: Path, expected: list[tuple]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(expected):
        raise RuntimeError("compact H cache has too many rows")
    for record, task in zip(records, expected):
        index, t_left, t_right, mode_left, mode_right = task
        if (
            record.get("index") != index
            or record.get("target_t_left") != str(t_left)
            or record.get("target_t_right") != str(t_right)
            or record.get("mode_left") != str(mode_left)
            or record.get("mode_right") != str(mode_right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"compact H cache mismatch at tile {index}")
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
                        f"compact H tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 500 == 0:
                    handle.flush()
                    print(
                        f"order-seven compact H2-H22 tiles: "
                        f"{len(records)}/{stop} ({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def write_manifest(path: Path, cache_path: Path, records: list[dict]) -> dict:
    manifest = {
        "kind": "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache",
        "date": "2026-07-13",
        "status": "rigorous t-uniform H2-H22 interval cache; computational input, not a curvature theorem",
        "proof_boundary": (
            f"This cache encloses first-summand H derivatives on "
            f"{START_T}<=t<={END_T}. "
            "It does not itself prove any stable-gap or curvature inequality."
        ),
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "tile_width_t": str(TILE_WIDTH_T),
            "mode_bisections": MODE_BISECTIONS,
            "max_moment": MAX_MOMENT,
            "precision_bits": PRECISION_BITS,
            "panels": compact.PANELS,
            "window_y": compact.WINDOW_Y,
            "eighth_envelope": str(compact.EIGHTH_ENVELOPE),
        },
        "cache": {
            "path": cache_path.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [2, MAX_MOMENT],
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
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    expected = deterministic_tasks()
    records = build_cache(
        args.cache,
        expected,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_tiles=args.max_tiles,
    )
    print(f"order-seven compact H2-H22 cache rows: {len(records)}/{len(expected)}")
    if args.cache_only:
        return 0
    if len(records) != len(expected):
        raise RuntimeError("complete the compact H cache before writing its manifest")
    manifest = write_manifest(args.manifest, args.cache, records)
    print(
        "wrote order-seven compact H2-H22 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
