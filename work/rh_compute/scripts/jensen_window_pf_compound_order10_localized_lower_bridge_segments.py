#!/usr/bin/env python3
"""Build adaptive order-ten localized lower-bridge segment certificates."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
from functools import lru_cache
import math
import mmap
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
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    localized_seventh_formula_continuation_row,
)


NEAR_H_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_lower_bridge_h2_h24_hundredth_tiles.jsonl"
)
NEAR_H_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_lower_bridge_h2_h24_hundredth_cache.json"
)
BROAD_H_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_lower_bridge_h2_h24_tenth_tiles.jsonl"
)
BROAD_H_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_lower_bridge_h2_h24_tenth_cache.json"
)
POINT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
)
POINT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.json"
)
POINT_EXTENSION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_point_h0_h8_half_grid_extension.jsonl"
)
POINT_EXTENSION_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_point_h0_h8_half_grid_extension_cache.json"
)
DEFAULT_SEGMENT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_localized_lower_bridge_segments.jsonl"
)
DEFAULT_RUN_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_localized_lower_bridge_run_contract.json"
)

START_T = Fraction(1251)
NEAR_END_T = Fraction(1500)
MIDDLE_END_T = Fraction(1800)
END_T = Fraction(5700)
CELL_WIDTH = Fraction(1, 2)
QUARTER_WIDTH = Fraction(1, 4)
ROOT_SEGMENT_WIDTH = Fraction(16)
POINT_ORDER = 7
REMAINDER_ORDER = 8
CURVATURE_CONSTANT = 5500

NEAR_H_START = Fraction(1243)
NEAR_H_END = Fraction(1509)
NEAR_H_STEP = Fraction(1, 100)
NEAR_H_ROWS = 26600
BROAD_H_START = Fraction(1491)
BROAD_H_END = Fraction(5709)
BROAD_H_STEP = Fraction(1, 10)
BROAD_H_ROWS = 42180
H_ORDERS = tuple(range(2, 25))

POINT_START = Fraction(1243)
POINT_END = Fraction(5707)
POINT_STEP = Fraction(1, 2)
POINT_ROWS = 8929
POINT_EXTENSION_START = Fraction(11415, 2)
POINT_EXTENSION_END = Fraction(5708)
POINT_EXTENSION_ROWS = 2
POINT_ORDERS = tuple(range(9))

DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))


_HANDLES = {}
_VIEWS = {}
_OFFSETS = {}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def _open_jsonl_mmap(path: Path):
    handle = path.open("rb")
    view = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
    offsets = []
    position = 0
    size = len(view)
    while position < size:
        offsets.append(position)
        newline = view.find(b"\n", position)
        position = size if newline < 0 else newline + 1
    return handle, view, offsets


def _jsonl_record(kind: str, index: int) -> dict:
    offsets = _OFFSETS[kind]
    view = _VIEWS[kind]
    if not 0 <= index < len(offsets):
        raise RuntimeError(f"{kind} JSONL row {index} is unavailable")
    end = offsets[index + 1] if index + 1 < len(offsets) else len(view)
    return json.loads(view[offsets[index] : end].strip())


def initialize_worker(
    near_h: str,
    broad_h: str,
    point: str,
    point_extension: str,
) -> None:
    flint.ctx.prec = PRECISION_BITS
    for kind, raw_path, expected in (
        ("near", near_h, NEAR_H_ROWS),
        ("broad", broad_h, BROAD_H_ROWS),
        ("point", point, POINT_ROWS),
        ("point_extension", point_extension, POINT_EXTENSION_ROWS),
    ):
        handle, view, offsets = _open_jsonl_mmap(Path(raw_path))
        if len(offsets) != expected:
            raise RuntimeError(
                f"{kind} source has {len(offsets)} rows, expected {expected}"
            )
        _HANDLES[kind] = handle
        _VIEWS[kind] = view
        _OFFSETS[kind] = offsets


def _aligned_index(value: Fraction, start: Fraction, step: Fraction) -> int:
    quotient = (value - start) / step
    if quotient.denominator != 1:
        raise RuntimeError(f"unaligned source target {value}")
    return quotient.numerator


@lru_cache(maxsize=5000)
def _load_h_row(kind: str, index: int) -> dict:
    if kind == "near":
        start, step, expected_rows = NEAR_H_START, NEAR_H_STEP, NEAR_H_ROWS
    elif kind == "broad":
        start, step, expected_rows = BROAD_H_START, BROAD_H_STEP, BROAD_H_ROWS
    else:
        raise RuntimeError(f"unknown H source kind: {kind}")
    if not 0 <= index < expected_rows:
        raise RuntimeError(f"{kind} H row {index} leaves its cache")
    record = _jsonl_record(kind, index)
    left = start + index * step
    right = left + step
    derivatives = record.get("h_derivatives", {})
    if (
        record.get("index") != index
        or record.get("target_t_left") != str(left)
        or record.get("target_t_right") != str(right)
        or record.get("passed") is not True
        or set(derivatives) != {str(order) for order in H_ORDERS}
    ):
        raise RuntimeError(f"invalid {kind} H source row {index}")
    return {
        "target_t_left": left,
        "target_t_right": right,
        "H": {
            order: compact.interval_from_text(derivatives[str(order)])
            for order in H_ORDERS
        },
    }


@lru_cache(maxsize=512)
def _load_point(target: Fraction) -> tuple[list, dict]:
    if target <= POINT_END:
        kind = "point"
        index = _aligned_index(target, POINT_START, POINT_STEP)
        expected_target = POINT_START + index * POINT_STEP
    else:
        kind = "point_extension"
        index = _aligned_index(target, POINT_EXTENSION_START, POINT_STEP)
        expected_target = POINT_EXTENSION_START + index * POINT_STEP
    record = _jsonl_record(kind, index)
    derivatives = record.get("h_derivatives", {})
    if (
        record.get("index") != index
        or record.get("target_t") != str(expected_target)
        or expected_target != target
        or record.get("passed") is not True
        or set(derivatives) != {str(order) for order in POINT_ORDERS}
    ):
        raise RuntimeError(f"invalid exact point source at t={target}")
    return (
        [
            compact.interval_from_text(derivatives[str(order)])
            / math.factorial(order)
            for order in POINT_ORDERS
        ],
        {
            "target_t": str(target),
            "mode_bracket": [record["mode_left"], record["mode_right"]],
            "maximum_panel_error_upper": record["maximum_panel_error_upper"],
            "maximum_tail_moment_upper": record[
                "maximum_tail_moment_upper"
            ],
            "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
        },
    )


def point_source_for(anchor: Fraction) -> dict:
    return {
        anchor + shift: _load_point(anchor + shift)
        for shift in range(-8, 9)
    }


def h_rows_for_cell(regime: str, left: Fraction, right: Fraction) -> list[dict]:
    if regime == "near":
        kind, start, step, row_count = (
            "near",
            NEAR_H_START,
            NEAR_H_STEP,
            NEAR_H_ROWS,
        )
    else:
        kind, start, step, row_count = (
            "broad",
            BROAD_H_START,
            BROAD_H_STEP,
            BROAD_H_ROWS,
        )
    support_left = left - 8
    support_right = right + 8
    first = _aligned_index(support_left, start, step)
    last = _aligned_index(support_right, start, step)
    if not 0 <= first < last <= row_count:
        raise RuntimeError(
            f"{regime} H collar {support_left}..{support_right} leaves its cache"
        )
    return [_load_h_row(kind, index) for index in range(first, last)]


def deterministic_segments() -> list[tuple[int, str, Fraction, Fraction]]:
    segments = []
    index = 0
    for regime, start, end in (
        ("near", START_T, NEAR_END_T),
        ("middle", NEAR_END_T, MIDDLE_END_T),
        ("far", MIDDLE_END_T, END_T),
    ):
        left = start
        while left < end:
            right = min(left + ROOT_SEGMENT_WIDTH, end)
            segments.append((index, regime, left, right))
            index += 1
            left = right
    return segments


def segment_task(task: tuple[int, str, Fraction, Fraction]) -> dict:
    index, regime, segment_left, segment_right = task
    flint.ctx.prec = PRECISION_BITS
    blocks = []
    cell_left = segment_left
    while cell_left < segment_right:
        cell_right = cell_left + CELL_WIDTH
        if cell_right > segment_right:
            raise RuntimeError(f"segment {index} cuts a half-unit cell")
        h_rows = h_rows_for_cell(regime, cell_left, cell_right)
        if regime in {"near", "middle"}:
            specs = (
                (cell_left, cell_left, cell_left + QUARTER_WIDTH),
                (cell_right, cell_left + QUARTER_WIDTH, cell_right),
            )
        else:
            specs = ((cell_left, cell_left, cell_right),)
        for expansion, left, right in specs:
            block = localized_seventh_formula_continuation_row(
                expansion,
                right,
                h_rows,
                point_order=POINT_ORDER,
                remainder_order=REMAINDER_ORDER,
                point_h_source=point_source_for(expansion),
                block_left=left,
            )
            blocks.append({**block, "regime": regime})
        cell_left = cell_right
    if not blocks or any(row.get("passed") is not True for row in blocks):
        raise RuntimeError(f"segment {index} did not produce only passing blocks")
    maximum = max(blocks, key=lambda row: float(row["scaled_curvature_upper"]))
    minimum = min(blocks, key=lambda row: float(row["curvature_margin_lower"]))
    return {
        "kind": "order10_localized_lower_bridge_segment",
        "index": index,
        "regime": regime,
        "segment_left": str(segment_left),
        "segment_right": str(segment_right),
        "blocks": blocks,
        "block_count": len(blocks),
        "largest_scaled_curvature_upper": maximum["scaled_curvature_upper"],
        "largest_scaled_curvature_anchor": maximum["anchor"],
        "smallest_margin_lower": minimum["curvature_margin_lower"],
        "smallest_margin_anchor": minimum["anchor"],
        "passed": True,
    }


def validate_source_manifest(
    manifest_path: Path,
    cache_path: Path,
    *,
    start: str,
    end: str,
    step: str,
    rows: int,
    max_moment: int,
) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    parameters = manifest.get("parameters", {})
    cache = manifest.get("cache", {})
    if (
        parameters.get("start_t") != start
        or parameters.get("end_t") != end
        or parameters.get("tile_width_t") != step
        or parameters.get("max_moment") != max_moment
        or cache.get("row_count") != rows
        or cache.get("all_rows_passed") is not True
        or cache.get("sha256") != sha256(cache_path)
    ):
        raise RuntimeError(f"source manifest changed: {manifest_path}")
    return manifest


def canonical_run_contract() -> dict:
    validate_source_manifest(
        NEAR_H_MANIFEST,
        NEAR_H_CACHE,
        start=str(NEAR_H_START),
        end=str(NEAR_H_END),
        step=str(NEAR_H_STEP),
        rows=NEAR_H_ROWS,
        max_moment=24,
    )
    validate_source_manifest(
        BROAD_H_MANIFEST,
        BROAD_H_CACHE,
        start=str(BROAD_H_START),
        end=str(BROAD_H_END),
        step=str(BROAD_H_STEP),
        rows=BROAD_H_ROWS,
        max_moment=24,
    )
    point_manifest = json.loads(POINT_MANIFEST.read_text(encoding="utf-8"))
    extension_manifest = json.loads(
        POINT_EXTENSION_MANIFEST.read_text(encoding="utf-8")
    )
    for manifest, cache_path, rows in (
        (point_manifest, POINT_CACHE, POINT_ROWS),
        (extension_manifest, POINT_EXTENSION, POINT_EXTENSION_ROWS),
    ):
        cache = manifest.get("cache", {})
        if (
            cache.get("row_count") != rows
            or cache.get("all_rows_passed") is not True
            or cache.get("sha256") != sha256(cache_path)
        ):
            raise RuntimeError(f"point source changed: {cache_path}")
    paths = (
        NEAR_H_CACHE,
        NEAR_H_MANIFEST,
        BROAD_H_CACHE,
        BROAD_H_MANIFEST,
        POINT_CACHE,
        POINT_MANIFEST,
        POINT_EXTENSION,
        POINT_EXTENSION_MANIFEST,
    )
    return {
        "kind": "order10_localized_lower_bridge_run_contract",
        "date": "2026-07-16",
        "parameters": {
            "start_t": str(START_T),
            "near_end_t": str(NEAR_END_T),
            "middle_end_t": str(MIDDLE_END_T),
            "end_t": str(END_T),
            "cell_width": str(CELL_WIDTH),
            "quarter_width": str(QUARTER_WIDTH),
            "root_segment_width": str(ROOT_SEGMENT_WIDTH),
            "point_order": POINT_ORDER,
            "remainder_order": REMAINDER_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "precision_bits": PRECISION_BITS,
            "segment_count": len(deterministic_segments()),
        },
        "sources": [
            {"path": relative(path), "sha256": sha256(path)}
            for path in paths
        ],
    }


def ensure_run_contract(path: Path, canonical: dict, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != canonical:
            raise RuntimeError(
                "lower-bridge source/parameter contract changed; use --overwrite"
            )
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(canonical, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_segment_cache(path: Path, tasks: list[tuple]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("lower-bridge segment cache has too many rows")
    for record, task in zip(records, tasks):
        index, regime, left, right = task
        if (
            record.get("kind") != "order10_localized_lower_bridge_segment"
            or record.get("index") != index
            or record.get("regime") != regime
            or record.get("segment_left") != str(left)
            or record.get("segment_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"invalid lower-bridge segment row {index}")
    return records


def build_segment_cache(
    path: Path,
    tasks: list[tuple],
    *,
    workers: int,
    overwrite: bool,
    max_segments: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_segment_cache(path, tasks)
    stop = len(tasks) if max_segments is None else min(len(tasks), max_segments)
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker(
            str(NEAR_H_CACHE),
            str(BROAD_H_CACHE),
            str(POINT_CACHE),
            str(POINT_EXTENSION),
        )
        iterator = map(segment_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(
                str(NEAR_H_CACHE),
                str(BROAD_H_CACHE),
                str(POINT_CACHE),
                str(POINT_EXTENSION),
            ),
        )
        iterator = executor.map(segment_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 5 == 0:
                    handle.flush()
                    print(
                        "order-ten lower-bridge segments: "
                        f"{len(records)}/{stop} ({perf_counter() - started:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_SEGMENT_CACHE)
    parser.add_argument("--run-contract", type=Path, default=DEFAULT_RUN_CONTRACT)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-segments", type=int)
    args = parser.parse_args()
    contract = canonical_run_contract()
    ensure_run_contract(args.run_contract, contract, overwrite=args.overwrite)
    tasks = deterministic_segments()
    records = build_segment_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_segments=args.max_segments,
    )
    print(f"order-ten lower-bridge segment rows: {len(records)}/{len(tasks)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
