#!/usr/bin/env python3
"""Retained finite-H8-stencil experiment for the order-eleven lower bridge.

This route closes low blocks but loses cancellation at larger t.  It is not the
canonical order-eleven driver; use the sparse-H23 driver for proof work.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
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

import jensen_window_pf_compound_order10_localized_lower_bridge_segments as base  # noqa: E402
from jensen_window_pf_compound_order11_point_h9_h16_stencil_core import (  # noqa: E402
    BASE_H_ORDER,
    MAXIMUM_H_ORDER as STENCIL_MAXIMUM_H_ORDER,
    NODE_COUNT as STENCIL_NODE_COUNT,
    REMAINDER_H_ORDER as STENCIL_REMAINDER_H_ORDER,
    higher_derivative_stencil_jet,
)
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    DERIVATIVE_MODEL_DEGREES,
    MAXIMUM_H_ORDER as MODEL_MAXIMUM_H_ORDER,
    PRECISION_BITS,
    STABLE_TAYLOR_SURPLUS,
    shifted_taylor_model_curvature_row,
)


NEAR_H_CACHE = base.NEAR_H_CACHE
NEAR_H_MANIFEST = base.NEAR_H_MANIFEST
BROAD_H_CACHE = base.BROAD_H_CACHE
BROAD_H_MANIFEST = base.BROAD_H_MANIFEST
POINT_CACHE = base.POINT_CACHE
POINT_MANIFEST = base.POINT_MANIFEST
POINT_EXTENSION = base.POINT_EXTENSION
POINT_EXTENSION_MANIFEST = base.POINT_EXTENSION_MANIFEST
POINT_FINAL_EXTENSION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_shifted_point_h0_h8_half_grid_extension_5708p5.jsonl"
)
POINT_FINAL_EXTENSION_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_shifted_point_h0_h8_half_grid_extension_5708p5_cache.json"
)
POINT_FURTHER_EXTENSION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_shifted_point_h0_h8_half_grid_extension_5709_5710.jsonl"
)
POINT_FURTHER_EXTENSION_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_shifted_point_h0_h8_half_grid_extension_5709_5710_cache.json"
)
DEFAULT_SEGMENT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_localized_lower_bridge_segments.jsonl"
)
DEFAULT_RUN_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_localized_lower_bridge_run_contract.json"
)
TAYLOR_MODEL_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_shifted_taylor_model_core.py"
)
STENCIL_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_point_h9_h16_stencil_core.py"
)

START_T = Fraction(1252)
NEAR_END_T = Fraction(1500)
MIDDLE_END_T = Fraction(1800)
END_T = Fraction(5700)
CELL_WIDTH = Fraction(1, 2)
QUARTER_WIDTH = Fraction(1, 4)
ROOT_SEGMENT_WIDTH = Fraction(16)
FINAL_POINT_T = Fraction(11417, 2)
FINAL_POINT_ROWS = 1
FURTHER_POINT_START = Fraction(5709)
FURTHER_POINT_END = Fraction(5710)
FURTHER_POINT_ROWS = 3
POINT_ORDERS = tuple(range(9))
STENCIL_POINT_START = base.POINT_START
STENCIL_POINT_END = FURTHER_POINT_START
STENCIL_POINT_STEP = Fraction(1, 2)
STENCIL_POINT_ROWS = int(
    (STENCIL_POINT_END - STENCIL_POINT_START) / STENCIL_POINT_STEP
) + 1
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))


_ADDITIONAL_POINTS: dict[Fraction, dict] = {}


def sha256(path: Path) -> str:
    return base.sha256(path)


def relative(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def initialize_worker(
    near_h: str,
    broad_h: str,
    point: str,
    point_extension: str,
    point_final_extension: str,
    point_further_extension: str,
) -> None:
    global _ADDITIONAL_POINTS
    flint.ctx.prec = PRECISION_BITS
    base.initialize_worker(near_h, broad_h, point, point_extension)
    expected = [FINAL_POINT_T] + [
        FURTHER_POINT_START + index * STENCIL_POINT_STEP
        for index in range(FURTHER_POINT_ROWS)
    ]
    paths = [Path(point_final_extension), Path(point_further_extension)]
    records = []
    for path in paths:
        records.extend(
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
    if len(records) != len(expected):
        raise RuntimeError("order-eleven terminal point sources have the wrong size")
    additional = {}
    for target, record in zip(expected, records):
        derivatives = record.get("h_derivatives", {})
        if (
            record.get("target_t") != str(target)
            or record.get("passed") is not True
            or set(derivatives) != {str(order) for order in POINT_ORDERS}
        ):
            raise RuntimeError(f"invalid order-eleven terminal point source at {target}")
        additional[target] = record
    _ADDITIONAL_POINTS = additional
    _load_point.cache_clear()
    _augmented_point.cache_clear()


@lru_cache(maxsize=6000)
def _load_point(target: Fraction) -> tuple[list, dict]:
    if target <= base.POINT_EXTENSION_END:
        return base._load_point(target)
    try:
        record = _ADDITIONAL_POINTS[target]
    except KeyError as exc:
        raise RuntimeError(f"order-eleven exact point source misses t={target}")
    derivatives = record["h_derivatives"]
    return (
        [
            base.compact.interval_from_text(derivatives[str(order)])
            / math.factorial(order)
            for order in POINT_ORDERS
        ],
        {
            "target_t": str(target),
            "mode_bracket": [
                record["mode_left"],
                record["mode_right"],
            ],
            "maximum_panel_error_upper": record[
                "maximum_panel_error_upper"
            ],
            "maximum_tail_moment_upper": record[
                "maximum_tail_moment_upper"
            ],
            "minimum_tail_slope_lower": record[
                "minimum_tail_slope_lower"
            ],
        },
    )


def _stencil_nodes(target: Fraction) -> list[Fraction]:
    index = base._aligned_index(target, STENCIL_POINT_START, STENCIL_POINT_STEP)
    if not 0 <= index < STENCIL_POINT_ROWS:
        raise RuntimeError(f"higher-derivative stencil target leaves grid: {target}")
    start = min(max(index - 7, 0), STENCIL_POINT_ROWS - STENCIL_NODE_COUNT)
    return [
        STENCIL_POINT_START + (start + offset) * STENCIL_POINT_STEP
        for offset in range(STENCIL_NODE_COUNT)
    ]


def _h24_rows_for_stencil(nodes: list[Fraction]) -> list[dict]:
    left, right = nodes[0], nodes[-1]
    if left < base.BROAD_H_START:
        kind, start, step, row_count = (
            "near",
            base.NEAR_H_START,
            base.NEAR_H_STEP,
            base.NEAR_H_ROWS,
        )
    else:
        kind, start, step, row_count = (
            "broad",
            base.BROAD_H_START,
            base.BROAD_H_STEP,
            base.BROAD_H_ROWS,
        )
    first = base._aligned_index(left, start, step)
    last = base._aligned_index(right, start, step)
    if not 0 <= first < last <= row_count:
        raise RuntimeError(f"H24 stencil support leaves {kind} cache: {left}..{right}")
    return [base._load_h_row(kind, index) for index in range(first, last)]


@lru_cache(maxsize=12000)
def _augmented_point(target: Fraction) -> tuple[list, dict]:
    nodes = _stencil_nodes(target)
    raw_source = {node: _load_point(node) for node in nodes}
    return higher_derivative_stencil_jet(
        target,
        nodes,
        _h24_rows_for_stencil(nodes),
        point_h_source=raw_source,
    )


def point_source_for(anchor: Fraction) -> dict:
    return {
        anchor + shift: _augmented_point(anchor + shift)
        for shift in range(-9, 10)
    }


def h_rows_for_cell(regime: str, left: Fraction, right: Fraction) -> list[dict]:
    if regime == "near":
        kind, start, step, row_count = (
            "near",
            base.NEAR_H_START,
            base.NEAR_H_STEP,
            base.NEAR_H_ROWS,
        )
    else:
        kind, start, step, row_count = (
            "broad",
            base.BROAD_H_START,
            base.BROAD_H_STEP,
            base.BROAD_H_ROWS,
        )
    support_left = left - 9
    support_right = right + 9
    first = base._aligned_index(support_left, start, step)
    last = base._aligned_index(support_right, start, step)
    if not 0 <= first < last <= row_count:
        raise RuntimeError(
            f"{regime} H collar {support_left}..{support_right} leaves its cache"
        )
    return [base._load_h_row(kind, index) for index in range(first, last)]


def deterministic_segments() -> list[tuple[int, str, Fraction, Fraction]]:
    result = []
    index = 0
    for regime, start, end in (
        ("near", START_T, NEAR_END_T),
        ("middle", NEAR_END_T, MIDDLE_END_T),
        ("far", MIDDLE_END_T, END_T),
    ):
        left = start
        while left < end:
            right = min(left + ROOT_SEGMENT_WIDTH, end)
            result.append((index, regime, left, right))
            index += 1
            left = right
    return result


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
        specs = (
            (cell_left, cell_left, cell_left + QUARTER_WIDTH),
            (cell_right, cell_left + QUARTER_WIDTH, cell_right),
        )
        for expansion, left, right in specs:
            block = shifted_taylor_model_curvature_row(
                expansion,
                left,
                right,
                h_rows,
                point_h_source=point_source_for(expansion),
            )
            blocks.append({**block, "regime": regime})
        cell_left = cell_right
    if not blocks or any(row.get("passed") is not True for row in blocks):
        raise RuntimeError(f"segment {index} did not produce only passing blocks")
    maximum = max(blocks, key=lambda row: float(row["scaled_curvature_upper"]))
    minimum = min(blocks, key=lambda row: float(row["curvature_margin_lower"]))
    return {
        "kind": "order11_localized_lower_bridge_segment",
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


def _validate_point_manifest(
    manifest_path: Path,
    cache_path: Path,
    *,
    rows: int,
    start: str,
    end: str,
) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cache = manifest.get("cache", {})
    parameters = manifest.get("parameters", {})
    if (
        cache.get("row_count") != rows
        or cache.get("all_rows_passed") is not True
        or cache.get("sha256") != sha256(cache_path)
        or parameters.get("start_t") != start
        or parameters.get("end_t") != end
        or parameters.get("step_t") != "1/2"
        or parameters.get("max_moment") != 8
        or parameters.get("precision_bits") != PRECISION_BITS
    ):
        raise RuntimeError(f"point source manifest changed: {manifest_path}")


def canonical_run_contract() -> dict:
    base.validate_source_manifest(
        NEAR_H_MANIFEST,
        NEAR_H_CACHE,
        start=str(base.NEAR_H_START),
        end=str(base.NEAR_H_END),
        step=str(base.NEAR_H_STEP),
        rows=base.NEAR_H_ROWS,
        max_moment=24,
    )
    base.validate_source_manifest(
        BROAD_H_MANIFEST,
        BROAD_H_CACHE,
        start=str(base.BROAD_H_START),
        end=str(base.BROAD_H_END),
        step=str(base.BROAD_H_STEP),
        rows=base.BROAD_H_ROWS,
        max_moment=24,
    )
    _validate_point_manifest(
        POINT_MANIFEST,
        POINT_CACHE,
        rows=base.POINT_ROWS,
        start=str(base.POINT_START),
        end=str(base.POINT_END),
    )
    _validate_point_manifest(
        POINT_EXTENSION_MANIFEST,
        POINT_EXTENSION,
        rows=base.POINT_EXTENSION_ROWS,
        start=str(base.POINT_EXTENSION_START),
        end=str(base.POINT_EXTENSION_END),
    )
    _validate_point_manifest(
        POINT_FINAL_EXTENSION_MANIFEST,
        POINT_FINAL_EXTENSION,
        rows=FINAL_POINT_ROWS,
        start=str(FINAL_POINT_T),
        end=str(FINAL_POINT_T),
    )
    _validate_point_manifest(
        POINT_FURTHER_EXTENSION_MANIFEST,
        POINT_FURTHER_EXTENSION,
        rows=FURTHER_POINT_ROWS,
        start=str(FURTHER_POINT_START),
        end=str(FURTHER_POINT_END),
    )
    paths = (
        NEAR_H_CACHE,
        NEAR_H_MANIFEST,
        BROAD_H_CACHE,
        BROAD_H_MANIFEST,
        POINT_CACHE,
        POINT_MANIFEST,
        POINT_EXTENSION,
        POINT_EXTENSION_MANIFEST,
        POINT_FINAL_EXTENSION,
        POINT_FINAL_EXTENSION_MANIFEST,
        POINT_FURTHER_EXTENSION,
        POINT_FURTHER_EXTENSION_MANIFEST,
        TAYLOR_MODEL_CORE,
        STENCIL_CORE,
    )
    return {
        "kind": "order11_localized_lower_bridge_run_contract",
        "date": "2026-07-17",
        "parameters": {
            "start_t": str(START_T),
            "near_end_t": str(NEAR_END_T),
            "middle_end_t": str(MIDDLE_END_T),
            "end_t": str(END_T),
            "cell_width": str(CELL_WIDTH),
            "quarter_width": str(QUARTER_WIDTH),
            "root_segment_width": str(ROOT_SEGMENT_WIDTH),
            "derivative_model_degrees": list(DERIVATIVE_MODEL_DEGREES),
            "stable_taylor_surplus": STABLE_TAYLOR_SURPLUS,
            "model_maximum_h_order": MODEL_MAXIMUM_H_ORDER,
            "stencil_base_h_order": BASE_H_ORDER,
            "stencil_maximum_h_order": STENCIL_MAXIMUM_H_ORDER,
            "stencil_node_count": STENCIL_NODE_COUNT,
            "stencil_remainder_h_order": STENCIL_REMAINDER_H_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "precision_bits": PRECISION_BITS,
            "segment_count": len(deterministic_segments()),
        },
        "sources": [
            {"path": relative(path), "sha256": sha256(path)} for path in paths
        ],
    }


def ensure_run_contract(path: Path, canonical: dict, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != canonical:
            raise RuntimeError(
                "order-eleven lower-bridge contract changed; use --overwrite"
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
        raise RuntimeError("order-eleven lower-bridge cache has too many rows")
    for record, task in zip(records, tasks):
        index, regime, left, right = task
        if (
            record.get("kind") != "order11_localized_lower_bridge_segment"
            or record.get("index") != index
            or record.get("regime") != regime
            or record.get("segment_left") != str(left)
            or record.get("segment_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"invalid order-eleven lower-bridge row {index}")
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
    initargs = (
        str(NEAR_H_CACHE),
        str(BROAD_H_CACHE),
        str(POINT_CACHE),
        str(POINT_EXTENSION),
        str(POINT_FINAL_EXTENSION),
        str(POINT_FURTHER_EXTENSION),
    )
    if workers == 1:
        initialize_worker(*initargs)
        iterator = map(segment_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=initargs,
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
                        "order-eleven lower-bridge segments: "
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
    parser.add_argument("--allow-experimental-stencil", action="store_true")
    args = parser.parse_args()
    if not args.allow_experimental_stencil:
        parser.error(
            "finite-H8 stencil mode is experimental and fails at high t; "
            "use jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments.py"
        )
    contract = canonical_run_contract()
    ensure_run_contract(args.run_contract, contract, overwrite=args.overwrite)
    tasks = deterministic_segments()
    records = build_segment_cache(
        args.cache,
        tasks,
        workers=max(1, min(4, args.workers)),
        overwrite=args.overwrite,
        max_segments=args.max_segments,
    )
    print(f"order-eleven lower-bridge segment rows: {len(records)}/{len(tasks)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
