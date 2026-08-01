#!/usr/bin/env python3
"""Build order-twelve lower-bridge segments from the inherited sparse H23 source."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from concurrent.futures import ProcessPoolExecutor
import ctypes
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
import mmap
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

import flint  # noqa: E402

import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as compact_h_source  # noqa: E402
import jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache as sparse_builder  # noqa: E402
import jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments as inherited  # noqa: E402
import jensen_window_pf_compound_order12_sparse_h0_h23_propagation_core as propagation  # noqa: E402
import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
from jensen_window_pf_compound_order12_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    DERIVATIVE_MODEL_DEGREES,
    MAXIMUM_H_ORDER as MODEL_MAXIMUM_H_ORDER,
    PRECISION_BITS,
    STABLE_TAYLOR_SURPLUS,
    shifted_taylor_model_curvature_row,
)


BROAD_H_CACHE = inherited.BROAD_H_CACHE
BROAD_H_MANIFEST = inherited.BROAD_H_MANIFEST
SPARSE_CACHE = inherited.SPARSE_CACHE
SPARSE_MANIFEST = inherited.SPARSE_MANIFEST
COMPACT_H_CACHE = compact_h_source.DEFAULT_CACHE
COMPACT_H_MANIFEST = compact_h_source.DEFAULT_MANIFEST
DEFAULT_SEGMENT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.jsonl"
)
DEFAULT_RUN_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_run_contract.json"
)
DEFAULT_PILOT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_pilot.json"
)
PROPAGATION_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order12_sparse_h0_h23_propagation_core.py"
)
PROPAGATION_CHECKER = (
    SCRIPT_DIR / "check_jensen_window_pf_compound_order12_sparse_h0_h23_propagation.py"
)
TAYLOR_MODEL_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order12_shifted_taylor_model_core.py"
)
TAYLOR_MODEL_BASE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_shifted_taylor_model_core.py"
)
TAYLOR_MODEL_CHECKER = (
    SCRIPT_DIR / "check_jensen_window_pf_compound_order12_shifted_taylor_model_algebra.py"
)
PROPAGATION_INTERVAL_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order10_localized_final_gap_interval_core.py"
)
ARB_BOUND_HELPER = (
    SCRIPT_DIR
    / "jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py"
)
SPARSE_GENERATOR = Path(sparse_builder.__file__).resolve()
BROAD_H_GENERATOR = Path(inherited.base.__file__).resolve()
INHERITED_LOWER_DRIVER = Path(inherited.__file__).resolve()
COMPACT_H_GENERATOR = Path(compact_h_source.__file__).resolve()
COMPACT_H_CHECKER = (
    SCRIPT_DIR / "check_jensen_window_pf_compound_order10_compact_h2_h24_unit_cache.py"
)
DRIVER = Path(__file__).resolve()
CHECKER = SCRIPT_DIR / "check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py"

START_T = Fraction(1503)
END_T = Fraction(5700)
CELL_WIDTH = Fraction(1, 2)
QUARTER_WIDTH = Fraction(1, 4)
ROOT_SEGMENT_WIDTH = Fraction(16)
STENCIL_RADIUS = 10
MAXIMUM_PROPAGATION_DISTANCE = Fraction(2)
EXPECTED_SEGMENTS = 263
EXPECTED_HALF_CELLS = 8394
EXPECTED_QUARTER_BLOCKS = 16788
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))

SPARSE_TASKS = sparse_builder.deterministic_tasks(
    sparse_builder.DEFAULT_START_T,
    sparse_builder.DEFAULT_END_T,
    sparse_builder.DEFAULT_STEP_T,
)
SPARSE_ANCHORS = [task[1] for task in SPARSE_TASKS]

_SPARSE_RECORDS: dict[Fraction, dict] = {}
_COMPACT_H_HANDLE = None
_COMPACT_H_VIEW = None
_COMPACT_H_OFFSETS: list[int] = []


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


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def _set_below_normal_priority() -> None:
    if os.name != "nt":
        return
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    kernel32.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
    kernel32.SetPriorityClass.restype = ctypes.c_int
    if not kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), 0x00004000):
        raise ctypes.WinError(ctypes.get_last_error())


def _open_jsonl_mmap(path: Path) -> tuple[object, mmap.mmap, list[int]]:
    handle = path.open("rb")
    view = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
    offsets = []
    position = 0
    while position < len(view):
        offsets.append(position)
        newline = view.find(b"\n", position)
        position = len(view) if newline < 0 else newline + 1
    return handle, view, offsets


def _mmap_record(view: mmap.mmap, offsets: list[int], index: int) -> dict:
    if not 0 <= index < len(offsets):
        raise RuntimeError(f"JSONL row {index} is unavailable")
    end = offsets[index + 1] if index + 1 < len(offsets) else len(view)
    return json.loads(view[offsets[index] : end].strip())


def initialize_worker(broad_h: str, sparse_cache: str, compact_h: str) -> None:
    global _SPARSE_RECORDS
    global _COMPACT_H_HANDLE, _COMPACT_H_VIEW, _COMPACT_H_OFFSETS
    _set_below_normal_priority()
    flint.ctx.prec = PRECISION_BITS

    handle, view, offsets = inherited.base._open_jsonl_mmap(Path(broad_h))
    if len(offsets) != inherited.base.BROAD_H_ROWS:
        raise RuntimeError("broad H2-H24 source row count changed")
    inherited.base._HANDLES["broad"] = handle
    inherited.base._VIEWS["broad"] = view
    inherited.base._OFFSETS["broad"] = offsets
    inherited.base._load_h_row.cache_clear()

    records = sparse_builder.load_cache(Path(sparse_cache), SPARSE_TASKS)
    if len(records) != len(SPARSE_TASKS):
        raise RuntimeError("inherited sparse H0-H23 source is incomplete")
    _SPARSE_RECORDS = {
        Fraction(record["target_t"]): record for record in records
    }
    _load_sparse_anchor.cache_clear()

    (
        _COMPACT_H_HANDLE,
        _COMPACT_H_VIEW,
        _COMPACT_H_OFFSETS,
    ) = _open_jsonl_mmap(Path(compact_h))
    expected = int(
        (compact_h_source.DEFAULT_END_T - compact_h_source.DEFAULT_START_T)
        / compact_h_source.DEFAULT_TILE_WIDTH_T
    )
    if len(_COMPACT_H_OFFSETS) != expected:
        raise RuntimeError("compact H2-H24 source row count changed")
    _load_compact_h_tile.cache_clear()


@lru_cache(maxsize=2400)
def _load_sparse_anchor(anchor: Fraction) -> tuple[list[flint.arb], dict]:
    try:
        record = _SPARSE_RECORDS[anchor]
    except KeyError as exc:
        raise RuntimeError(f"sparse exact H source misses anchor t={anchor}") from exc
    derivatives = record["h_derivatives"]
    series = [
        compact.interval_from_text(derivatives[str(order)]) / math.factorial(order)
        for order in range(propagation.ANCHOR_MAXIMUM_H_ORDER + 1)
    ]
    diagnostics = {
        "target_t": str(anchor),
        "contract_id": record["contract_id"],
        "source_row_sha256": record_sha256(record),
        "mode_bracket": [record["mode_left"], record["mode_right"]],
        "maximum_panel_error_upper": record["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": record["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
    }
    return series, diagnostics


@lru_cache(maxsize=64)
def _load_compact_h_tile(left: Fraction) -> dict:
    index_fraction = (
        (left - compact_h_source.DEFAULT_START_T)
        / compact_h_source.DEFAULT_TILE_WIDTH_T
    )
    if index_fraction.denominator != 1:
        raise RuntimeError(f"unaligned compact H tile {left}")
    index = index_fraction.numerator
    record = _mmap_record(_COMPACT_H_VIEW, _COMPACT_H_OFFSETS, index)
    right = left + compact_h_source.DEFAULT_TILE_WIDTH_T
    derivatives = record.get("h_derivatives", {})
    if (
        record.get("contract_id") != compact_h_source.ROW_CONTRACT
        or record.get("index") != index
        or record.get("target_t_left") != str(left)
        or record.get("target_t_right") != str(right)
        or record.get("passed") is not True
        or set(derivatives) != {str(order) for order in range(2, 25)}
    ):
        raise RuntimeError(f"invalid compact H2-H24 tile at t={left}")
    return {
        "target_t_left": left,
        "target_t_right": right,
        "H": {
            order: compact.interval_from_text(derivatives[str(order)])
            for order in range(2, 25)
        },
        "source_kind": "compact_unit",
    }


def anchor_for_target(target: Fraction) -> Fraction:
    index = bisect_left(SPARSE_ANCHORS, target)
    candidates = SPARSE_ANCHORS[
        max(0, index - 1) : min(len(SPARSE_ANCHORS), index + 2)
    ]
    anchor = min(candidates, key=lambda value: (abs(value - target), value))
    if abs(anchor - target) > MAXIMUM_PROPAGATION_DISTANCE:
        raise RuntimeError(f"sparse H0-H23 anchor is too far from {target}: {anchor}")
    return anchor


def _floor(value: Fraction) -> int:
    return value.numerator // value.denominator


def _ceil(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def h_rows_for_cell(left: Fraction, right: Fraction) -> list[dict]:
    targets = [
        expansion + shift
        for expansion in (left, right)
        for shift in range(-STENCIL_RADIUS, STENCIL_RADIUS + 1)
    ]
    anchors = [anchor_for_target(target) for target in targets]
    support_left = min(left - STENCIL_RADIUS, *anchors)
    support_right = max(right + STENCIL_RADIUS, *anchors)
    if support_left < inherited.base.BROAD_H_START:
        raise RuntimeError(f"lower H collar starts before the broad source: {support_left}")
    if support_right > Fraction(5710):
        raise RuntimeError(f"lower H collar exceeds the certified handoff: {support_right}")

    broad_right = min(support_right, inherited.base.BROAD_H_END)
    first = inherited.base._aligned_index(
        support_left,
        inherited.base.BROAD_H_START,
        inherited.base.BROAD_H_STEP,
    )
    stop = inherited.base._aligned_index(
        broad_right,
        inherited.base.BROAD_H_START,
        inherited.base.BROAD_H_STEP,
    )
    rows = [
        {
            **inherited.base._load_h_row("broad", index),
            "source_kind": "broad_tenth",
        }
        for index in range(first, stop)
    ]
    if support_right > inherited.base.BROAD_H_END:
        compact_first = _floor(
            max(inherited.base.BROAD_H_END, support_left)
        )
        compact_stop = _ceil(support_right)
        rows.extend(
            _load_compact_h_tile(Fraction(tile_left))
            for tile_left in range(compact_first, compact_stop)
        )
    rows.sort(key=lambda row: row["target_t_left"])
    if (
        not rows
        or rows[0]["target_t_left"] > support_left
        or rows[-1]["target_t_right"] < support_right
        or any(
            row["target_t_right"] < next_row["target_t_left"]
            for row, next_row in zip(rows, rows[1:])
        )
    ):
        raise RuntimeError(f"mixed H2-H24 collar is incomplete on {support_left}..{support_right}")
    return rows


def point_source_for(
    expansion_anchor: Fraction,
    h_rows: list[dict],
) -> tuple[dict, dict]:
    targets = [
        expansion_anchor + shift
        for shift in range(-STENCIL_RADIUS, STENCIL_RADIUS + 1)
    ]
    expected_anchors = [anchor_for_target(target) for target in targets]
    anchor_source = {}
    anchor_diagnostics = {}
    for anchor in sorted(set(expected_anchors)):
        series, diagnostics = _load_sparse_anchor(anchor)
        anchor_source[anchor] = series
        anchor_diagnostics[anchor] = diagnostics
    propagated = propagation.propagated_point_source(
        targets,
        anchor_source,
        h_rows,
        anchor_diagnostic_source=anchor_diagnostics,
    )
    actual_anchors = [
        Fraction(propagated[target][1]["sparse_h23_propagation"]["anchor_t"])
        for target in targets
    ]
    if actual_anchors != expected_anchors:
        raise RuntimeError(f"propagated anchor selection changed at {expansion_anchor}")
    metadata = {
        "point_source_targets": [str(target) for target in targets],
        "point_source_anchors": [str(anchor) for anchor in actual_anchors],
        "point_source_contract_ids": [
            _SPARSE_RECORDS[anchor]["contract_id"] for anchor in actual_anchors
        ],
        "maximum_actual_propagation_distance": str(
            max(abs(target - anchor) for target, anchor in zip(targets, actual_anchors))
        ),
        "h_source_kinds": sorted({row["source_kind"] for row in h_rows}),
    }
    return propagated, metadata


def deterministic_segments() -> list[tuple[int, Fraction, Fraction]]:
    tasks = []
    left = START_T
    while left < END_T:
        right = min(left + ROOT_SEGMENT_WIDTH, END_T)
        tasks.append((len(tasks), left, right))
        left = right
    if (
        len(tasks) != EXPECTED_SEGMENTS
        or tasks[-1][2] != END_T
        or sum((right - left for _, left, right in tasks), Fraction())
        != END_T - START_T
    ):
        raise RuntimeError("order-twelve lower segmentation changed")
    return tasks


def segment_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, segment_left, segment_right = task
    flint.ctx.prec = PRECISION_BITS
    blocks = []
    cell_left = segment_left
    while cell_left < segment_right:
        cell_right = cell_left + CELL_WIDTH
        if cell_right > segment_right:
            raise RuntimeError(f"lower segment {index} cuts a half-cell")
        h_rows = h_rows_for_cell(cell_left, cell_right)
        for expansion, left, right in (
            (cell_left, cell_left, cell_left + QUARTER_WIDTH),
            (cell_right, cell_left + QUARTER_WIDTH, cell_right),
        ):
            propagated, metadata = point_source_for(expansion, h_rows)
            block = shifted_taylor_model_curvature_row(
                expansion,
                left,
                right,
                h_rows,
                point_h_source=propagated,
            )
            blocks.append({**block, **metadata, "regime": "lower"})
        cell_left = cell_right
    expected_blocks = int((segment_right - segment_left) / QUARTER_WIDTH)
    if len(blocks) != expected_blocks:
        raise RuntimeError(f"lower segment {index} block count changed")
    if any(block.get("passed") is not True for block in blocks):
        raise RuntimeError(f"lower segment {index} contains a failing block")
    maximum = max(blocks, key=lambda block: float(block["scaled_curvature_upper"]))
    minimum = min(blocks, key=lambda block: float(block["curvature_margin_lower"]))
    return {
        "kind": "order12_sparse_h23_lower_bridge_segment",
        "index": index,
        "segment_left": str(segment_left),
        "segment_right": str(segment_right),
        "blocks": blocks,
        "block_count": len(blocks),
        "point_source_contract_ids": sorted(
            {
                contract
                for block in blocks
                for contract in block["point_source_contract_ids"]
            }
        ),
        "h_source_kinds": sorted(
            {kind for block in blocks for kind in block["h_source_kinds"]}
        ),
        "largest_scaled_curvature_upper": maximum["scaled_curvature_upper"],
        "largest_scaled_curvature_anchor": maximum["anchor"],
        "smallest_margin_lower": minimum["curvature_margin_lower"],
        "smallest_margin_anchor": minimum["anchor"],
        "passed": True,
    }


def validate_broad_h_source() -> None:
    inherited.base.validate_source_manifest(
        BROAD_H_MANIFEST,
        BROAD_H_CACHE,
        start=str(inherited.base.BROAD_H_START),
        end=str(inherited.base.BROAD_H_END),
        step=str(inherited.base.BROAD_H_STEP),
        rows=inherited.base.BROAD_H_ROWS,
        max_moment=24,
    )


def validate_compact_h_source() -> dict:
    import check_jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as checker

    manifest = json.loads(COMPACT_H_MANIFEST.read_text(encoding="utf-8"))
    cache = manifest.get("cache", {})
    expected_rows = int(
        (compact_h_source.DEFAULT_END_T - compact_h_source.DEFAULT_START_T)
        / compact_h_source.DEFAULT_TILE_WIDTH_T
    )
    if (
        manifest.get("parameters") != checker.expected_parameters()
        or cache.get("path") != relative(COMPACT_H_CACHE)
        or cache.get("row_count") != expected_rows
        or cache.get("all_rows_passed") is not True
        or cache.get("h_derivative_orders") != [2, 24]
        or cache.get("sha256") != sha256(COMPACT_H_CACHE)
    ):
        raise RuntimeError("compact H2-H24 source manifest changed")
    return manifest


def canonical_run_contract() -> dict:
    validate_broad_h_source()
    sparse_manifest = inherited.validate_sparse_manifest()
    compact_manifest = validate_compact_h_source()
    paths = (
        BROAD_H_CACHE,
        BROAD_H_MANIFEST,
        SPARSE_CACHE,
        SPARSE_MANIFEST,
        COMPACT_H_CACHE,
        COMPACT_H_MANIFEST,
        BROAD_H_GENERATOR,
        INHERITED_LOWER_DRIVER,
        SPARSE_GENERATOR,
        COMPACT_H_GENERATOR,
        COMPACT_H_CHECKER,
        PROPAGATION_CORE,
        PROPAGATION_CHECKER,
        TAYLOR_MODEL_CORE,
        TAYLOR_MODEL_BASE,
        TAYLOR_MODEL_CHECKER,
        PROPAGATION_INTERVAL_CORE,
        ARB_BOUND_HELPER,
        DRIVER,
        CHECKER,
    )
    return {
        "kind": "order12_sparse_h23_lower_bridge_run_contract",
        "date": "2026-07-22",
        "status": "hash-bound mixed-source order-twelve lower-bridge contract",
        "proof_boundary": (
            "This contract binds the physical lower first-summand computation only. "
            "It does not certify a segment, the compact range, PF-infinity, or RH."
        ),
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "cell_width": str(CELL_WIDTH),
            "quarter_width": str(QUARTER_WIDTH),
            "root_segment_width": str(ROOT_SEGMENT_WIDTH),
            "segment_count": EXPECTED_SEGMENTS,
            "half_cell_count": EXPECTED_HALF_CELLS,
            "quarter_block_count": EXPECTED_QUARTER_BLOCKS,
            "stencil_radius": STENCIL_RADIUS,
            "maximum_propagation_distance": str(MAXIMUM_PROPAGATION_DISTANCE),
            "broad_to_compact_handoff_t": str(inherited.base.BROAD_H_END),
            "derivative_model_degrees": list(DERIVATIVE_MODEL_DEGREES),
            "stable_taylor_surplus": STABLE_TAYLOR_SURPLUS,
            "model_maximum_h_order": MODEL_MAXIMUM_H_ORDER,
            "anchor_maximum_h_order": propagation.ANCHOR_MAXIMUM_H_ORDER,
            "propagated_maximum_h_order": propagation.OUTPUT_MAXIMUM_H_ORDER,
            "remainder_h_order": propagation.REMAINDER_H_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "precision_bits": PRECISION_BITS,
        },
        "source_completion": {
            "sparse_rows": sparse_manifest["cache"]["row_count"],
            "sparse_cache_sha256": sparse_manifest["cache"]["sha256"],
            "broad_h_rows": inherited.base.BROAD_H_ROWS,
            "compact_h_rows": compact_manifest["cache"]["row_count"],
        },
        "sources": [
            {"path": relative(path), "sha256": sha256(path)} for path in paths
        ],
    }


def ensure_run_contract(path: Path, canonical: dict, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != canonical:
            raise RuntimeError("order-twelve lower contract changed; use --overwrite")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(canonical, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_segment_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if count >= len(tasks):
                raise RuntimeError("order-twelve lower cache has too many rows")
            record = json.loads(line)
            index, left, right = tasks[count]
            if (
                record.get("kind") != "order12_sparse_h23_lower_bridge_segment"
                or record.get("index") != index
                or record.get("segment_left") != str(left)
                or record.get("segment_right") != str(right)
                or record.get("passed") is not True
            ):
                raise RuntimeError(
                    f"invalid order-twelve lower segment {index} at line {line_number}"
                )
            count += 1
    return count


def _initializer_args() -> tuple[str, str, str]:
    return str(BROAD_H_CACHE), str(SPARSE_CACHE), str(COMPACT_H_CACHE)


def build_segment_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    overwrite: bool,
    max_segments: int | None,
    runtime_seconds: float | None,
) -> int:
    if overwrite and path.exists():
        path.unlink()
    count = load_segment_cache(path, tasks)
    stop = len(tasks) if max_segments is None else min(len(tasks), max_segments)
    remaining = tasks[count:stop]
    if not remaining:
        return count
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    executor = None
    try:
        if workers > 1:
            executor = ProcessPoolExecutor(
                max_workers=workers,
                initializer=initialize_worker,
                initargs=_initializer_args(),
            )
        else:
            initialize_worker(*_initializer_args())
        with path.open("a", encoding="utf-8") as handle:
            cursor = 0
            while cursor < len(remaining):
                batch = remaining[cursor : cursor + workers]
                if executor is None:
                    completed = [segment_task(task) for task in batch]
                else:
                    completed = list(executor.map(segment_task, batch, chunksize=1))
                for record in completed:
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    count += 1
                handle.flush()
                cursor += len(batch)
                elapsed = perf_counter() - started
                print(
                    f"order-twelve sparse lower segments: {count}/{stop} ({elapsed:.1f}s)",
                    flush=True,
                )
                if runtime_seconds is not None and elapsed >= runtime_seconds:
                    print(
                        "order-twelve lower runtime checkpoint: parking after the completed batch",
                        flush=True,
                    )
                    break
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return count


def build_pilot(indices: list[int], *, workers: int) -> dict:
    canonical = canonical_run_contract()
    tasks = deterministic_segments()
    unique = sorted(set(indices))
    if not unique or any(not 0 <= index < len(tasks) for index in unique):
        raise ValueError("pilot segment index leaves the lower bridge")
    selected = [tasks[index] for index in unique]
    if workers > 1 and len(selected) > 1:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(selected)),
            initializer=initialize_worker,
            initargs=_initializer_args(),
        ) as executor:
            records = list(executor.map(segment_task, selected, chunksize=1))
    else:
        initialize_worker(*_initializer_args())
        records = [segment_task(task) for task in selected]
    maximum = max(
        (block for record in records for block in record["blocks"]),
        key=lambda block: float(block["scaled_curvature_upper"]),
    )
    minimum = min(
        (block for record in records for block in record["blocks"]),
        key=lambda block: float(block["curvature_margin_lower"]),
    )
    return {
        "kind": "order12_sparse_h23_lower_bridge_pilot",
        "date": "2026-07-22",
        "status": "rigorous named lower-segment pilots over hash-bound mixed sources",
        "proof_boundary": (
            "This certifies only the named lower segments. It is not contiguous "
            "lower coverage, an order-twelve theorem, or a proof of RH."
        ),
        "run_contract": canonical,
        "run_contract_sha256": record_sha256(canonical),
        "segment_indices": unique,
        "segments": records,
        "summary": {
            "segments": len(records),
            "quarter_blocks": sum(record["block_count"] for record in records),
            "largest_scaled_curvature_upper": maximum["scaled_curvature_upper"],
            "smallest_margin_lower": minimum["curvature_margin_lower"],
            "failed_segments": 0,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_SEGMENT_CACHE)
    parser.add_argument("--run-contract", type=Path, default=DEFAULT_RUN_CONTRACT)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-segments", type=int)
    parser.add_argument("--runtime-seconds", type=float)
    parser.add_argument("--pilot-segment", action="append", type=int, default=[])
    parser.add_argument("--pilot-output", type=Path, default=DEFAULT_PILOT)
    args = parser.parse_args()
    if args.runtime_seconds is not None and args.runtime_seconds <= 0:
        parser.error("--runtime-seconds must be positive")
    if args.max_segments is not None and args.max_segments < 0:
        parser.error("--max-segments cannot be negative")
    workers = max(1, min(4, args.workers))
    if args.pilot_segment:
        artifact = build_pilot(args.pilot_segment, workers=workers)
        args.pilot_output.parent.mkdir(parents=True, exist_ok=True)
        args.pilot_output.write_text(
            json.dumps(artifact, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(
            "wrote order-twelve sparse lower pilots: "
            f"{artifact['summary']['segments']} segments, "
            f"{artifact['summary']['quarter_blocks']} quarter blocks"
        )
        return 0
    canonical = canonical_run_contract()
    ensure_run_contract(args.run_contract, canonical, overwrite=args.overwrite)
    tasks = deterministic_segments()
    count = build_segment_cache(
        args.cache,
        tasks,
        workers=workers,
        overwrite=args.overwrite,
        max_segments=args.max_segments,
        runtime_seconds=args.runtime_seconds,
    )
    print(f"order-twelve sparse lower segment rows: {count}/{len(tasks)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
