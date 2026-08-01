#!/usr/bin/env python3
"""Build canonical compact order-eleven segments from the adaptive H0-H23 source."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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

import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache as point_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension as boundary_source  # noqa: E402
import jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core as propagation  # noqa: E402
import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    DERIVATIVE_MODEL_DEGREES,
    MAXIMUM_H_ORDER as MODEL_MAXIMUM_H_ORDER,
    PRECISION_BITS,
    STABLE_TAYLOR_SURPLUS,
    shifted_taylor_model_curvature_row,
)


POINT_CACHE = point_source.DEFAULT_CACHE
POINT_MANIFEST = point_source.DEFAULT_MANIFEST
H_CACHE = h_source.DEFAULT_CACHE
H_MANIFEST = h_source.DEFAULT_MANIFEST
H_BOUNDARY = boundary_source.DEFAULT_OUT
DEFAULT_SEGMENT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_segments.jsonl"
)
DEFAULT_RUN_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_run_contract.json"
)
DEFAULT_PILOT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_adaptive_h23_segment_pilot.json"
)
PROPAGATION_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core.py"
)
TAYLOR_MODEL_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_shifted_taylor_model_core.py"
)
POINT_CHECKER = (
    SCRIPT_DIR
    / "check_jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache.py"
)
H_CHECKER = (
    SCRIPT_DIR / "check_jensen_window_pf_compound_order10_compact_h2_h24_unit_cache.py"
)
BOUNDARY_CHECKER = (
    SCRIPT_DIR / "check_jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension.py"
)
DRIVER = Path(__file__).resolve()
CHECKER = SCRIPT_DIR / "check_jensen_window_pf_compound_order11_compact_adaptive_h23_segments.py"

START_T = Fraction(5700)
END_T = Fraction(38020)
CELL_WIDTH = Fraction(1, 2)
QUARTER_WIDTH = Fraction(1, 4)
ROOT_SEGMENT_WIDTH = Fraction(16)
STENCIL_RADIUS = 9
MAXIMUM_PROPAGATION_DISTANCE = Fraction(2)
EXPECTED_SEGMENTS = 2020
EXPECTED_HALF_CELLS = 64640
EXPECTED_QUARTER_BLOCKS = 129280
POINT_TASKS = point_source.deterministic_tasks()
POINT_ROWS = len(POINT_TASKS)
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))


_POINT_HANDLE = None
_POINT_VIEW = None
_POINT_OFFSETS: list[int] = []
_H_HANDLE = None
_H_VIEW = None
_H_OFFSETS: list[int] = []
_BOUNDARY_H: dict[Fraction, dict] = {}


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


def _parse_h_record(record: dict, *, expected_left: Fraction, index: int | None) -> dict:
    right = expected_left + 1
    derivatives = record.get("h_derivatives", {})
    if (
        record.get("contract_id") != h_source.ROW_CONTRACT
        or record.get("target_t_left") != str(expected_left)
        or record.get("target_t_right") != str(right)
        or record.get("passed") is not True
        or set(derivatives) != {str(order) for order in range(2, 25)}
        or (index is not None and record.get("index") != index)
    ):
        raise RuntimeError(f"invalid compact H2-H24 tile at t={expected_left}")
    return {
        "target_t_left": expected_left,
        "target_t_right": right,
        "H": {
            order: compact.interval_from_text(derivatives[str(order)])
            for order in range(2, 25)
        },
    }


def _load_boundary_rows(path: Path) -> dict[Fraction, dict]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    rows = artifact.get("rows", [])
    expected = (Fraction(5691), Fraction(38028))
    if len(rows) != 2:
        raise RuntimeError("compact H boundary extension does not have two rows")
    result = {}
    for record, left in zip(rows, expected):
        result[left] = _parse_h_record(record, expected_left=left, index=None)
    if set(result) != set(expected):
        raise RuntimeError("compact H boundary extension coverage changed")
    return result


def initialize_worker(point_cache: str, h_cache: str, boundary_path: str) -> None:
    global _POINT_HANDLE, _POINT_VIEW, _POINT_OFFSETS
    global _H_HANDLE, _H_VIEW, _H_OFFSETS, _BOUNDARY_H
    flint.ctx.prec = PRECISION_BITS
    propagation.MAXIMUM_PROPAGATION_DISTANCE = MAXIMUM_PROPAGATION_DISTANCE
    _POINT_HANDLE, _POINT_VIEW, _POINT_OFFSETS = _open_jsonl_mmap(Path(point_cache))
    _H_HANDLE, _H_VIEW, _H_OFFSETS = _open_jsonl_mmap(Path(h_cache))
    if len(_POINT_OFFSETS) != POINT_ROWS:
        raise RuntimeError("adaptive H0-H23 source row count changed")
    if len(_H_OFFSETS) != int(h_source.DEFAULT_END_T - h_source.DEFAULT_START_T):
        raise RuntimeError("compact H2-H24 source row count changed")
    _BOUNDARY_H = _load_boundary_rows(Path(boundary_path))
    _load_point_anchor.cache_clear()
    _load_h_tile.cache_clear()


def anchor_for_target(target: Fraction) -> Fraction:
    quotient = (target - point_source.START_T) / point_source.STEP_T
    floor_index = quotient.numerator // quotient.denominator
    candidates = {
        max(0, min(POINT_ROWS - 1, floor_index)),
        max(0, min(POINT_ROWS - 1, floor_index + 1)),
    }
    anchors = [
        point_source.START_T + index * point_source.STEP_T
        for index in candidates
    ]
    anchor = min(anchors, key=lambda value: (abs(value - target), value))
    if abs(anchor - target) > MAXIMUM_PROPAGATION_DISTANCE:
        raise RuntimeError(f"adaptive H0-H23 anchor is too far from {target}: {anchor}")
    return anchor


def profile_for_anchor(anchor: Fraction) -> str:
    return point_source.profile_for_target(anchor)


@lru_cache(maxsize=1024)
def _load_point_anchor(anchor: Fraction) -> tuple[list[flint.arb], dict]:
    index_fraction = (anchor - point_source.START_T) / point_source.STEP_T
    if index_fraction.denominator != 1:
        raise RuntimeError(f"unaligned adaptive point anchor {anchor}")
    index = index_fraction.numerator
    record = _mmap_record(_POINT_VIEW, _POINT_OFFSETS, index)
    expected = POINT_TASKS[index]
    point_source.validate_record(record, expected)
    derivatives = record["h_derivatives"]
    series = [
        compact.interval_from_text(derivatives[str(order)]) / math.factorial(order)
        for order in range(propagation.ANCHOR_MAXIMUM_H_ORDER + 1)
    ]
    diagnostics = {
        "target_t": str(anchor),
        "profile": record["profile"],
        "source_row_sha256": record_sha256(record),
        "mode_bracket": [record["mode_left"], record["mode_right"]],
        "maximum_panel_error_upper": record["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": record["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
    }
    return series, diagnostics


@lru_cache(maxsize=512)
def _load_h_tile(left: Fraction) -> dict:
    if left in _BOUNDARY_H:
        return _BOUNDARY_H[left]
    index_fraction = left - h_source.DEFAULT_START_T
    if index_fraction.denominator != 1:
        raise RuntimeError(f"unaligned compact H tile {left}")
    index = index_fraction.numerator
    record = _mmap_record(_H_VIEW, _H_OFFSETS, index)
    return _parse_h_record(record, expected_left=left, index=index)


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
    first = _floor(support_left)
    stop = _ceil(support_right)
    if first < 5691 or stop > 38029:
        raise RuntimeError(f"compact H collar {support_left}..{support_right} leaves source")
    return [_load_h_tile(Fraction(tile_left)) for tile_left in range(first, stop)]


def point_source_for(
    expansion_anchor: Fraction,
    h_rows: list[dict],
) -> tuple[dict, dict]:
    targets = [
        expansion_anchor + shift
        for shift in range(-STENCIL_RADIUS, STENCIL_RADIUS + 1)
    ]
    expected_anchors = [anchor_for_target(target) for target in targets]
    chosen = sorted(set(expected_anchors))
    anchor_source = {}
    anchor_diagnostics = {}
    for anchor in chosen:
        series, diagnostics = _load_point_anchor(anchor)
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
        "point_source_profiles": [
            profile_for_anchor(anchor) for anchor in actual_anchors
        ],
        "maximum_actual_propagation_distance": str(
            max(abs(target - anchor) for target, anchor in zip(targets, actual_anchors))
        ),
    }
    return propagated, metadata


def deterministic_segments() -> list[tuple[int, Fraction, Fraction]]:
    width = END_T - START_T
    quotient = width / ROOT_SEGMENT_WIDTH
    if quotient.denominator != 1:
        raise RuntimeError("compact root segmentation does not close")
    tasks = [
        (
            index,
            START_T + index * ROOT_SEGMENT_WIDTH,
            START_T + (index + 1) * ROOT_SEGMENT_WIDTH,
        )
        for index in range(quotient.numerator)
    ]
    if len(tasks) != EXPECTED_SEGMENTS or tasks[-1][2] != END_T:
        raise RuntimeError("compact root segment count changed")
    return tasks


def segment_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, segment_left, segment_right = task
    flint.ctx.prec = PRECISION_BITS
    propagation.MAXIMUM_PROPAGATION_DISTANCE = MAXIMUM_PROPAGATION_DISTANCE
    blocks = []
    cell_left = segment_left
    while cell_left < segment_right:
        cell_right = cell_left + CELL_WIDTH
        if cell_right > segment_right:
            raise RuntimeError(f"compact segment {index} cuts a half-cell")
        h_rows = h_rows_for_cell(cell_left, cell_right)
        specs = (
            (cell_left, cell_left, cell_left + QUARTER_WIDTH),
            (cell_right, cell_left + QUARTER_WIDTH, cell_right),
        )
        for expansion, left, right in specs:
            propagated, metadata = point_source_for(expansion, h_rows)
            block = shifted_taylor_model_curvature_row(
                expansion,
                left,
                right,
                h_rows,
                point_h_source=propagated,
            )
            blocks.append({**block, **metadata, "regime": "compact"})
        cell_left = cell_right
    if len(blocks) != int((segment_right - segment_left) / QUARTER_WIDTH):
        raise RuntimeError(f"compact segment {index} block count changed")
    if any(block.get("passed") is not True for block in blocks):
        raise RuntimeError(f"compact segment {index} contains a failing block")
    maximum = max(blocks, key=lambda block: float(block["scaled_curvature_upper"]))
    minimum = min(blocks, key=lambda block: float(block["curvature_margin_lower"]))
    return {
        "kind": "order11_compact_adaptive_h23_segment",
        "index": index,
        "segment_left": str(segment_left),
        "segment_right": str(segment_right),
        "blocks": blocks,
        "block_count": len(blocks),
        "source_profiles": sorted(
            {profile for block in blocks for profile in block["point_source_profiles"]}
        ),
        "largest_scaled_curvature_upper": maximum["scaled_curvature_upper"],
        "largest_scaled_curvature_anchor": maximum["anchor"],
        "smallest_margin_lower": minimum["curvature_margin_lower"],
        "smallest_margin_anchor": minimum["anchor"],
        "passed": True,
    }


def validate_point_source() -> dict:
    import check_jensen_window_pf_compound_order11_compact_adaptive_point_h0_h23_cache as checker

    records = point_source.load_cache(POINT_CACHE, POINT_TASKS)
    manifest = json.loads(POINT_MANIFEST.read_text(encoding="utf-8"))
    issues = checker.manifest_issues(manifest, records, POINT_TASKS, POINT_CACHE)
    if len(records) != POINT_ROWS or not manifest.get("cache", {}).get("complete"):
        issues.append(f"adaptive point source is incomplete: {len(records)}/{POINT_ROWS}")
    if issues:
        raise RuntimeError("adaptive point source failed: " + "; ".join(issues))
    return manifest


def validate_h_source() -> dict:
    import check_jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as checker

    manifest = json.loads(H_MANIFEST.read_text(encoding="utf-8"))
    cache = manifest.get("cache", {})
    expected_rows = int(h_source.DEFAULT_END_T - h_source.DEFAULT_START_T)
    if (
        manifest.get("parameters") != checker.expected_parameters()
        or cache.get("path") != relative(H_CACHE)
        or cache.get("row_count") != expected_rows
        or cache.get("all_rows_passed") is not True
        or cache.get("h_derivative_orders") != [2, 24]
        or cache.get("sha256") != sha256(H_CACHE)
    ):
        raise RuntimeError("compact H2-H24 source manifest changed")
    return manifest


def validate_boundary_source() -> dict:
    artifact = json.loads(H_BOUNDARY.read_text(encoding="utf-8"))
    _load_boundary_rows(H_BOUNDARY)
    if (
        artifact.get("kind")
        != "jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension"
        or artifact.get("summary", {}).get("tiles") != 2
        or artifact.get("summary", {}).get("failed_tiles") != 0
    ):
        raise RuntimeError("compact H boundary artifact changed")
    return artifact


def canonical_run_contract() -> dict:
    point_manifest = validate_point_source()
    h_manifest = validate_h_source()
    validate_boundary_source()
    paths = (
        POINT_CACHE,
        POINT_MANIFEST,
        H_CACHE,
        H_MANIFEST,
        H_BOUNDARY,
        Path(point_source.__file__).resolve(),
        POINT_CHECKER,
        Path(h_source.__file__).resolve(),
        H_CHECKER,
        Path(boundary_source.__file__).resolve(),
        BOUNDARY_CHECKER,
        PROPAGATION_CORE,
        TAYLOR_MODEL_CORE,
        DRIVER,
        CHECKER,
    )
    return {
        "kind": "order11_compact_adaptive_h23_run_contract",
        "date": "2026-07-20",
        "proof_boundary": (
            "This contract binds the compact first-summand quarter-block computation. "
            "It does not by itself certify any block or prove RH."
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
            "derivative_model_degrees": list(DERIVATIVE_MODEL_DEGREES),
            "stable_taylor_surplus": STABLE_TAYLOR_SURPLUS,
            "model_maximum_h_order": MODEL_MAXIMUM_H_ORDER,
            "anchor_maximum_h_order": propagation.ANCHOR_MAXIMUM_H_ORDER,
            "propagated_maximum_h_order": propagation.OUTPUT_MAXIMUM_H_ORDER,
            "remainder_h_order": propagation.REMAINDER_H_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "precision_bits": PRECISION_BITS,
            "transition_segment": [1207, "25012", "25028"],
        },
        "source_completion": {
            "point_rows": point_manifest["cache"]["row_count"],
            "point_cache_sha256": point_manifest["cache"]["sha256"],
            "lower_profile_rows": point_manifest["cache"]["lower_profile_rows"],
            "upper_profile_rows": point_manifest["cache"]["upper_profile_rows"],
            "h_rows": h_manifest["cache"]["row_count"],
            "h_cache_sha256": h_manifest["cache"]["sha256"],
            "h_boundary_tiles": 2,
        },
        "sources": [
            {"path": relative(path), "sha256": sha256(path)} for path in paths
        ],
    }


def ensure_run_contract(path: Path, canonical: dict, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != canonical:
            raise RuntimeError("compact adaptive run contract changed; use --overwrite")
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
                raise RuntimeError("compact segment cache has too many rows")
            record = json.loads(line)
            index, left, right = tasks[count]
            if (
                record.get("kind") != "order11_compact_adaptive_h23_segment"
                or record.get("index") != index
                or record.get("segment_left") != str(left)
                or record.get("segment_right") != str(right)
                or record.get("block_count") != 64
                or record.get("passed") is not True
            ):
                raise RuntimeError(
                    f"invalid compact segment row {index} at JSONL line {line_number}"
                )
            count += 1
    return count


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
    record_count = load_segment_cache(path, tasks)
    stop = len(tasks) if max_segments is None else min(len(tasks), max_segments)
    remaining = tasks[record_count:stop]
    if not remaining:
        return record_count
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    initargs = (str(POINT_CACHE), str(H_CACHE), str(H_BOUNDARY))
    executor = None
    try:
        if workers > 1:
            executor = ProcessPoolExecutor(
                max_workers=workers,
                initializer=initialize_worker,
                initargs=initargs,
            )
        else:
            initialize_worker(*initargs)
        with path.open("a", encoding="utf-8", newline="\n") as handle:
            cursor = 0
            while cursor < len(remaining):
                batch = remaining[cursor : cursor + workers]
                completed = (
                    [segment_task(task) for task in batch]
                    if executor is None
                    else list(executor.map(segment_task, batch, chunksize=1))
                )
                for record in completed:
                    handle.write(json.dumps(record, sort_keys=True) + "\n")
                    record_count += 1
                handle.flush()
                os.fsync(handle.fileno())
                cursor += len(batch)
                elapsed = perf_counter() - started
                print(
                    f"compact adaptive H23 segments: {record_count}/{stop} ({elapsed:.1f}s)",
                    flush=True,
                )
                if runtime_seconds is not None and elapsed >= runtime_seconds:
                    print(
                        "compact segment runtime checkpoint: parking after completed batch",
                        flush=True,
                    )
                    break
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return record_count


def build_pilot(indices: list[int], *, workers: int) -> dict:
    tasks = deterministic_segments()
    if not indices or len(set(indices)) != len(indices):
        raise ValueError("pilot segment indices must be nonempty and distinct")
    if any(not 0 <= index < len(tasks) for index in indices):
        raise ValueError("pilot segment index leaves the compact task list")
    canonical = canonical_run_contract()
    chosen = [tasks[index] for index in indices]
    initargs = (str(POINT_CACHE), str(H_CACHE), str(H_BOUNDARY))
    if workers > 1:
        with ProcessPoolExecutor(
            max_workers=min(workers, len(chosen)),
            initializer=initialize_worker,
            initargs=initargs,
        ) as executor:
            records = list(executor.map(segment_task, chosen, chunksize=1))
    else:
        initialize_worker(*initargs)
        records = [segment_task(task) for task in chosen]
    maximum = max(
        (block for record in records for block in record["blocks"]),
        key=lambda block: float(block["scaled_curvature_upper"]),
    )
    minimum = min(
        (block for record in records for block in record["blocks"]),
        key=lambda block: float(block["curvature_margin_lower"]),
    )
    return {
        "kind": "order11_compact_adaptive_h23_segment_pilot",
        "date": "2026-07-20",
        "status": "rigorous named compact segment pilots over hash-bound inputs",
        "proof_boundary": (
            "This certifies only the named segments. It is not contiguous compact "
            "coverage, a global order-eleven theorem, or a proof of RH."
        ),
        "run_contract": canonical,
        "run_contract_sha256": record_sha256(canonical),
        "segment_indices": indices,
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
            "wrote compact adaptive H23 pilots: "
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
    print(f"compact adaptive H23 segment rows: {count}/{len(tasks)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
