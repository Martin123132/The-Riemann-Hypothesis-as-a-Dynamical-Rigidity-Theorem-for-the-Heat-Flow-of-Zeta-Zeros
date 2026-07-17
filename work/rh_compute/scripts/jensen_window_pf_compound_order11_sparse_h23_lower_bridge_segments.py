#!/usr/bin/env python3
"""Build order-eleven lower-bridge segments from sparse exact H0-H23 anchors."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
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

import jensen_window_pf_compound_order10_localized_lower_bridge_segments as base  # noqa: E402
import jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache as sparse_builder  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    DERIVATIVE_MODEL_DEGREES,
    MAXIMUM_H_ORDER as MODEL_MAXIMUM_H_ORDER,
    PRECISION_BITS,
    STABLE_TAYLOR_SURPLUS,
    shifted_taylor_model_curvature_row,
)
from jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core import (  # noqa: E402
    ANCHOR_MAXIMUM_H_ORDER,
    MAXIMUM_PROPAGATION_DISTANCE,
    OUTPUT_MAXIMUM_H_ORDER,
    REMAINDER_H_ORDER,
    propagated_point_source,
)


NEAR_H_CACHE = base.NEAR_H_CACHE
NEAR_H_MANIFEST = base.NEAR_H_MANIFEST
BROAD_H_CACHE = base.BROAD_H_CACHE
BROAD_H_MANIFEST = base.BROAD_H_MANIFEST
SPARSE_CACHE = sparse_builder.DEFAULT_CACHE
SPARSE_MANIFEST = sparse_builder.DEFAULT_MANIFEST
DEFAULT_SEGMENT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_sparse_h23_lower_bridge_segments.jsonl"
)
DEFAULT_RUN_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_sparse_h23_lower_bridge_run_contract.json"
)
TAYLOR_MODEL_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_shifted_taylor_model_core.py"
)
PROPAGATION_CORE = (
    SCRIPT_DIR / "jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core.py"
)
GENERATOR = (
    SCRIPT_DIR
    / "jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache.py"
)
DRIVER = Path(__file__).resolve()

START_T = Fraction(1252)
NEAR_END_T = Fraction(1500)
MIDDLE_END_T = Fraction(1800)
END_T = Fraction(5700)
CELL_WIDTH = Fraction(1, 2)
QUARTER_WIDTH = Fraction(1, 4)
ROOT_SEGMENT_WIDTH = Fraction(16)
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))


_SPARSE_RECORDS: dict[Fraction, dict] = {}
_SPARSE_ANCHORS: list[Fraction] = []


def sha256(path: Path) -> str:
    return base.sha256(path)


def record_sha256(record: dict) -> str:
    payload = json.dumps(
        record,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def _initialize_h_sources(near_h: Path, broad_h: Path) -> None:
    for kind, path, expected in (
        ("near", near_h, base.NEAR_H_ROWS),
        ("broad", broad_h, base.BROAD_H_ROWS),
    ):
        handle, view, offsets = base._open_jsonl_mmap(path)
        if len(offsets) != expected:
            raise RuntimeError(
                f"{kind} H source has {len(offsets)} rows, expected {expected}"
            )
        base._HANDLES[kind] = handle
        base._VIEWS[kind] = view
        base._OFFSETS[kind] = offsets
    base._load_h_row.cache_clear()


def load_sparse_prefix(path: Path) -> list[dict]:
    expected = sparse_builder.deterministic_tasks(
        sparse_builder.DEFAULT_START_T,
        sparse_builder.DEFAULT_END_T,
        sparse_builder.DEFAULT_STEP_T,
    )
    return sparse_builder.load_cache(path, expected)


def load_additional_pilot_records(paths: list[Path]) -> list[dict]:
    derivative_keys = {
        str(order) for order in range(ANCHOR_MAXIMUM_H_ORDER + 1)
    }
    records = []
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                record = json.loads(line)
                if (
                    record.get("passed") is not True
                    or set(record.get("h_derivatives", {})) != derivative_keys
                    or "target_t" not in record
                ):
                    raise RuntimeError(
                        f"invalid additional pilot source {path}:{line_number}"
                    )
                records.append(record)
    return records


def initialize_worker(near_h: str, broad_h: str, sparse_cache: str) -> None:
    global _SPARSE_RECORDS, _SPARSE_ANCHORS
    flint.ctx.prec = PRECISION_BITS
    _initialize_h_sources(Path(near_h), Path(broad_h))
    records = load_sparse_prefix(Path(sparse_cache))
    _SPARSE_RECORDS = {
        Fraction(record["target_t"]): record for record in records
    }
    _SPARSE_ANCHORS = sorted(_SPARSE_RECORDS)
    _load_sparse_anchor.cache_clear()


@lru_cache(maxsize=2400)
def _load_sparse_anchor(anchor: Fraction) -> tuple[list[flint.arb], dict]:
    try:
        record = _SPARSE_RECORDS[anchor]
    except KeyError as exc:
        raise RuntimeError(f"sparse exact H source misses anchor t={anchor}") from exc
    derivatives = record["h_derivatives"]
    series = [
        base.compact.interval_from_text(derivatives[str(order)])
        / math.factorial(order)
        for order in range(ANCHOR_MAXIMUM_H_ORDER + 1)
    ]
    diagnostics = {
        "target_t": str(anchor),
        "mode_bracket": [record["mode_left"], record["mode_right"]],
        "maximum_panel_error_upper": record["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": record["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
    }
    return series, diagnostics


def regime_for_cell(left: Fraction) -> str:
    if left < NEAR_END_T:
        return "near"
    if left < MIDDLE_END_T:
        return "middle"
    return "far"


def anchor_for_target(target: Fraction, regime: str) -> Fraction:
    if regime == "near":
        source_left, source_right = base.NEAR_H_START, base.NEAR_H_END
    else:
        source_left, source_right = base.BROAD_H_START, base.BROAD_H_END
    index = bisect_left(_SPARSE_ANCHORS, target)
    candidates = _SPARSE_ANCHORS[
        max(0, index - 1) : min(len(_SPARSE_ANCHORS), index + 2)
    ]
    eligible = [
        anchor
        for anchor in candidates
        if source_left <= anchor <= source_right
        and abs(anchor - target) <= MAXIMUM_PROPAGATION_DISTANCE
    ]
    if not eligible:
        raise RuntimeError(
            f"{regime} sparse H source has no admissible anchor for t={target}"
        )
    return min(eligible, key=lambda anchor: (abs(anchor - target), anchor))


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
    point_targets = [
        expansion + shift
        for expansion in (left, right)
        for shift in range(-9, 10)
    ]
    propagation_anchors = [
        anchor_for_target(target, regime) for target in point_targets
    ]
    support_left = min(left - 9, *propagation_anchors)
    support_right = max(right + 9, *propagation_anchors)
    first = base._aligned_index(support_left, start, step)
    last = base._aligned_index(support_right, start, step)
    if not 0 <= first < last <= row_count:
        raise RuntimeError(
            f"{regime} H collar {support_left}..{support_right} leaves its cache"
        )
    return [base._load_h_row(kind, index) for index in range(first, last)]


def point_source_for(
    expansion_anchor: Fraction,
    h_rows: list[dict],
    regime: str,
) -> dict[Fraction, tuple[list[flint.arb], dict]]:
    targets = [expansion_anchor + shift for shift in range(-9, 10)]
    chosen = {anchor_for_target(target, regime) for target in targets}
    anchor_source = {}
    anchor_diagnostics = {}
    for anchor in chosen:
        series, diagnostics = _load_sparse_anchor(anchor)
        anchor_source[anchor] = series
        anchor_diagnostics[anchor] = diagnostics
    return propagated_point_source(
        targets,
        anchor_source,
        h_rows,
        anchor_diagnostic_source=anchor_diagnostics,
    )


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
                point_h_source=point_source_for(expansion, h_rows, regime),
            )
            blocks.append({**block, "regime": regime})
        cell_left = cell_right
    if not blocks or any(row.get("passed") is not True for row in blocks):
        raise RuntimeError(f"segment {index} did not produce only passing blocks")
    maximum = max(blocks, key=lambda row: float(row["scaled_curvature_upper"]))
    minimum = min(blocks, key=lambda row: float(row["curvature_margin_lower"]))
    return {
        "kind": "order11_sparse_h23_lower_bridge_segment",
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


def validate_h_sources() -> None:
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


def validate_sparse_manifest() -> dict:
    manifest = json.loads(SPARSE_MANIFEST.read_text(encoding="utf-8"))
    parameters = manifest.get("parameters", {})
    cache = manifest.get("cache", {})
    expected_rows = len(
        sparse_builder.deterministic_tasks(
            sparse_builder.DEFAULT_START_T,
            sparse_builder.DEFAULT_END_T,
            sparse_builder.DEFAULT_STEP_T,
        )
    )
    if (
        parameters.get("start_t") != str(sparse_builder.DEFAULT_START_T)
        or parameters.get("end_t") != str(sparse_builder.DEFAULT_END_T)
        or parameters.get("step_t") != str(sparse_builder.DEFAULT_STEP_T)
        or parameters.get("precision_bits") != sparse_builder.PRECISION_BITS
        or parameters.get("max_moment") != sparse_builder.MAX_MOMENT
        or parameters.get("row_contract") != sparse_builder.ROW_CONTRACT
        or cache.get("row_count") != expected_rows
        or cache.get("all_rows_passed") is not True
        or cache.get("sha256") != sha256(SPARSE_CACHE)
    ):
        raise RuntimeError("order-eleven sparse H23 source manifest changed")
    return manifest


def canonical_run_contract() -> dict:
    validate_h_sources()
    validate_sparse_manifest()
    paths = (
        NEAR_H_CACHE,
        NEAR_H_MANIFEST,
        BROAD_H_CACHE,
        BROAD_H_MANIFEST,
        SPARSE_CACHE,
        SPARSE_MANIFEST,
        GENERATOR,
        PROPAGATION_CORE,
        TAYLOR_MODEL_CORE,
        DRIVER,
    )
    return {
        "kind": "order11_sparse_h23_lower_bridge_run_contract",
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
            "anchor_maximum_h_order": ANCHOR_MAXIMUM_H_ORDER,
            "propagated_maximum_h_order": OUTPUT_MAXIMUM_H_ORDER,
            "remainder_h_order": REMAINDER_H_ORDER,
            "maximum_propagation_distance": str(MAXIMUM_PROPAGATION_DISTANCE),
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
                "order-eleven sparse lower-bridge contract changed; use --overwrite"
            )
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(canonical, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_segment_cache(path: Path, tasks: list[tuple]) -> int:
    if not path.exists():
        return 0
    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            if count >= len(tasks):
                raise RuntimeError(
                    "order-eleven sparse segment cache has too many rows"
                )
            record = json.loads(line)
            index, regime, left, right = tasks[count]
            if (
                record.get("kind") != "order11_sparse_h23_lower_bridge_segment"
                or record.get("index") != index
                or record.get("regime") != regime
                or record.get("segment_left") != str(left)
                or record.get("segment_right") != str(right)
                or record.get("passed") is not True
            ):
                raise RuntimeError(
                    f"invalid order-eleven sparse segment row {index} "
                    f"at JSONL line {line_number}"
                )
            count += 1
    return count


def build_segment_cache(
    path: Path,
    tasks: list[tuple],
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
    initargs = (str(NEAR_H_CACHE), str(BROAD_H_CACHE), str(SPARSE_CACHE))
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
                    record_count += 1
                handle.flush()
                cursor += len(batch)
                elapsed = perf_counter() - started
                print(
                    "order-eleven sparse lower-bridge segments: "
                    f"{record_count}/{stop} ({elapsed:.1f}s)",
                    flush=True,
                )
                if runtime_seconds is not None and elapsed >= runtime_seconds:
                    print(
                        "order-eleven sparse segment runtime checkpoint: "
                        "parking after the completed worker batch",
                        flush=True,
                    )
                    break
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return record_count


def pilot_cell(
    cell_left: Fraction,
    sparse_cache: Path,
    additional_sparse_caches: list[Path],
) -> list[dict]:
    global _SPARSE_RECORDS, _SPARSE_ANCHORS
    if not START_T <= cell_left < END_T:
        raise ValueError("pilot cell leaves the order-eleven lower bridge")
    if ((cell_left - START_T) / CELL_WIDTH).denominator != 1:
        raise ValueError("pilot cell is not aligned to the half-grid")
    regime = regime_for_cell(cell_left)
    initialize_worker(str(NEAR_H_CACHE), str(BROAD_H_CACHE), str(sparse_cache))
    for record in load_additional_pilot_records(additional_sparse_caches):
        target = Fraction(record["target_t"])
        if target in _SPARSE_RECORDS and _SPARSE_RECORDS[target] != record:
            raise RuntimeError(f"conflicting pilot H23 sources at t={target}")
        _SPARSE_RECORDS[target] = record
    _SPARSE_ANCHORS = sorted(_SPARSE_RECORDS)
    _load_sparse_anchor.cache_clear()
    cell_right = cell_left + CELL_WIDTH
    h_rows = h_rows_for_cell(regime, cell_left, cell_right)
    blocks = []
    for expansion, left, right in (
        (cell_left, cell_left, cell_left + QUARTER_WIDTH),
        (cell_right, cell_left + QUARTER_WIDTH, cell_right),
    ):
        blocks.append(
            shifted_taylor_model_curvature_row(
                expansion,
                left,
                right,
                h_rows,
                point_h_source=point_source_for(expansion, h_rows, regime),
            )
        )
    return blocks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_SEGMENT_CACHE)
    parser.add_argument("--run-contract", type=Path, default=DEFAULT_RUN_CONTRACT)
    parser.add_argument("--sparse-cache", type=Path, default=SPARSE_CACHE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-segments", type=int)
    parser.add_argument("--runtime-seconds", type=float)
    parser.add_argument("--pilot-cell-left", type=Fraction)
    parser.add_argument("--pilot-output", type=Path)
    parser.add_argument(
        "--pilot-additional-cache",
        type=Path,
        action="append",
        default=[],
    )
    args = parser.parse_args()
    if args.runtime_seconds is not None and args.runtime_seconds <= 0:
        parser.error("--runtime-seconds must be positive")
    if args.max_segments is not None and args.max_segments < 0:
        parser.error("--max-segments cannot be negative")
    if args.pilot_cell_left is not None:
        validate_h_sources()
        rows = load_sparse_prefix(args.sparse_cache)
        blocks = pilot_cell(
            args.pilot_cell_left,
            args.sparse_cache,
            args.pilot_additional_cache,
        )
        required_targets = [
            expansion + shift
            for expansion in (
                args.pilot_cell_left,
                args.pilot_cell_left + CELL_WIDTH,
            )
            for shift in range(-9, 10)
        ]
        pilot_regime = regime_for_cell(args.pilot_cell_left)
        selected_anchors = sorted(
            {
                anchor_for_target(target, pilot_regime)
                for target in required_targets
            }
        )
        summary = {
            "kind": "order11_sparse_h23_lower_bridge_pilot",
            "cell_left": str(args.pilot_cell_left),
            "sparse_source_rows": len(rows),
            "sparse_source_path": relative(args.sparse_cache),
            "sparse_source_prefix_sha256": sha256(args.sparse_cache),
            "selected_anchor_records": [
                {
                    "target_t": str(anchor),
                    "sha256": record_sha256(_SPARSE_RECORDS[anchor]),
                }
                for anchor in selected_anchors
            ],
            "additional_pilot_sources": [
                {
                    "path": relative(path),
                    "sha256": sha256(path),
                }
                for path in args.pilot_additional_cache
            ],
            "scaled_curvature_upper": [
                block["scaled_curvature_upper"] for block in blocks
            ],
            "passed": all(block["passed"] for block in blocks),
        }
        if args.pilot_output is not None:
            artifact = {
                **summary,
                "date": "2026-07-17",
                "status": "rigorous cell-local interval pilot over hashed inputs",
                "proof_boundary": (
                    "This certifies two quarter blocks in one named half-cell. "
                    "It does not certify the full lower bridge, the continuum "
                    "curvature target, order eleven, PF-infinity, Lambda<=0, or RH."
                ),
                "h_sources": [
                    {"path": relative(path), "sha256": sha256(path)}
                    for path in (
                        NEAR_H_CACHE,
                        NEAR_H_MANIFEST,
                        BROAD_H_CACHE,
                        BROAD_H_MANIFEST,
                    )
                ],
                "code_sources": [
                    {"path": relative(path), "sha256": sha256(path)}
                    for path in (
                        GENERATOR,
                        PROPAGATION_CORE,
                        TAYLOR_MODEL_CORE,
                        DRIVER,
                    )
                ],
                "blocks": blocks,
            }
            args.pilot_output.parent.mkdir(parents=True, exist_ok=True)
            args.pilot_output.write_text(
                json.dumps(artifact, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            summary["output"] = relative(args.pilot_output)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0

    if args.pilot_output is not None:
        parser.error("--pilot-output requires --pilot-cell-left")
    if args.pilot_additional_cache:
        parser.error("additional sparse caches are allowed only in pilot mode")
    if args.sparse_cache.resolve() != SPARSE_CACHE.resolve():
        parser.error("canonical segment mode requires the canonical sparse cache")
    contract = canonical_run_contract()
    ensure_run_contract(args.run_contract, contract, overwrite=args.overwrite)
    tasks = deterministic_segments()
    record_count = build_segment_cache(
        args.cache,
        tasks,
        workers=max(1, min(4, args.workers)),
        overwrite=args.overwrite,
        max_segments=args.max_segments,
        runtime_seconds=args.runtime_seconds,
    )
    print(
        "order-eleven sparse lower-bridge segment rows: "
        f"{record_count}/{len(tasks)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
