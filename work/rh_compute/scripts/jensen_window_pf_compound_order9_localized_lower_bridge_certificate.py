#!/usr/bin/env python3
"""Certify the localized order-nine first-summand bridge on 1250<=t<=5700."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
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
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
from jensen_window_pf_compound_order9_localized_final_gap_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
    localized_sixth_continuation_row,
)


DEFAULT_H_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_jet_compact_h2_h22_tenth_tiles.jsonl"
)
DEFAULT_H_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_jet_compact_h2_h22_cache.json"
)
DEFAULT_POINT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
)
DEFAULT_POINT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.json"
)
DEFAULT_SEGMENT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_localized_lower_bridge_segments.jsonl"
)
DEFAULT_RUN_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_localized_lower_bridge_run_contract.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_localized_lower_bridge_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order9_localized_lower_bridge_certificate.md"
)

START_T = Fraction(1250)
END_T = Fraction(5700)
ROOT_SEGMENT_WIDTH = Fraction(16)
MINIMUM_BLOCK_WIDTH = Fraction(1, 2)
COLLAR_T = Fraction(7)
POINT_ORDER = 6
REMAINDER_ORDER = 7

H_START_T = Fraction(1243)
H_END_T = Fraction(5707)
H_STEP_T = Fraction(1, 10)
H_ROW_COUNT = 44640
H_ORDERS = tuple(range(2, 23))

POINT_START_T = Fraction(1243)
POINT_END_T = Fraction(5707)
POINT_STEP_T = Fraction(1, 2)
POINT_ROW_COUNT = 8929
POINT_ORDERS = tuple(range(9))

DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


_H_HANDLE = None
_H_MMAP = None
_H_OFFSETS: list[int] = []
_POINT_HANDLE = None
_POINT_MMAP = None
_POINT_OFFSETS: list[int] = []


@dataclass(frozen=True)
class BridgeRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_path(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def jsonl_row_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(block.count(b"\n") for block in iter(lambda: handle.read(1 << 20), b""))


def deterministic_segments() -> list[tuple[int, Fraction, Fraction]]:
    segments = []
    left = START_T
    index = 0
    while left < END_T:
        right = min(left + ROOT_SEGMENT_WIDTH, END_T)
        segments.append((index, left, right))
        left = right
        index += 1
    return segments


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


def _jsonl_record(view: mmap.mmap, offsets: list[int], index: int) -> dict:
    if not 0 <= index < len(offsets):
        raise RuntimeError(f"JSONL row index {index} is unavailable")
    end = offsets[index + 1] if index + 1 < len(offsets) else len(view)
    return json.loads(view[offsets[index] : end].strip())


def initialize_worker(h_cache: str, point_cache: str) -> None:
    global _H_HANDLE, _H_MMAP, _H_OFFSETS
    global _POINT_HANDLE, _POINT_MMAP, _POINT_OFFSETS
    flint.ctx.prec = PRECISION_BITS
    _H_HANDLE, _H_MMAP, _H_OFFSETS = _open_jsonl_mmap(Path(h_cache))
    _POINT_HANDLE, _POINT_MMAP, _POINT_OFFSETS = _open_jsonl_mmap(
        Path(point_cache)
    )
    if len(_H_OFFSETS) != H_ROW_COUNT:
        raise RuntimeError(
            f"exact-t H source has {len(_H_OFFSETS)} rows, expected {H_ROW_COUNT}"
        )
    if len(_POINT_OFFSETS) != POINT_ROW_COUNT:
        raise RuntimeError(
            "exact-point H source has "
            f"{len(_POINT_OFFSETS)} rows, expected {POINT_ROW_COUNT}"
        )


@lru_cache(maxsize=1024)
def _load_h_row(index: int) -> dict:
    if _H_MMAP is None:
        raise RuntimeError("lower-bridge worker is not initialized")
    record = _jsonl_record(_H_MMAP, _H_OFFSETS, index)
    left = H_START_T + index * H_STEP_T
    right = left + H_STEP_T
    if (
        record.get("index") != index
        or record.get("target_t_left") != str(left)
        or record.get("target_t_right") != str(right)
        or record.get("passed") is not True
        or set(record.get("h_derivatives", {}))
        != {str(order) for order in H_ORDERS}
    ):
        raise RuntimeError(f"invalid exact-t H source row {index}")
    return {
        "target_t_left": left,
        "target_t_right": right,
        "H": {
            order: compact.interval_from_text(
                record["h_derivatives"][str(order)]
            )
            for order in H_ORDERS
        },
    }


@lru_cache(maxsize=512)
def _load_point_row(index: int) -> tuple[Fraction, tuple[list, dict]]:
    if _POINT_MMAP is None:
        raise RuntimeError("lower-bridge worker is not initialized")
    record = _jsonl_record(_POINT_MMAP, _POINT_OFFSETS, index)
    target = POINT_START_T + index * POINT_STEP_T
    derivatives = record.get("h_derivatives", {})
    if (
        record.get("index") != index
        or record.get("target_t") != str(target)
        or record.get("passed") is not True
        or set(derivatives) != {str(order) for order in POINT_ORDERS}
    ):
        raise RuntimeError(f"invalid exact-point H source row {index}")
    return target, (
        [
            compact.interval_from_text(derivatives[str(order)])
            / math.factorial(order)
            for order in POINT_ORDERS
        ],
        {
            "target_t": str(target),
            "mode_bracket": [record["mode_left"], record["mode_right"]],
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


def _aligned_index(value: Fraction, start: Fraction, step: Fraction) -> int:
    quotient = (value - start) / step
    if quotient.denominator != 1:
        raise RuntimeError(f"unaligned cache target {value}")
    return quotient.numerator


def _block_sources(
    anchor: Fraction,
    right: Fraction,
) -> tuple[list[dict], dict[Fraction, tuple[list, dict]], dict]:
    support_left = anchor - COLLAR_T
    support_right = right + COLLAR_T
    first = _aligned_index(support_left, H_START_T, H_STEP_T)
    last = _aligned_index(support_right, H_START_T, H_STEP_T)
    if not 0 <= first < last <= H_ROW_COUNT:
        raise RuntimeError(
            f"H collar {support_left}..{support_right} leaves the source cache"
        )
    h_rows = [_load_h_row(index) for index in range(first, last)]

    point_source = {}
    point_indices = []
    for shift in range(-7, 8):
        target = anchor + shift
        index = _aligned_index(target, POINT_START_T, POINT_STEP_T)
        loaded_target, value = _load_point_row(index)
        if loaded_target != target:
            raise RuntimeError(f"point cache target mismatch at {target}")
        point_source[target] = value
        point_indices.append(index)
    return h_rows, point_source, {
        "h_row_index_first": first,
        "h_row_index_last": last - 1,
        "h_row_count": last - first,
        "h_t_left": str(support_left),
        "h_t_right": str(support_right),
        "point_row_index_first": min(point_indices),
        "point_row_index_last": max(point_indices),
        "point_row_count": len(point_indices),
    }


def _certify_block(anchor: Fraction, right: Fraction) -> dict:
    h_rows, point_source, source_support = _block_sources(anchor, right)
    result = localized_sixth_continuation_row(
        anchor,
        right,
        h_rows,
        point_order=POINT_ORDER,
        remainder_order=REMAINDER_ORDER,
        point_h_source=point_source,
    )
    return {**result, "source_support": source_support}


def _split_reason(exc: Exception) -> str:
    message = str(exc)
    if "localized continuation failed" in message:
        return "curvature_enclosure"
    if "nonpositive" in message:
        return "coordinate_enclosure"
    if "disjoint" in message:
        return "floor_disjoint_enclosure"
    if "incomplete" in message or "misses" in message:
        return "support_enclosure"
    return "other_interval_enclosure"


def _adaptive_blocks(
    anchor: Fraction,
    right: Fraction,
    *,
    depth: int,
    statistics: dict,
) -> list[dict]:
    statistics["attempt_count"] += 1
    statistics["maximum_depth"] = max(statistics["maximum_depth"], depth)
    try:
        return [_certify_block(anchor, right)]
    except (RuntimeError, ValueError) as exc:
        if right - anchor <= MINIMUM_BLOCK_WIDTH:
            raise RuntimeError(
                f"irreducible lower-bridge block {anchor}..{right} failed"
            ) from exc
        statistics["split_count"] += 1
        statistics["split_reasons"][_split_reason(exc)] += 1
        midpoint = (anchor + right) / 2
        if _aligned_index(midpoint, START_T, MINIMUM_BLOCK_WIDTH) < 0:
            raise RuntimeError("adaptive midpoint precedes the bridge")
        return _adaptive_blocks(
            anchor,
            midpoint,
            depth=depth + 1,
            statistics=statistics,
        ) + _adaptive_blocks(
            midpoint,
            right,
            depth=depth + 1,
            statistics=statistics,
        )


def segment_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, anchor, right = task
    flint.ctx.prec = PRECISION_BITS
    statistics = {
        "attempt_count": 0,
        "split_count": 0,
        "maximum_depth": 0,
        "split_reasons": Counter(),
    }
    blocks = _adaptive_blocks(
        anchor,
        right,
        depth=0,
        statistics=statistics,
    )
    if blocks[0]["anchor"] != str(anchor) or blocks[-1]["right"] != str(right):
        raise RuntimeError(f"adaptive segment {index} misses an endpoint")
    for previous, current in zip(blocks, blocks[1:]):
        if previous["right"] != current["anchor"]:
            raise RuntimeError(f"adaptive segment {index} has a coverage gap")
    return {
        "kind": "order9_localized_lower_bridge_segment",
        "segment_index": index,
        "root_anchor": str(anchor),
        "root_right": str(right),
        "attempt_count": statistics["attempt_count"],
        "split_count": statistics["split_count"],
        "maximum_depth": statistics["maximum_depth"],
        "split_reasons": dict(sorted(statistics["split_reasons"].items())),
        "accepted_block_count": len(blocks),
        "blocks": blocks,
        "passed": True,
    }


def _validate_block(block: dict, expected_anchor: Fraction | None = None) -> None:
    if block.get("passed") is not True:
        raise RuntimeError("lower-bridge cache contains a failed block")
    anchor = Fraction(block["anchor"])
    right = Fraction(block["right"])
    if expected_anchor is not None and anchor != expected_anchor:
        raise RuntimeError(f"lower-bridge cache gap at {expected_anchor}")
    if not anchor < right or right - anchor < MINIMUM_BLOCK_WIDTH:
        raise RuntimeError(f"invalid lower-bridge block {anchor}..{right}")
    if Decimal(block["scaled_curvature_upper"]) >= Decimal(
        CURVATURE_CONSTANT
    ):
        raise RuntimeError(f"failed curvature row {anchor}..{right}")
    if Decimal(block["curvature_margin_lower"]) <= 0:
        raise RuntimeError(f"nonpositive curvature margin {anchor}..{right}")
    for name in ("B", "J", "R", "S", "T", "U", "V"):
        if Decimal(block["point_coordinates"][name]) <= 0:
            raise RuntimeError(f"nonpositive point {name} at {anchor}")
        if Decimal(block["common_coordinates"][name]["dimensionless_common_lower"]) <= 0:
            raise RuntimeError(f"nonpositive common {name} at {anchor}")
    if block.get("quadrature", {}).get("shift_count") != 15:
        raise RuntimeError(f"bad point-quadrature support at {anchor}")


def load_segment_cache(
    path: Path,
    expected: list[tuple[int, Fraction, Fraction]],
) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(expected):
        raise RuntimeError("lower-bridge segment cache has too many rows")
    for record, task in zip(records, expected):
        index, anchor, right = task
        blocks = record.get("blocks", [])
        if (
            record.get("kind") != "order9_localized_lower_bridge_segment"
            or record.get("segment_index") != index
            or record.get("root_anchor") != str(anchor)
            or record.get("root_right") != str(right)
            or record.get("passed") is not True
            or record.get("accepted_block_count") != len(blocks)
            or not blocks
        ):
            raise RuntimeError(f"invalid lower-bridge segment row {index}")
        cursor = anchor
        for block in blocks:
            _validate_block(block, cursor)
            cursor = Fraction(block["right"])
        if cursor != right:
            raise RuntimeError(f"lower-bridge segment {index} misses its right edge")
    return records


def build_segment_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    h_cache: Path,
    point_cache: Path,
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
    start = perf_counter()
    if workers == 1:
        initialize_worker(str(h_cache), str(point_cache))
        iterator = map(segment_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(str(h_cache), str(point_cache)),
        )
        iterator = executor.map(segment_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                handle.flush()
                records.append(record)
                if completed % 5 == 0 or len(records) == stop:
                    blocks = sum(row["accepted_block_count"] for row in records)
                    print(
                        "order-nine localized lower bridge: "
                        f"{len(records)}/{stop} root segments, {blocks} accepted "
                        f"blocks ({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def canonical_run_contract(
    *,
    h_cache: Path,
    h_manifest: Path,
    point_cache: Path,
    point_manifest: Path,
) -> dict:
    sources = (
        ("exact_t_h2_h22", h_cache, h_manifest, H_ROW_COUNT),
        ("exact_point_h0_h8", point_cache, point_manifest, POINT_ROW_COUNT),
    )
    source_rows = []
    for role, cache, manifest, expected_rows in sources:
        if not cache.exists() or not manifest.exists():
            raise RuntimeError(f"lower-bridge source is missing: {cache}")
        actual_rows = jsonl_row_count(cache)
        if actual_rows != expected_rows:
            raise RuntimeError(
                f"{role} has {actual_rows} rows, expected {expected_rows}"
            )
        source_rows.append(
            {
                "role": role,
                "cache": relative_path(cache),
                "cache_sha256": sha256(cache),
                "rows": actual_rows,
                "manifest": relative_path(manifest),
                "manifest_sha256": sha256(manifest),
            }
        )
    return {
        "kind": "order9_localized_lower_bridge_run_contract",
        "date": "2026-07-14",
        "sources": source_rows,
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "root_segment_width": str(ROOT_SEGMENT_WIDTH),
            "minimum_block_width": str(MINIMUM_BLOCK_WIDTH),
            "collar_t": str(COLLAR_T),
            "point_order": POINT_ORDER,
            "remainder_order": REMAINDER_ORDER,
            "precision_bits": PRECISION_BITS,
            "curvature_constant": CURVATURE_CONSTANT,
        },
    }


def ensure_run_contract(path: Path, canonical: dict, *, overwrite: bool) -> dict:
    if path.exists() and not overwrite:
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored != canonical:
            raise RuntimeError(
                "lower-bridge source/parameter contract changed; use --overwrite "
                "only after intentionally discarding the segment cache"
            )
        return stored
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(canonical, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return canonical


def flatten_blocks(segments: list[dict]) -> list[dict]:
    return [block for segment in segments for block in segment["blocks"]]


def build_artifact(
    segments: list[dict],
    *,
    run_contract: dict,
    run_contract_path: Path,
    segment_cache_path: Path,
) -> dict:
    blocks = flatten_blocks(segments)
    if not blocks:
        raise RuntimeError("cannot build an empty lower-bridge artifact")
    cursor = START_T
    for block in blocks:
        _validate_block(block, cursor)
        cursor = Fraction(block["right"])
    if cursor != END_T:
        raise RuntimeError("lower-bridge blocks do not cover the theorem range")

    largest = max(blocks, key=lambda row: Decimal(row["scaled_curvature_upper"]))
    weakest = min(
        (
            (Decimal(block["common_coordinates"][name]["dimensionless_common_lower"]), name, block)
            for block in blocks
            for name in ("B", "J", "R", "S", "T", "U", "V")
        ),
        key=lambda item: item[0],
    )
    width_counts = Counter(block["width"] for block in blocks)
    total_attempts = sum(segment["attempt_count"] for segment in segments)
    total_splits = sum(segment["split_count"] for segment in segments)
    theorem = (
        "w_1''(t)<=4200/t^2 for every real 1250<=t<=5700"
    )
    rows = [
        BridgeRow(
            "co9llb_01_exact_sources",
            "interval_input",
            "ready_to_apply",
            "Exact-t H2-H22 collars and exact-point H0-H8 jets cover the bridge.",
            "H^(2),...,H^(22) on t+-7 and point H^(0),...,H^(8)",
            "First Newman summand at lambda=-100 only.",
        ),
        BridgeRow(
            "co9llb_02_localized_sixth_layer",
            "exact_interval_algebra",
            "ready_to_apply",
            "The sixth stable logarithm is localized with positive B,J,R,S,T,U,V.",
            "point Taylor order 6 plus localized order-7 remainder",
            "Outward-rounded arithmetic on each accepted rational block.",
        ),
        BridgeRow(
            "co9llb_03_contiguous_bridge",
            "interval_theorem",
            "ready_to_apply",
            "Adaptive dyadic blocks cover the complete lower bridge without gaps.",
            theorem,
            "The displayed finite real-t interval only.",
            {
                "root_segments": len(segments),
                "accepted_blocks": len(blocks),
                "largest_scaled_upper": largest["scaled_curvature_upper"],
            },
        ),
        BridgeRow(
            "co9llb_04_upper_composition",
            "proved_handoff",
            "ready_to_apply",
            "The separate compact and saddle-ray certificates cover every real t>=5700.",
            "compose at t=5700",
            "Those source theorems remain separate hash-bound artifacts.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_localized_lower_bridge_certificate",
        "date": "2026-07-14",
        "status": "rigorous order-nine first-summand curvature theorem on 1250<=t<=5700",
        "proof_boundary": (
            "This artifact proves only the displayed localized first-summand "
            "curvature bridge. It does not by itself prove the upper ranges, "
            "full-kernel transfer, all-shift order-nine entry, PF-infinity, RH, "
            "or Lambda<=0."
        ),
        "source_contract": {
            **run_contract,
            "run_contract": relative_path(run_contract_path),
            "run_contract_sha256": sha256(run_contract_path),
            "segment_cache": relative_path(segment_cache_path),
            "segment_cache_sha256": sha256(segment_cache_path),
            "segment_rows": len(segments),
        },
        "theorem": theorem,
        "coverage": {
            "t_range": [str(START_T), str(END_T)],
            "root_segment_count": len(segments),
            "accepted_block_count": len(blocks),
            "all_blocks_passed": True,
            "width_counts": dict(sorted(width_counts.items(), key=lambda item: Fraction(item[0]))),
            "total_attempts": total_attempts,
            "total_splits": total_splits,
            "largest_scaled_curvature_upper": largest[
                "scaled_curvature_upper"
            ],
            "largest_scaled_block": [largest["anchor"], largest["right"]],
            "weakest_common_coordinate_lower": str(weakest[0]),
            "weakest_common_coordinate_name": weakest[1],
            "weakest_common_coordinate_block": [
                weakest[2]["anchor"],
                weakest[2]["right"],
            ],
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(row.readiness != "ready_to_apply" for row in rows),
            "curvature_theorems": 1,
            "root_segments": len(segments),
            "accepted_blocks": len(blocks),
        },
        "generator": relative_path(Path(__file__)),
    }


def write_note(path: Path, artifact: dict) -> None:
    coverage = artifact["coverage"]
    widths = ", ".join(
        f"`{width}`: {count}"
        for width, count in coverage["width_counts"].items()
    )
    lines = [
        "# Order-nine localized lower-bridge certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous first-summand curvature theorem on `1250<=t<=5700`; This is not a proof of RH, `Lambda <= 0`, PF-infinity, or even the",
        "full order-nine entry theorem by itself.",
        "",
        "## Certified statement",
        "",
        "The hash-bound interval computation proves",
        "",
        f"`{artifact['theorem']}`.",
        "",
        "It uses exact-t `H2-H22` collars, exact-point `H0-H8` jets, a strict",
        "`t+-7` support collar, positive `B,J,R,S,T,U,V`, a sixth-order point",
        "Taylor polynomial, and a localized seventh-order remainder.",
        "",
        "## Coverage",
        "",
        f"- Root segments: {coverage['root_segment_count']}",
        f"- Accepted blocks: {coverage['accepted_block_count']}",
        f"- Accepted widths: {widths}",
        "- Largest scaled curvature upper: "
        + coverage["largest_scaled_curvature_upper"]
        + "<4200",
        "- Weakest common coordinate lower: "
        + coverage["weakest_common_coordinate_lower"]
        + " ("
        + coverage["weakest_common_coordinate_name"]
        + ")",
        "",
        "The adaptive algorithm starts with exact 16-unit rational segments and",
        "bisects failed interval enclosures dyadically. Every accepted block is",
        "checked independently, and the accepted endpoints join exactly from",
        "`t=1250` through `t=5700`.",
        "",
        "## Boundary",
        "",
        "Separate compact, finite-ray, and asymptotic-ray artifacts cover the",
        "first-summand range above `t=5700`. Their composition and the existing",
        "full-kernel transfer still require an explicit umbrella certificate.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h-cache", type=Path, default=DEFAULT_H_CACHE)
    parser.add_argument("--h-manifest", type=Path, default=DEFAULT_H_MANIFEST)
    parser.add_argument("--point-cache", type=Path, default=DEFAULT_POINT_CACHE)
    parser.add_argument(
        "--point-manifest", type=Path, default=DEFAULT_POINT_MANIFEST
    )
    parser.add_argument(
        "--segment-cache", type=Path, default=DEFAULT_SEGMENT_CACHE
    )
    parser.add_argument(
        "--run-contract", type=Path, default=DEFAULT_RUN_CONTRACT
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--max-segments", type=int)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    h_cache = args.h_cache.resolve()
    h_manifest = args.h_manifest.resolve()
    point_cache = args.point_cache.resolve()
    point_manifest = args.point_manifest.resolve()
    if args.overwrite and args.segment_cache.exists():
        args.segment_cache.unlink()
    canonical = canonical_run_contract(
        h_cache=h_cache,
        h_manifest=h_manifest,
        point_cache=point_cache,
        point_manifest=point_manifest,
    )
    contract = ensure_run_contract(
        args.run_contract,
        canonical,
        overwrite=args.overwrite,
    )
    tasks = deterministic_segments()
    segments = build_segment_cache(
        args.segment_cache,
        tasks,
        h_cache=h_cache,
        point_cache=point_cache,
        workers=max(1, args.workers),
        overwrite=False,
        max_segments=args.max_segments,
    )
    print(
        "order-nine lower-bridge segment cache rows: "
        f"{len(segments)}/{len(tasks)}"
    )
    if args.cache_only:
        return 0
    if len(segments) != len(tasks):
        raise RuntimeError(
            "complete the lower-bridge segment cache before writing the theorem"
        )
    artifact = build_artifact(
        segments,
        run_contract=contract,
        run_contract_path=args.run_contract,
        segment_cache_path=args.segment_cache,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-nine localized lower-bridge certificate: "
        f"{artifact['coverage']['accepted_block_count']} blocks, largest "
        f"scaled upper {artifact['coverage']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
