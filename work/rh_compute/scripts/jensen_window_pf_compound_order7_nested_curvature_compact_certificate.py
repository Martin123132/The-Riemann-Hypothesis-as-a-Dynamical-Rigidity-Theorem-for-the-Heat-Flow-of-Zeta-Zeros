#!/usr/bin/env python3
"""Certify order-seven first-summand curvature from t=1000 to mode u=2."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order7_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    evaluate_nested_curvature_from_h_cover,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_EXTENSION_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_compact_h11_h12_tiles.jsonl"
)
DEFAULT_RIGHT_COLLAR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_compact_right_collar.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.json"
)
TARGET_T = 1000
COLLAR_T = 5
INITIAL_CENTRAL_TILE_COUNT = 1000
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))


@dataclass(frozen=True)
class CompactRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order6_compact.deterministic_tasks()


def initialize_worker() -> None:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS


def extension_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4_compact.PANELS,
        window_y=order4_compact.WINDOW_Y,
        eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
        max_moment=12,
    )
    if not result.get("passed"):
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": result.get("failure"),
        }
    return {
        "kind": "order7_compact_h11_h12_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "h_derivatives": {
            str(order): order4_compact.interval_text(
                result["h_derivatives"][order]
            )
            for order in (11, 12)
        },
    }


def load_extension_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("order-seven H11/H12 cache has too many rows")
    for record, (index, left, right) in zip(records, tasks):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != {"11", "12"}
        ):
            raise RuntimeError(f"order-seven H11/H12 mismatch at tile {index}")
    return records


def build_extension_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    overwrite: bool,
    max_tiles: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_extension_cache(path, tasks)
    stop = len(tasks) if max_tiles is None else min(len(tasks), max_tiles)
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    if workers == 1:
        iterator = map(extension_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(extension_task, remaining, chunksize=16)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"H11/H12 tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 1000 == 0:
                    handle.flush()
                    print(
                        f"order-seven compact H11/H12 tiles: "
                        f"{len(records)}/{stop} ({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def build_right_collar(path: Path, *, overwrite: bool = False) -> dict:
    if path.exists() and not overwrite:
        return load_json(path)
    tasks = deterministic_tasks()
    index = len(tasks)
    left = order4_compact.OUTER_END
    right = left + order4_compact.TILE_WIDTH
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4_compact.PANELS,
        window_y=order4_compact.WINDOW_Y,
        eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
        max_moment=12,
    )
    if not result.get("passed"):
        raise RuntimeError(f"order-seven right collar failed: {result}")
    record = {
        "kind": "order7_compact_full_right_collar_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "t_left": order4_compact.interval_text(
            potential_jet_arb(arb_rational(left), 1)[1]
        ),
        "t_right": order4_compact.interval_text(
            potential_jet_arb(arb_rational(right), 1)[1]
        ),
        "h_derivatives": {
            str(order): order4_compact.interval_text(
                result["h_derivatives"][order]
            )
            for order in range(2, 13)
        },
        "normalizer_lower": result["normalizer_lower"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return record


def append_right_collar(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    collar: dict,
) -> tuple[list[dict], list[dict], list[dict]]:
    expected = len(base)
    if collar.get("index") != expected or collar.get("passed") is not True:
        raise RuntimeError("order-seven right collar is not cache-aligned")

    def derivative_row(orders: tuple[int, ...], kind: str) -> dict:
        row = {
            key: value
            for key, value in collar.items()
            if key not in {"kind", "h_derivatives"}
        }
        row["kind"] = kind
        row["h_derivatives"] = {
            str(order): collar["h_derivatives"][str(order)]
            for order in orders
        }
        return row

    return (
        [*base, derivative_row(tuple(range(2, 9)), "order4_compact_h_tile")],
        [*high, derivative_row((9, 10), "order6_compact_h9_h10_tile")],
        [*top, derivative_row((11, 12), "order7_compact_h11_h12_tile")],
    )


def first_target_tile(records: list[dict]) -> int:
    target = flint.arb(TARGET_T)
    for index, record in enumerate(records):
        left = order4_compact.interval_from_text(record["t_left"])
        right = order4_compact.interval_from_text(record["t_right"])
        if bool(right > target):
            if not bool(left < target):
                raise RuntimeError("compact target tile does not straddle t=1000")
            return index
    raise RuntimeError("compact H cache does not reach t=1000")


def collar_indices(
    records: list[dict],
    central_left: int,
    central_right: int,
) -> tuple[int, int, dict]:
    lower_t = order4_compact.interval_from_text(records[central_left]["t_left"])
    upper_t = order4_compact.interval_from_text(
        records[central_right - 1]["t_right"]
    )
    left_index = central_left
    while left_index > 0:
        candidate = order4_compact.interval_from_text(
            records[left_index]["t_left"]
        )
        if bool(candidate < lower_t - COLLAR_T):
            break
        left_index -= 1
    right_index = central_right - 1
    while right_index + 1 < len(records):
        candidate = order4_compact.interval_from_text(
            records[right_index]["t_right"]
        )
        if bool(candidate > upper_t + COLLAR_T):
            break
        right_index += 1
    outer_lower = order4_compact.interval_from_text(
        records[left_index]["t_left"]
    )
    outer_upper = order4_compact.interval_from_text(
        records[right_index]["t_right"]
    )
    if not bool(
        outer_lower < lower_t - COLLAR_T
        and outer_upper > upper_t + COLLAR_T
    ):
        raise RuntimeError("cached H tiles do not contain the required t+-5 collar")
    return left_index, right_index, {
        "left_t_collar_ball": (lower_t - outer_lower).str(40).replace("e", "E"),
        "right_t_collar_ball": (outer_upper - upper_t).str(40).replace("e", "E"),
    }


def derivative_cover(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    left: int,
    right: int,
) -> tuple[dict[int, flint.arb], dict]:
    derivatives, diagnostics = order4_compact.derivative_cover(base, left, right)
    for records, orders in ((high, (9, 10)), (top, (11, 12))):
        selected = records[left : right + 1]
        for order in orders:
            derivatives[order] = order4_compact.hull(
                [
                    order4_compact.interval_from_text(
                        row["h_derivatives"][str(order)]
                    )
                    for row in selected
                ]
            )
    diagnostics["extended_derivative_orders"] = [9, 10, 11, 12]
    return derivatives, diagnostics


def compact_block(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    central_left: int,
    central_right: int,
) -> dict:
    cover_left, cover_right, collar = collar_indices(
        base,
        central_left,
        central_right,
    )
    derivatives, diagnostics = derivative_cover(
        base,
        high,
        top,
        cover_left,
        cover_right,
    )
    diagnostics.update(collar)
    left = Fraction(base[central_left]["mode_left"])
    right = Fraction(base[central_right - 1]["mode_right"])
    result = evaluate_nested_curvature_from_h_cover(
        left,
        right,
        derivatives,
        cover_diagnostics=diagnostics,
    )
    result["central_tile_index_first"] = central_left
    result["central_tile_index_last"] = central_right - 1
    result["central_tile_count"] = central_right - central_left
    result["cover_tile_count"] = cover_right - cover_left + 1
    return result


def certify_adaptive_block(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    left: int,
    right: int,
) -> list[dict]:
    result = compact_block(base, high, top, left, right)
    if result.get("passed"):
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            f"order-seven compact failed on irreducible tile {left}: "
            f"{result.get('failure')} {result.get('detail', '')}"
        )
    midpoint = (left + right) // 2
    return certify_adaptive_block(
        base, high, top, left, midpoint
    ) + certify_adaptive_block(base, high, top, midpoint, right)


def source_contract(extension_path: Path, collar_path: Path) -> dict:
    bridge = load_json(SOURCE_BRIDGE)
    if bridge.get("summary", {}).get("all_blocks_passed") is not True:
        raise RuntimeError("shifted-jet compact bridge is not certified")
    base_hash = sha256(order4_compact.DEFAULT_CACHE)
    high_hash = sha256(order6_compact.DEFAULT_EXTENSION_CACHE)
    return {
        "base_cache": order4_compact.DEFAULT_CACHE.relative_to(REPO_ROOT).as_posix(),
        "base_cache_sha256": base_hash,
        "base_cache_rows": 107452,
        "h9_h10_cache": order6_compact.DEFAULT_EXTENSION_CACHE.relative_to(
            REPO_ROOT
        ).as_posix(),
        "h9_h10_cache_sha256": high_hash,
        "h9_h10_cache_rows": 107452,
        "h11_h12_cache": extension_path.relative_to(REPO_ROOT).as_posix(),
        "h11_h12_cache_sha256": sha256(extension_path),
        "h11_h12_cache_rows": 107452,
        "right_collar": collar_path.relative_to(REPO_ROOT).as_posix(),
        "right_collar_sha256": sha256(collar_path),
        "compact_bridge": SOURCE_BRIDGE.relative_to(REPO_ROOT).as_posix(),
        "compact_bridge_sha256": sha256(SOURCE_BRIDGE),
    }


def compact_certificate(
    base: list[dict],
    high: list[dict],
    top: list[dict],
) -> dict:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    central_left = first_target_tile(base)
    central_right = order4_compact.mode_index(Fraction(2))
    accepted = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(cursor + INITIAL_CENTRAL_TILE_COUNT, central_right)
        accepted.extend(
            certify_adaptive_block(base, high, top, cursor, endpoint)
        )
        cursor = endpoint
    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("order-seven compact cover has a mode gap")
    if accepted[-1]["central_mode"][1] != "2":
        raise RuntimeError("order-seven compact cover does not reach mode u=2")

    start_mode = Fraction(base[central_left]["mode_left"])
    start_t = potential_jet_arb(arb_rational(start_mode), 1)[1]
    first_right_t = order4_compact.interval_from_text(
        base[central_left]["t_right"]
    )
    if not bool(start_t < TARGET_T < first_right_t):
        raise RuntimeError("first compact block does not cover t=1000")

    largest = max(
        accepted,
        key=lambda row: Decimal(row["scaled_curvature_upper"]),
    )
    weakest_t = min(
        accepted,
        key=lambda row: Decimal(row["T_lower"]),
    )
    return {
        "target_t": TARGET_T,
        "mode_range": [str(start_mode), "2"],
        "start_t_ball": start_t.str(40).replace("e", "E"),
        "end_t_ball": potential_jet_arb(flint.arb(2), 1)[1]
        .str(40)
        .replace("e", "E"),
        "blocks": accepted,
        "block_count": len(accepted),
        "all_blocks_passed": True,
        "largest_scaled_curvature_upper": largest[
            "scaled_curvature_upper"
        ],
        "largest_scaled_mode": largest["central_mode"],
        "weakest_T_lower": weakest_t["T_lower"],
        "weakest_T_mode": weakest_t["central_mode"],
        "theorem": "r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2)",
    }


def build_artifact(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    extension_path: Path,
    collar_path: Path,
) -> dict:
    compact = compact_certificate(base, high, top)
    rows = [
        CompactRow(
            "co7ncc_01_aligned_h_cache",
            "interval_input",
            "ready_to_apply",
            "Aligned quadrature caches enclose H derivatives through order twelve on the compact mode range.",
            "H^(2),...,H^(12) on a strict t+-5 collar",
            "First Newman summand at lambda=-100 only.",
            {"base_rows": 107452, "extension_rows": 107452},
        ),
        CompactRow(
            "co7ncc_02_common_nested_core",
            "exact_interval_algebra",
            "ready_to_apply",
            "Four cancellation-preserving stable logarithms convert every common derivative cover into B,J,R,S,T and r'' intervals.",
            "B,J,R,S,T>0 and t^2 r_1''(t)<600",
            "Outward-rounded common-collar arithmetic; no point sampling.",
        ),
        CompactRow(
            "co7ncc_03_compact_theorem",
            "interval_theorem",
            "ready_to_apply",
            "Adaptive rational mode blocks prove the order-seven first-summand curvature ceiling from t=1000 through saddle mode two.",
            compact["theorem"],
            "Compact saddle range only.",
            {
                "blocks": compact["block_count"],
                "largest_scaled_upper": compact[
                    "largest_scaled_curvature_upper"
                ],
                "weakest_T_lower": compact["weakest_T_lower"],
            },
        ),
        CompactRow(
            "co7ncc_04_finite_ray_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the same order-seven curvature theorem on the finite saddle ray beyond mode two.",
            "r_1''(t)<=600/t^2 for 2<=u<=20",
            "The finite and asymptotic saddle rays remain separate tasks.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_nested_curvature_compact_certificate",
        "date": "2026-07-14",
        "status": "rigorous order-seven nested first-summand curvature theorem on 1000<=t<=V'(2)",
        "proof_boundary": (
            "This artifact proves the first-summand curvature only from t=1000 "
            "through saddle mode u=2. It does not prove the finite or asymptotic "
            "saddle rays, full-kernel order seven, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source_contract(extension_path, collar_path),
        "parameters": {
            "target_t": TARGET_T,
            "collar_t": COLLAR_T,
            "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
            "precision_bits": order4_compact.DEFAULT_PRECISION_BITS,
            "panels": order4_compact.PANELS,
            "curvature_constant": CURVATURE_CONSTANT,
        },
        "compact": compact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 3,
            "open_rows": 1,
            "compact_blocks": compact["block_count"],
            "compact_curvature_theorems": 1,
            "open_finite_rays": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_nested_curvature_compact_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    compact = artifact["compact"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Nested Curvature Compact Certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous first-summand order-seven curvature theorem on",
        "`1000<=t<=V'(2)`. This is not a proof of the finite or asymptotic",
        "saddle rays, full-kernel order seven, PF-infinity, RH, or `Lambda<=0`.",
        "",
        "## Common Collar",
        "",
        "The aligned compact caches enclose `H^(2),...,H^(12)` on a strict",
        "`t+-5` collar. Four outward-rounded stable-log layers preserve the",
        "common dependence and certify `B,J,R,S,T>0` before evaluating `r''`.",
        "",
        "Adaptive rational mode blocks cover the range from the tile containing",
        "`t=1000` through `u=2` without point sampling.",
        "",
        "```text",
        compact["theorem"],
        f"blocks={compact['block_count']}",
        f"largest scaled upper={compact['largest_scaled_curvature_upper']}",
        f"weakest T lower={compact['weakest_T_lower']}",
        "```",
        "",
        "## Remaining Ray",
        "",
        "The next certificate starts at saddle mode `u=2`. The finite mode ray",
        "through `u=20`, the asymptotic ray, and the full-kernel transfer remain",
        "separate proof tasks.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension-cache", type=Path, default=DEFAULT_EXTENSION_CACHE)
    parser.add_argument("--right-collar", type=Path, default=DEFAULT_RIGHT_COLLAR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks()
    extension = build_extension_cache(
        args.extension_cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite_cache,
        max_tiles=args.max_tiles,
    )
    print(f"order-seven compact H11/H12 cache: {len(extension)}/{len(tasks)}")
    if args.cache_only or len(extension) != len(tasks):
        return 0

    collar = build_right_collar(
        args.right_collar,
        overwrite=args.overwrite_cache,
    )
    base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
    high = order6_compact.load_extension_cache(
        order6_compact.DEFAULT_EXTENSION_CACHE,
        tasks,
    )
    if len(base) != len(tasks) or len(high) != len(tasks):
        raise RuntimeError("existing compact source caches are incomplete")
    base, high, extension = append_right_collar(base, high, extension, collar)
    artifact = build_artifact(
        base,
        high,
        extension,
        args.extension_cache,
        args.right_collar,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-seven nested compact certificate: "
        f"{artifact['summary']['compact_blocks']} blocks, scaled upper "
        f"{artifact['compact']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
