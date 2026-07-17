#!/usr/bin/env python3
"""Certify the localized order-four curvature ceiling on the compact mode range."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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

from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    DEFAULT_PRECISION_BITS,
    evaluate_localized_curvature_from_h_cover,
    hull,
    integrate_h_derivatives,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_localized_curvature_compact_h_tiles.jsonl"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md"
)

OUTER_START = Fraction(1851, 2000)
CENTRAL_START = Fraction(1159, 1250)
CENTRAL_END = Fraction(2)
OUTER_END = Fraction(100001, 50000)
TILE_WIDTH = Fraction(1, 100000)
INITIAL_CENTRAL_TILE_COUNT = 100
PANELS = 200
WINDOW_Y = 6
EIGHTH_ENVELOPE = Fraction(1, 50000)
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


def fraction_grid(start: Fraction, end: Fraction, width: Fraction) -> list[tuple[Fraction, Fraction]]:
    rows: list[tuple[Fraction, Fraction]] = []
    left = start
    while left < end:
        right = min(left + width, end)
        rows.append((left, right))
        left = right
    return rows


def interval_text(value: flint.arb) -> list[str]:
    return [arb_lower_text(value), arb_upper_text(value)]


def interval_from_text(bounds: list[str]) -> flint.arb:
    if len(bounds) != 2:
        raise ValueError(f"bad interval bounds: {bounds!r}")
    lower = flint.arb(bounds[0])
    upper = flint.arb(bounds[1])
    return (lower + upper) / 2 + flint.arb(0, (upper - lower) / 2)


def tile_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        PANELS,
        window_y=WINDOW_Y,
        eighth_envelope_bound=EIGHTH_ENVELOPE,
        max_moment=8,
    )
    if not result.get("passed"):
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": result.get("failure"),
        }
    t_left = potential_jet_arb(arb_rational(left), 1)[1]
    t_right = potential_jet_arb(arb_rational(right), 1)[1]
    return {
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "parameters": {
            "precision_bits": DEFAULT_PRECISION_BITS,
            "panels": PANELS,
            "window_y": WINDOW_Y,
            "eighth_envelope": str(EIGHTH_ENVELOPE),
        },
        "t_left": interval_text(t_left),
        "t_right": interval_text(t_right),
        "h_derivatives": {
            str(order): interval_text(result["h_derivatives"][order])
            for order in range(2, 9)
        },
        "normalizer_lower": result["normalizer_lower"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
    }


def load_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
    if not path.exists():
        return []
    records: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RuntimeError(f"invalid cache JSON on line {line_number}") from exc
            records.append(record)
    if len(records) > len(tasks):
        raise RuntimeError("cache has more rows than the deterministic tile grid")
    for record, (index, left, right) in zip(records, tasks):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"cache prefix mismatch at tile {index}")
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    overwrite: bool,
    max_tiles: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, tasks)
    stop = len(tasks) if max_tiles is None else min(len(tasks), max_tiles)
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records

    path.parent.mkdir(parents=True, exist_ok=True)
    start_time = perf_counter()
    if workers == 1:
        iterator = map(tile_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(tile_task, remaining, chunksize=8)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if not record.get("passed"):
                    raise RuntimeError(
                        f"H-derivative tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 1000 == 0:
                    handle.flush()
                    elapsed = perf_counter() - start_time
                    print(
                        f"localized-curvature H tiles: {len(records)}/{stop} "
                        f"({elapsed:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def cache_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mode_index(value: Fraction) -> int:
    quotient = (value - OUTER_START) / TILE_WIDTH
    if quotient.denominator != 1:
        raise RuntimeError(f"mode endpoint is not on the tile grid: {value}")
    return quotient.numerator


def derivative_cover(records: list[dict], left_index: int, right_index: int) -> tuple[dict[int, flint.arb], dict]:
    selected = records[left_index : right_index + 1]
    derivatives = {
        order: hull(
            [interval_from_text(row["h_derivatives"][str(order)]) for row in selected]
        )
        for order in range(2, 9)
    }
    diagnostics = {
        "tile_index_first": left_index,
        "tile_index_last": right_index,
        "tile_count": len(selected),
        "mode_interval": [selected[0]["mode_left"], selected[-1]["mode_right"]],
        "minimum_normalizer_lower": min(
            (row["normalizer_lower"] for row in selected), key=Decimal
        ),
        "minimum_tail_slope_lower": min(
            (row["minimum_tail_slope_lower"] for row in selected), key=Decimal
        ),
        "maximum_tail_moment_upper": max(
            (row["maximum_tail_moment_upper"] for row in selected), key=Decimal
        ),
        "maximum_simpson_error_upper": max(
            (row["maximum_simpson_error_upper"] for row in selected), key=Decimal
        ),
    }
    return derivatives, diagnostics


def cover_indices(records: list[dict], central_left: int, central_right: int) -> tuple[int, int]:
    lower_t = interval_from_text(records[central_left]["t_left"])
    upper_t = interval_from_text(records[central_right - 1]["t_right"])
    left_index = central_left
    while left_index > 0:
        candidate = interval_from_text(records[left_index]["t_left"])
        if bool(candidate < lower_t - 2):
            break
        left_index -= 1
    right_index = central_right - 1
    while right_index + 1 < len(records):
        candidate = interval_from_text(records[right_index]["t_right"])
        if bool(candidate > upper_t + 2):
            break
        right_index += 1
    outer_lower = interval_from_text(records[left_index]["t_left"])
    outer_upper = interval_from_text(records[right_index]["t_right"])
    if not bool(outer_lower < lower_t - 2 and outer_upper > upper_t + 2):
        raise RuntimeError("cached H tiles do not contain the required t collar")
    return left_index, right_index


def compact_block(records: list[dict], central_left: int, central_right: int) -> dict:
    cover_left, cover_right = cover_indices(records, central_left, central_right)
    derivatives, diagnostics = derivative_cover(records, cover_left, cover_right)
    left = Fraction(records[central_left]["mode_left"])
    right = Fraction(records[central_right - 1]["mode_right"])
    result = evaluate_localized_curvature_from_h_cover(
        left,
        right,
        derivatives,
        cover_diagnostics=diagnostics,
    )
    result["central_tile_index_first"] = central_left
    result["central_tile_index_last"] = central_right - 1
    result["central_tile_count"] = central_right - central_left
    return result


def certify_adaptive_block(records: list[dict], left: int, right: int) -> list[dict]:
    result = compact_block(records, left, right)
    if result.get("passed"):
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            f"localized curvature failed on irreducible tile {left}: "
            f"{result.get('failure')}"
        )
    midpoint = (left + right) // 2
    return certify_adaptive_block(records, left, midpoint) + certify_adaptive_block(
        records, midpoint, right
    )


def compact_certificate(
    records: list[dict], cache_path: Path, *, progress: bool = True
) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    required_count = mode_index(OUTER_END)
    if len(records) != required_count:
        raise RuntimeError(
            f"compact assembly needs {required_count} cache rows, found {len(records)}"
        )
    central_left = mode_index(CENTRAL_START)
    central_right = mode_index(CENTRAL_END)
    accepted: list[dict] = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(cursor + INITIAL_CENTRAL_TILE_COUNT, central_right)
        accepted.extend(certify_adaptive_block(records, cursor, endpoint))
        cursor = endpoint
        if progress and len(accepted) % 100 == 0:
            print(f"localized-curvature central blocks accepted: {len(accepted)}")

    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("adaptive central cover has a mode gap")
    if accepted[0]["central_mode"][0] != str(CENTRAL_START):
        raise RuntimeError("central cover has the wrong lower endpoint")
    if accepted[-1]["central_mode"][1] != str(CENTRAL_END):
        raise RuntimeError("central cover has the wrong upper endpoint")

    start_t = potential_jet_arb(arb_rational(CENTRAL_START), 1)[1]
    end_t = potential_jet_arb(arb_rational(CENTRAL_END), 1)[1]
    if not bool(start_t < 319 and end_t > 319):
        raise RuntimeError("compact mode handoff does not cover t=319")

    worst_margin_index = min(
        range(len(accepted)), key=lambda index: Decimal(accepted[index]["margin_lower"])
    )
    largest_scaled_index = max(
        range(len(accepted)),
        key=lambda index: Decimal(accepted[index]["scaled_localized_upper"]),
    )
    selected_indices = sorted(
        {0, len(accepted) // 4, len(accepted) // 2, 3 * len(accepted) // 4, len(accepted) - 1,
         worst_margin_index, largest_scaled_index}
    )

    def summary_row(index: int) -> dict:
        row = accepted[index]
        return {
            "index": index,
            "central_mode": row["central_mode"],
            "central_t_ball": row["central_t_ball"],
            "central_tile_count": row["central_tile_count"],
            "cover_tile_count": row["cover_diagnostics"]["tile_count"],
            "j_lower": row["j_lower"],
            "localized_upper": row["localized_upper"],
            "target_lower": row["target_lower"],
            "margin_lower": row["margin_lower"],
            "scaled_localized_upper": row["scaled_localized_upper"],
        }

    return {
        "kind": "jensen_window_pf_compound_order4_localized_curvature_compact_certificate",
        "date": "2026-07-12",
        "status": "rigorous compact interval certificate; analytic ray remains open",
        "proof_boundary": (
            "This artifact certifies K_1(t)<=7/(2t^2) only from t=319 through "
            "the mode u=2 handoff. It does not certify u>=2, complete order-four "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "parameters": {
            "precision_bits": DEFAULT_PRECISION_BITS,
            "outer_mode": [str(OUTER_START), str(OUTER_END)],
            "central_mode": [str(CENTRAL_START), str(CENTRAL_END)],
            "tile_width": str(TILE_WIDTH),
            "tile_count": len(records),
            "panels_per_tile": PANELS,
            "window_y": WINDOW_Y,
            "eighth_envelope": str(EIGHTH_ENVELOPE),
            "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
        },
        "mode_handoff": {
            "central_start_t_ball": start_t.str(50).replace("e", "E"),
            "central_end_t_ball": end_t.str(50).replace("e", "E"),
            "certified_t_range": "319<=t<=V'(2)",
            "open_ray": "u>=2",
        },
        "cache": {
            "path": cache_path.relative_to(REPO_ROOT).as_posix(),
            "sha256": cache_sha256(cache_path),
            "row_count": len(records),
        },
        "summary": {
            "accepted_central_blocks": len(accepted),
            "minimum_central_tile_count": min(row["central_tile_count"] for row in accepted),
            "maximum_central_tile_count": max(row["central_tile_count"] for row in accepted),
            "worst_margin_index": worst_margin_index,
            "worst_margin_lower": accepted[worst_margin_index]["margin_lower"],
            "largest_scaled_index": largest_scaled_index,
            "largest_scaled_localized_upper": accepted[largest_scaled_index]["scaled_localized_upper"],
            "all_blocks_passed": True,
        },
        "selected_blocks": [summary_row(index) for index in selected_indices],
        "blocks": [summary_row(index) for index in range(len(accepted))],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_localized_curvature_compact_certificate.py"
        ),
        "remaining_target": "Prove K_1(t)<=7/(2t^2) on the analytic mode ray u>=2.",
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    handoff = artifact["mode_handoff"]
    selected = artifact["selected_blocks"]
    lines = [
        "# Jensen-Window PF Order-Four Localized Curvature Compact Certificate",
        "",
        "Date: 2026-07-12",
        "",
        "Status: rigorous compact interval certificate with an open analytic ray.",
        "This is not a proof of complete order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_localized_curvature_compact_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json",
        "work/rh_compute/results/jensen_window_pf_compound_order4_localized_curvature_compact_h_tiles.jsonl",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.py",
        "```",
        "",
        "## Certified Range",
        "",
        "The cached Arb quadrature tiles and localized same-point inequality prove",
        "",
        "```text",
        "K_1(t)<=7/(2*t^2), 319<=t<=V'(2).",
        "```",
        "",
        f"The lower mode endpoint has `t` ball `{handoff['central_start_t_ball']}` and",
        f"the handoff endpoint has `t` ball `{handoff['central_end_t_ball']}`.",
        "The lower endpoint lies below 319, so the cover contains the full required start.",
        "",
        "## Interval Construction",
        "",
        f"The deterministic cache contains `{artifact['parameters']['tile_count']}` adjacent mode tiles of width",
        f"`{artifact['parameters']['tile_width']}`. Each tile uses `{artifact['parameters']['panels_per_tile']}` paired",
        "composite-Simpson panels on `|y|<=6`, the proved `V^(8)/a^4<=1/50000`",
        "envelope, and outward-rounded Arb arithmetic through standardized moment order eight.",
        "The expensive `H^(2),...,H^(8)` enclosures are cached once and reused.",
        "",
        f"Adaptive assembly accepted `{summary['accepted_central_blocks']}` central blocks.",
        f"The least outward-rounded curvature margin is `{summary['worst_margin_lower']}`;",
        f"the largest certified value of `t^2 U(t)` is `{summary['largest_scaled_localized_upper']}`.",
        "",
        "| mode interval | t ball | t^2 U upper | margin lower | cover tiles |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in selected:
        lines.append(
            f"| `{row['central_mode'][0]}..{row['central_mode'][1]}` | "
            f"`{row['central_t_ball']}` | `{row['scaled_localized_upper']}` | "
            f"`{row['margin_lower']}` | `{row['cover_tile_count']}` |"
        )
    lines.extend(
        [
            "",
            "## Remaining Ray",
            "",
            "The sole remaining part of the first-summand curvature theorem is",
            "",
            "```text",
            "K_1(t)<=7/(2*t^2) on the mode ray u>=2.",
            "```",
            "",
            "The ray requires analytic higher-cumulant saddle bounds. Finite scouts or",
            "extrapolation from this compact certificate are not promoted to that theorem.",
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
            "outputs/formal_core.md",
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.workers <= 0:
        raise ValueError("workers must be positive")
    tasks = [
        (index, left, right)
        for index, (left, right) in enumerate(
            fraction_grid(OUTER_START, OUTER_END, TILE_WIDTH)
        )
    ]
    records = build_cache(
        args.cache,
        tasks,
        workers=args.workers,
        overwrite=args.overwrite_cache,
        max_tiles=args.max_tiles,
    )
    print(f"localized-curvature cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only or len(records) != len(tasks):
        return 0
    artifact = compact_certificate(records, args.cache)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "certified localized order-four compact curvature: "
        f"{artifact['summary']['accepted_central_blocks']} blocks, "
        f"worst margin {artifact['summary']['worst_margin_lower']}, "
        "open ray u>=2"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
