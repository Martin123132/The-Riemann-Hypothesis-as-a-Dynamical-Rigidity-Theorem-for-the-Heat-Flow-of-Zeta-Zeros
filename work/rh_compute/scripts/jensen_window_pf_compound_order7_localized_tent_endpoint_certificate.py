#!/usr/bin/env python3
"""Certify the order-seven localized-tent curvature near t=320."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
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
from jensen_window_pf_compound_order7_localized_tent_interval_core import (  # noqa: E402
    build_order7_taylor_hierarchy,
    curvature_row,
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
    "jensen_window_pf_compound_order7_localized_tent_endpoint_h2_h22_tiles.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_localized_tent_endpoint_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_localized_tent_endpoint_certificate.md"
)
MODE_START = Fraction(4619, 5000)
MODE_END = Fraction(933, 1000)
TILE_WIDTH = Fraction(1, 100_000)
TARGET_T = 320
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))


def tasks() -> list[tuple[int, Fraction, Fraction]]:
    return [
        (index, left, right)
        for index, (left, right) in enumerate(
            compact.fraction_grid(MODE_START, MODE_END, TILE_WIDTH)
        )
    ]


def initialize_worker() -> None:
    flint.ctx.prec = compact.DEFAULT_PRECISION_BITS


def tile_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = compact.DEFAULT_PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        compact.PANELS,
        window_y=compact.WINDOW_Y,
        eighth_envelope_bound=compact.EIGHTH_ENVELOPE,
        max_moment=22,
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
        "kind": "order7_localized_tent_endpoint_h2_h22_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "t_left": compact.interval_text(t_left),
        "t_right": compact.interval_text(t_right),
        "h_derivatives": {
            str(order): compact.interval_text(result["h_derivatives"][order])
            for order in range(2, 23)
        },
    }


def load_cache(path: Path, expected: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(expected):
        raise RuntimeError("endpoint cache has too many rows")
    for record, (index, left, right) in zip(records, expected):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"endpoint cache mismatch at tile {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction, Fraction]],
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
        iterator = executor.map(tile_task, remaining, chunksize=16)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"endpoint tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 1000 == 0:
                    handle.flush()
                    print(
                        f"order-seven endpoint tiles: {len(records)}/{stop} "
                        f"({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def arb_rows(records: list[dict]) -> list[dict]:
    return [
        {
            "mode_left": record["mode_left"],
            "mode_right": record["mode_right"],
            "t_left": compact.interval_from_text(record["t_left"]),
            "t_right": compact.interval_from_text(record["t_right"]),
            "H": {
                order: compact.interval_from_text(
                    record["h_derivatives"][str(order)]
                )
                for order in range(2, 23)
            },
        }
        for record in records
    ]


def build_artifact(records: list[dict]) -> dict:
    flint.ctx.prec = compact.DEFAULT_PRECISION_BITS
    hierarchy = build_order7_taylor_hierarchy(arb_rows(records))
    final_rows = hierarchy["T"]
    target_indices = [
        index
        for index, row in enumerate(final_rows)
        if bool(row["t_left"] < TARGET_T) and bool(row["t_right"] > TARGET_T)
    ]
    if len(target_indices) != 1:
        raise RuntimeError("endpoint cover does not contain a unique t=320 tile")
    final_rows = final_rows[target_indices[0] :]
    certified = []
    for row in final_rows:
        result = curvature_row(row)
        if not result["passed"]:
            raise RuntimeError(
                f"curvature failed on {row['mode_left']}: {result['scaled']}"
            )
        certified.append((row, result))
    target_rows = [
        (row, result)
        for row, result in certified
        if bool(row["t_left"] < TARGET_T < row["t_right"])
    ]
    if len(target_rows) != 1:
        raise RuntimeError("certified cover lost its t=320 tile")
    largest = max(certified, key=lambda item: item[1]["scaled"].upper())
    weakest_t = min(certified, key=lambda item: item[1]["T"].lower())
    first_row, _ = certified[0]
    last_row, _ = certified[-1]
    target_row, target_result = target_rows[0]
    return {
        "kind": "jensen_window_pf_compound_order7_localized_tent_endpoint_certificate",
        "date": "2026-07-13",
        "status": "rigorous Taylor-localized order-seven curvature theorem on a compact neighborhood of t=320",
        "proof_boundary": (
            "This artifact proves the fourth-nested first-summand curvature only "
            "on its recorded endpoint interval. It does not prove the remaining "
            "compact or saddle rays, full-kernel order seven, PF-infinity, RH, or Lambda<=0."
        ),
        "parameters": {
            "mode_start": str(MODE_START),
            "mode_end": str(MODE_END),
            "tile_width": str(TILE_WIDTH),
            "cache_rows": len(records),
            "precision_bits": compact.DEFAULT_PRECISION_BITS,
            "panels": compact.PANELS,
        },
        "stage_rows": {name: len(rows) for name, rows in hierarchy.items()},
        "certified": {
            "mode_start": first_row["mode_left"],
            "mode_end": last_row["mode_right"],
            "t_start_lower": arb_lower_text(first_row["t_left"]),
            "t_end_upper": arb_upper_text(last_row["t_right"]),
            "final_rows": len(certified),
            "all_rows_passed": True,
            "largest_scaled_upper": arb_upper_text(largest[1]["scaled"]),
            "weakest_T_lower": arb_lower_text(weakest_t[1]["T"]),
            "target_tile": {
                "mode": [target_row["mode_left"], target_row["mode_right"]],
                "t_ball": target_result["t"].str(40).replace("e", "E"),
                "J_lower": arb_lower_text(target_result["J"]),
                "R_lower": arb_lower_text(target_result["R"]),
                "S_lower": arb_lower_text(target_result["S"]),
                "T_lower": arb_lower_text(target_result["T"]),
                "scaled_upper": arb_upper_text(target_result["scaled"]),
                "margin_lower": arb_lower_text(target_result["margin"]),
            },
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_localized_tent_endpoint_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    certified = artifact["certified"]
    target = certified["target_tile"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Localized-Tent Endpoint Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous fourth-nested first-summand curvature theorem on the",
        "recorded compact neighborhood of `t=320`. This is not a proof of the",
        "remaining compact/ray ranges, full-kernel order seven, PF-infinity, RH,",
        "or `Lambda <= 0`.",
        "",
        "## Endpoint Theorem",
        "",
        "Each centered difference is evaluated through the exact symmetric",
        "Taylor inclusion `T(f) in f+Hull(f'')/12` on its unit collar.",
        "The resulting outward-rounded hierarchy proves all recorded final",
        "tiles satisfy `t^2*r_1''(t)<600`.",
        "",
        "```text",
        f"certified t start lower={certified['t_start_lower']}",
        f"certified t end upper={certified['t_end_upper']}",
        f"final rows={certified['final_rows']}",
        f"largest scaled upper={certified['largest_scaled_upper']}",
        f"weakest T lower={certified['weakest_T_lower']}",
        "```",
        "",
        "At the unique tile containing `t=320`:",
        "",
        "```text",
        f"T lower={target['T_lower']}",
        f"scaled curvature upper={target['scaled_upper']}",
        f"margin below 600={target['margin_lower']}",
        "```",
        "",
        "The remaining compact range begins at the certified right endpoint and",
        "requires a coarser continuation plus the finite and asymptotic rays.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    expected = tasks()
    records = build_cache(
        args.cache,
        expected,
        workers=max(1, args.workers),
        overwrite=args.overwrite_cache,
        max_tiles=args.max_tiles,
    )
    print(f"order-seven endpoint cache rows: {len(records)}/{len(expected)}")
    if args.cache_only:
        return 0
    if len(records) != len(expected):
        raise RuntimeError("complete the endpoint cache before certification")
    artifact = build_artifact(records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    certified = artifact["certified"]
    print(
        "wrote order-seven localized-tent endpoint certificate: "
        f"{certified['final_rows']} rows, "
        f"scaled upper {certified['largest_scaled_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
