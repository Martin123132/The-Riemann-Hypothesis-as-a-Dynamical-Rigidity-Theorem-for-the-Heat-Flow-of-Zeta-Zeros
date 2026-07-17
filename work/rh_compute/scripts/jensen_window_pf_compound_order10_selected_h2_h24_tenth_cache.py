#!/usr/bin/env python3
"""Build selected tenth-tile H2-H24 collars for order-ten scaling tests."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache as base


REPO_ROOT = base.REPO_ROOT
DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_selected_h2_h24_tenth_tiles.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_selected_h2_h24_tenth_cache.json"
)
ANCHORS = tuple(
    Fraction(value)
    for value in (
        1300,
        1350,
        1400,
        1450,
        1500,
        1600,
        1700,
        1800,
        2000,
        3000,
        4000,
        5000,
        5500,
        5600,
        5690,
    )
)
LEFT_COLLAR = Fraction(17, 2)
RIGHT_COLLAR = Fraction(9)
TILE_WIDTH = Fraction(1, 10)
ROWS_PER_ANCHOR = 175

base.MAX_MOMENT = 24
base.TILE_WIDTH_T = TILE_WIDTH
base.CACHE_TILE_LABEL = "order10_selected_h2_h24_tenth"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple]:
    base.flint.ctx.prec = base.PRECISION_BITS
    tasks = []
    index = 0
    for anchor in ANCHORS:
        start = anchor - LEFT_COLLAR
        endpoints = [start + row * TILE_WIDTH for row in range(ROWS_PER_ANCHOR + 1)]
        brackets = [base.mode_bracket(target) for target in endpoints]
        for row in range(ROWS_PER_ANCHOR):
            tasks.append(
                (
                    index,
                    endpoints[row],
                    endpoints[row + 1],
                    brackets[row][0],
                    brackets[row + 1][1],
                )
            )
            index += 1
    return tasks


def write_manifest(path: Path, cache: Path, records: list[dict]) -> dict:
    manifest = {
        "kind": "jensen_window_pf_compound_order10_selected_h2_h24_tenth_cache",
        "date": "2026-07-16",
        "status": "rigorous selected H2-H24 interval collars; computational input only",
        "proof_boundary": (
            "These disjoint collars enclose first-summand H derivatives only. "
            "They do not prove curvature between or beyond the selected blocks."
        ),
        "parameters": {
            "anchors": [str(anchor) for anchor in ANCHORS],
            "left_collar": str(LEFT_COLLAR),
            "right_collar": str(RIGHT_COLLAR),
            "tile_width": str(TILE_WIDTH),
            "rows_per_anchor": ROWS_PER_ANCHOR,
            "max_moment": 24,
            "precision_bits": base.PRECISION_BITS,
        },
        "cache": {
            "path": cache.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [2, 24],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_selected_h2_h24_tenth_cache.py"
        ),
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
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()
    tasks = deterministic_tasks()
    records = base.build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_tiles=args.max_tiles,
    )
    print(f"order-ten selected H2-H24 rows: {len(records)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(records) != len(tasks):
        raise RuntimeError("complete the selected H2-H24 cache first")
    manifest = write_manifest(args.manifest, args.cache, records)
    print(
        "wrote selected H2-H24 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
