#!/usr/bin/env python3
"""Build the aligned H15/H16 compact-mode cache for order nine."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as order8_cache  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_tiles.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_cache.json"
)
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))
PRECISION_BITS = order4_compact.DEFAULT_PRECISION_BITS


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order8_cache.deterministic_tasks()


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def tile_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4_compact.PANELS,
        window_y=order4_compact.WINDOW_Y,
        eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
        max_moment=16,
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
        "kind": "order9_compact_h15_h16_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "h_derivatives": {
            str(order): order4_compact.interval_text(
                result["h_derivatives"][order]
            )
            for order in (15, 16)
        },
    }


def load_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("order-nine H15/H16 cache has too many rows")
    for record, (index, left, right) in zip(records, tasks):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != {"15", "16"}
        ):
            raise RuntimeError(f"order-nine H15/H16 mismatch at tile {index}")
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
                        f"H15/H16 tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 1000 == 0:
                    handle.flush()
                    print(
                        "order-nine compact H15/H16 tiles: "
                        f"{len(records)}/{stop} ({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def write_manifest(path: Path, cache: Path, records: list[dict]) -> dict:
    manifest = {
        "kind": "jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_cache",
        "date": "2026-07-13",
        "status": "rigorous aligned H15/H16 compact-mode interval cache; computational input only",
        "proof_boundary": (
            "This cache encloses first-summand H derivatives 15 and 16 on the "
            "finite compact mode grid. It does not prove a curvature inequality."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "panels": order4_compact.PANELS,
            "window_y": order4_compact.WINDOW_Y,
            "h_derivative_orders": [15, 16],
        },
        "cache": {
            "path": cache.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache),
            "row_count": len(records),
            "all_rows_passed": True,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_cache.py"
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
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks()
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_tiles=args.max_tiles,
    )
    print(f"order-nine compact H15/H16 cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(records) != len(tasks):
        raise RuntimeError("complete H15/H16 cache before writing its manifest")
    manifest = write_manifest(args.manifest, args.cache, records)
    print(
        "wrote order-nine H15/H16 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
