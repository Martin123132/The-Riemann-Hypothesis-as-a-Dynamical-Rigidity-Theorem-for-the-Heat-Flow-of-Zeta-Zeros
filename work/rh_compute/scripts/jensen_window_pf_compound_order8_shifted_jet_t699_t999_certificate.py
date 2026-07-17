#!/usr/bin/env python3
"""Certify s_1''(t)<=2000/t^2 continuously on 699<=t<=999."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate as order7_shifted  # noqa: E402
from jensen_window_pf_compound_order8_shifted_jet_continuation_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
    localized_dimensionless_continuation_row,
)


DEFAULT_BLOCK_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_shifted_jet_t699_t999_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.md"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_first_summand_curvature_bridge.json"
)
SOURCE_SPEC = order7_shifted.CACHE_SPECS[1]
START_T = Fraction(699)
END_T = Fraction(999)
STEP_RATIO = Fraction(501, 500)
BLOCK_GRID = Fraction(1, 10)
COLLAR_RADIUS = Fraction(6)
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))
WORKER_RECORDS: list[dict] | None = None


@dataclass(frozen=True)
class CertificateRow:
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


def load_source() -> dict:
    return order7_shifted.load_source_cache(SOURCE_SPEC)


def geometric_blocks() -> list[tuple[Fraction, Fraction]]:
    blocks = []
    anchor = START_T
    while anchor < END_T:
        candidate = anchor * STEP_RATIO
        snapped = math.floor(candidate / BLOCK_GRID) * BLOCK_GRID
        right = min(snapped, END_T)
        if right <= anchor:
            raise RuntimeError("geometric block grid made no progress")
        blocks.append((anchor, right))
        anchor = right
    return blocks


def deterministic_tasks(source: dict) -> list[dict]:
    tasks = []
    for index, (anchor, right) in enumerate(geometric_blocks()):
        first, last = order7_shifted.source_slice(
            source,
            anchor - COLLAR_RADIUS,
            right + COLLAR_RADIUS,
        )
        tasks.append(
            {
                "index": index,
                "anchor": anchor,
                "right": right,
                "source_sha256": source["sha256"],
                "source_first_index": first,
                "source_last_index": last - 1,
                "collar_left": anchor - COLLAR_RADIUS,
                "collar_right": right + COLLAR_RADIUS,
            }
        )
    return tasks


def initialize_worker() -> None:
    global WORKER_RECORDS
    flint.ctx.prec = PRECISION_BITS
    WORKER_RECORDS = load_source()["records"]


def block_task(task: dict) -> dict:
    flint.ctx.prec = PRECISION_BITS
    records = WORKER_RECORDS
    if records is None:
        records = load_source()["records"]
    first = task["source_first_index"]
    last = task["source_last_index"] + 1
    h_rows = []
    for record in records[first:last]:
        h_rows.append(
            {
                "target_t_left": Fraction(record["target_t_left"]),
                "target_t_right": Fraction(record["target_t_right"]),
                "H": {
                    order: compact.interval_from_text(
                        record["h_derivatives"][str(order)]
                    )
                    for order in range(2, 23)
                },
            }
        )
    row = localized_dimensionless_continuation_row(
        task["anchor"],
        task["right"],
        h_rows,
    )
    row.update(
        {
            "index": task["index"],
            "source_sha256": task["source_sha256"],
            "source_first_index": first,
            "source_last_index": last - 1,
            "source_rows": last - first,
            "collar_left": str(task["collar_left"]),
            "collar_right": str(task["collar_right"]),
        }
    )
    return row


def load_block_cache(path: Path, tasks: list[dict]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("order-eight shifted-jet cache has too many rows")
    for record, task in zip(records, tasks):
        if (
            record.get("index") != task["index"]
            or record.get("anchor") != str(task["anchor"])
            or record.get("right") != str(task["right"])
            or record.get("source_sha256") != task["source_sha256"]
            or record.get("passed") is not True
        ):
            raise RuntimeError(
                f"order-eight shifted-jet cache mismatch at {task['index']}"
            )
    return records


def build_block_cache(
    path: Path,
    tasks: list[dict],
    *,
    workers: int,
    overwrite: bool,
    max_blocks: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_block_cache(path, tasks)
    stop = len(tasks) if max_blocks is None else min(len(tasks), max_blocks)
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    if workers == 1:
        initialize_worker()
        iterator = map(block_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(block_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                handle.flush()
                records.append(record)
                if completed % 10 == 0:
                    print(
                        f"order-eight shifted-jet blocks: "
                        f"{len(records)}/{stop} ({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def build_artifact(blocks: list[dict], source: dict) -> dict:
    expected = geometric_blocks()
    if len(blocks) != len(expected):
        raise RuntimeError("order-eight geometric block count changed")
    for index, (block, (anchor, right)) in enumerate(zip(blocks, expected)):
        if (
            block.get("index") != index
            or block.get("anchor") != str(anchor)
            or block.get("right") != str(right)
            or block.get("passed") is not True
        ):
            raise RuntimeError(f"order-eight block mismatch at {index}")
    for previous, current in zip(blocks, blocks[1:]):
        if previous["right"] != current["anchor"]:
            raise RuntimeError("order-eight shifted-jet blocks have a gap")

    bridge = json.loads(SOURCE_BRIDGE.read_text(encoding="utf-8"))
    if bridge.get("exact", {}).get("continuous_target") != (
        "s_1''(t)<=4000/t^2 for every real t>=999"
    ):
        raise RuntimeError("order-eight bridge target changed")

    largest = max(blocks, key=lambda row: Decimal(row["scaled_curvature_upper"]))
    weakest_margin = min(
        blocks, key=lambda row: Decimal(row["curvature_margin_lower"])
    )
    weakest_u = min(
        blocks,
        key=lambda row: Decimal(
            row["common_coordinates"]["U"]["dimensionless_common_lower"]
        ),
    )
    rows = [
        CertificateRow(
            "co8sj_01_source_cache",
            "interval_input",
            "ready_to_apply",
            "The hashed t-uniform H2-H22 cache covers every six-unit continuation collar.",
            "H^(j)(t) enclosed for 2<=j<=22 and 693<=t<=1005",
            "Inherited first-summand derivative cache only.",
            {
                "path": source["path"].relative_to(REPO_ROOT).as_posix(),
                "sha256": source["sha256"],
                "rows": len(source["records"]),
            },
        ),
        CertificateRow(
            "co8sj_02_point_jets",
            "interval_certificate",
            "ready_to_apply",
            "Thirteen exact shifted H jets define the fifth stable point jet at every block anchor.",
            "shifts -6,...,6; point jet through order 9",
            "Finite anchor set with rigorous quadrature.",
        ),
        CertificateRow(
            "co8sj_03_localized_remainder",
            "interval_certificate",
            "ready_to_apply",
            "Five staged localized tent convolutions preserve positive stable gaps and bound the tenth-order remainder.",
            "B,J,R,S,T,U>0 on every block; localized s coefficient order 10",
            "Continuous interval remainder theorem on the displayed range.",
        ),
        CertificateRow(
            "co8sj_04_curvature",
            "interval_theorem",
            "ready_to_apply",
            "Every continuation block satisfies the order-eight first-summand curvature ceiling.",
            "s_1''(t)<=2000/t^2 for every 699<=t<=999",
            "Rigorous compact t-range theorem only.",
            {
                "largest_scaled_block": largest["index"],
                "largest_scaled_upper": largest["scaled_curvature_upper"],
                "weakest_margin_block": weakest_margin["index"],
                "weakest_margin_lower": weakest_margin["curvature_margin_lower"],
                "weakest_U_block": weakest_u["index"],
                "weakest_U_lower": weakest_u["common_coordinates"]["U"]["dimensionless_common_lower"],
            },
        ),
        CertificateRow(
            "co8sj_05_coverage",
            "exact_composition",
            "ready_to_apply",
            "The rational blocks are contiguous and cover the complete target interval.",
            "union blocks=[699,999]",
            "Coverage arithmetic only.",
        ),
        CertificateRow(
            "co8sj_06_open_ray",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the bridge curvature theorem from t=999 through the saddle ray.",
            "prove s_1''(t)<=4000/t^2 for every real t>=999",
            "The t>=999 range is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate",
        "date": "2026-07-13",
        "status": "rigorous order-eight first-summand curvature theorem on 699<=t<=999 with open t>=999 ray",
        "proof_boundary": (
            "This artifact proves s_1''(t)<=2000/t^2 only on 699<=t<=999. "
            "It does not prove the remaining ray, order-eight entry, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_ratio": str(STEP_RATIO),
            "block_grid": str(BLOCK_GRID),
            "collar_radius": str(COLLAR_RADIUS),
            "curvature_constant": CURVATURE_CONSTANT,
            "precision_bits": PRECISION_BITS,
        },
        "source": {
            "cache": source["path"].relative_to(REPO_ROOT).as_posix(),
            "manifest": source["manifest"].relative_to(REPO_ROOT).as_posix(),
            "sha256": source["sha256"],
            "rows": len(source["records"]),
            "t_range": [str(source["start"]), str(source["end"])],
            "h_derivative_orders": [2, 22],
        },
        "blocks": blocks,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 5,
            "open_rows": 1,
            "blocks": len(blocks),
            "passed_blocks": len(blocks),
            "point_shift_count": 13,
            "stable_layers": 5,
            "continuous_curvature_theorems": 1,
            "open_ray_targets": 1,
            "largest_scaled_upper": largest["scaled_curvature_upper"],
            "weakest_U_lower": weakest_u["common_coordinates"]["U"]["dimensionless_common_lower"],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Compound Order-Eight Shifted-Jet Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous continuous first-summand curvature theorem on",
        "`699<=t<=999`, with the `t>=999` ray still open. This is not a proof",
        "of order-eight entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.py",
        "```",
        "",
        "## Interval Construction",
        "",
        "The inherited hashed H2-H22 cache covers every common t+-6 collar.",
        "At each rational anchor, thirteen rigorous shifted quadratures build",
        "the point jet through order nine. Five staged localized tent",
        "convolutions preserve B,J,R,S,T,U positivity and enclose the",
        "tenth-order Taylor remainder without collapsing the stable gaps.",
        "",
        "## Theorem",
        "",
        "```text",
        f"all {summary['blocks']} rational continuation blocks pass,",
        "s_1''(t)<=2000/t^2 for every real 699<=t<=999,",
        "largest scaled upper=" + summary["largest_scaled_upper"] + "<2000,",
        "weakest dimensionless U lower=" + summary["weakest_U_lower"] + ">0.",
        "```",
        "",
        "The blocks are exactly contiguous. No interpolation between sampled",
        "points is being assumed; every interior point is covered by the",
        "Taylor remainder certificate.",
        "",
        "## Remaining Range",
        "",
        "```text",
        "prove s_1''(t)<=4000/t^2 for every real t>=999.",
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md",
        "outputs/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_BLOCK_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    source = load_source()
    tasks = deterministic_tasks(source)
    blocks = build_block_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_blocks=args.max_blocks,
    )
    print(f"order-eight shifted-jet cache rows: {len(blocks)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(blocks) != len(tasks):
        raise RuntimeError("complete the shifted-jet cache before composition")
    artifact = build_artifact(blocks, source)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "wrote order-eight shifted-jet certificate: "
        f"{artifact['summary']['blocks']} blocks, "
        f"largest scaled upper {artifact['summary']['largest_scaled_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
