#!/usr/bin/env python3
"""Certify r_1''(t)<=600/t^2 continuously on 320<=t<=1000."""

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
from jensen_window_pf_compound_order7_shifted_jet_continuation_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    POINT_ORDER,
    PRECISION_BITS,
    QUADRATURE_PANELS,
    QUADRATURE_TAYLOR_ORDER,
    QUADRATURE_WINDOW_Y,
    REMAINDER_ORDER,
    dimensionless_continuation_row,
)


DEFAULT_BLOCK_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md"
)
START_T = Fraction(320)
END_T = Fraction(1000)
STEP_RATIO = Fraction(161, 160)
BLOCK_GRID = Fraction(1, 10)
COLLAR_RADIUS = Fraction(5)
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))

CACHE_SPECS = (
    {
        "label": "fine_hundredth",
        "path": REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_hundredth_tiles.jsonl",
        "manifest": REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_hundredth_cache.json",
        "start": Fraction(315),
        "end": Fraction(505),
        "width": Fraction(1, 100),
    },
    {
        "label": "coarse_tenth",
        "path": REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_tenth_tiles.jsonl",
        "manifest": REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache.json",
        "start": Fraction(315),
        "end": Fraction(1005),
        "width": Fraction(1, 10),
    },
)


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


def load_source_cache(spec: dict) -> dict:
    path = spec["path"]
    manifest_path = spec["manifest"]
    if not path.exists() or not manifest_path.exists():
        raise RuntimeError(f"missing {spec['label']} H cache or manifest")
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    expected_rows = int((spec["end"] - spec["start"]) / spec["width"])
    if len(records) != expected_rows:
        raise RuntimeError(
            f"incomplete {spec['label']} H cache: {len(records)}/{expected_rows}"
        )
    for index, record in enumerate(records):
        left = spec["start"] + index * spec["width"]
        right = left + spec["width"]
        if (
            record.get("index") != index
            or record.get("target_t_left") != str(left)
            or record.get("target_t_right") != str(right)
            or record.get("passed") is not True
            or any(
                str(order) not in record.get("h_derivatives", {})
                for order in range(2, 23)
            )
        ):
            raise RuntimeError(f"{spec['label']} H cache mismatch at {index}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    digest = sha256(path)
    if (
        manifest.get("cache", {}).get("row_count") != len(records)
        or manifest.get("cache", {}).get("sha256") != digest
        or manifest.get("cache", {}).get("all_rows_passed") is not True
    ):
        raise RuntimeError(f"{spec['label']} H manifest contract failed")
    return {**spec, "records": records, "sha256": digest}


def hull_text(records: list[dict], first: int, last: int, order: int) -> list[str]:
    bounds = [record["h_derivatives"][str(order)] for record in records[first:last]]
    lower = min((bound[0] for bound in bounds), key=Decimal)
    upper = max((bound[1] for bound in bounds), key=Decimal)
    return [lower, upper]


def source_slice(source: dict, left: Fraction, right: Fraction) -> tuple[int, int]:
    first = max(0, math.floor((left - source["start"]) / source["width"]))
    last = min(
        len(source["records"]),
        math.ceil((right - source["start"]) / source["width"]),
    )
    if first >= last:
        raise RuntimeError(f"empty {source['label']} collar slice")
    first_left = Fraction(source["records"][first]["target_t_left"])
    last_right = Fraction(source["records"][last - 1]["target_t_right"])
    if first_left > left or last_right < right:
        raise RuntimeError(f"{source['label']} cache misses collar {left}..{right}")
    return first, last


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


def deterministic_tasks(sources: dict[str, dict]) -> list[dict]:
    tasks = []
    fine = sources["fine_hundredth"]
    coarse = sources["coarse_tenth"]
    for index, (anchor, right) in enumerate(geometric_blocks()):
        collar_left = anchor - COLLAR_RADIUS
        collar_right = right + COLLAR_RADIUS
        source = fine if collar_right <= fine["end"] else coarse
        first, last = source_slice(source, collar_left, collar_right)
        tasks.append(
            {
                "index": index,
                "anchor": anchor,
                "right": right,
                "source_label": source["label"],
                "source_sha256": source["sha256"],
                "source_first_index": first,
                "source_last_index": last - 1,
                "collar_left": collar_left,
                "collar_right": collar_right,
                "h_bounds": {
                    order: hull_text(source["records"], first, last, order)
                    for order in range(2, REMAINDER_ORDER + 11)
                },
            }
        )
    return tasks


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def block_task(task: dict) -> dict:
    flint.ctx.prec = PRECISION_BITS
    h_bounds = {
        order: compact.interval_from_text(bounds)
        for order, bounds in task["h_bounds"].items()
    }
    row = dimensionless_continuation_row(
        task["anchor"],
        task["right"],
        h_bounds,
    )
    row.update(
        {
            "index": task["index"],
            "source_label": task["source_label"],
            "source_sha256": task["source_sha256"],
            "source_first_index": task["source_first_index"],
            "source_last_index": task["source_last_index"],
            "source_rows": task["source_last_index"]
            - task["source_first_index"]
            + 1,
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
        raise RuntimeError("shifted-jet block cache has too many rows")
    for record, task in zip(records, tasks):
        if (
            record.get("index") != task["index"]
            or record.get("anchor") != str(task["anchor"])
            or record.get("right") != str(task["right"])
            or record.get("source_label") != task["source_label"]
            or record.get("source_sha256") != task["source_sha256"]
            or record.get("passed") is not True
        ):
            raise RuntimeError(
                f"shifted-jet block cache mismatch at {task['index']}"
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
                        f"order-seven shifted-jet bridge blocks: "
                        f"{len(records)}/{stop} ({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def source_contract(sources: dict[str, dict]) -> dict:
    return {
        source["label"]: {
            "cache": source["path"].relative_to(REPO_ROOT).as_posix(),
            "manifest": source["manifest"].relative_to(REPO_ROOT).as_posix(),
            "sha256": source["sha256"],
            "rows": len(source["records"]),
            "t_range": [str(source["start"]), str(source["end"])],
            "tile_width": str(source["width"]),
            "h_derivative_orders": [2, 22],
        }
        for source in sources.values()
    }


def build_artifact(blocks: list[dict], sources: dict[str, dict]) -> dict:
    expected = geometric_blocks()
    if len(blocks) != len(expected):
        raise RuntimeError("geometric continuation block count changed")
    for index, (block, (anchor, right)) in enumerate(zip(blocks, expected)):
        if (
            block.get("index") != index
            or block.get("anchor") != str(anchor)
            or block.get("right") != str(right)
            or block.get("passed") is not True
        ):
            raise RuntimeError(f"geometric continuation mismatch at {index}")
    for previous, current in zip(blocks, blocks[1:]):
        if previous["right"] != current["anchor"]:
            raise RuntimeError("geometric continuation blocks have a gap")

    largest = max(blocks, key=lambda row: Decimal(row["scaled_curvature_upper"]))
    weakest_margin = min(
        blocks,
        key=lambda row: Decimal(row["curvature_margin_lower"]),
    )
    weakest_t = min(
        blocks,
        key=lambda row: Decimal(
            row["common_coordinates"]["T"]["dimensionless_common_lower"]
        ),
    )
    source_counts = {
        label: sum(block["source_label"] == label for block in blocks)
        for label in sources
    }
    rows = [
        CertificateRow(
            "co7sjb_01_exact_potential_point_jets",
            "interval_input",
            "ready_to_apply",
            "Exact-potential panel Taylor quadrature encloses eleven shifted H jets at every geometric anchor.",
            "H^(0),...,H^(10) at a_i+j for -5<=j<=5",
            "First Newman summand at lambda=-100 on the recorded anchor set.",
            {
                "precision_bits": PRECISION_BITS,
                "panels": QUADRATURE_PANELS,
                "taylor_order": QUADRATURE_TAYLOR_ORDER,
                "window_y": QUADRATURE_WINDOW_Y,
            },
        ),
        CertificateRow(
            "co7sjb_02_dimensionless_remainder",
            "analytic_interval_bridge",
            "ready_to_apply",
            "A five-shift H2-H21 collar gives cancellation-preserving normalized stable jets and an order-eleven curvature remainder on every block.",
            "rho_11=a^11 r^(11)/11!; R_2<=110 |rho_11| a^-2 ((t-a)/a)^9",
            "Uses the two hashed t-uniform H caches and exact outward-rounded series algebra.",
            {
                "remainder_order": REMAINDER_ORDER,
                "relative_step": str(STEP_RATIO - 1),
                "source_blocks": source_counts,
            },
        ),
        CertificateRow(
            "co7sjb_03_stable_coordinate_domain",
            "interval_theorem",
            "ready_to_apply",
            "Every common-collar B, J, R, S, and T coordinate stays strictly positive on every continuation block.",
            "B,J,R,S,T>0 on 320<=t<=1000",
            "This establishes the real stable-log domain for the first-summand hierarchy only.",
            {
                "weakest_dimensionless_T_lower": weakest_t["common_coordinates"]["T"]["dimensionless_common_lower"],
            },
        ),
        CertificateRow(
            "co7sjb_04_continuous_curvature",
            "interval_theorem",
            "ready_to_apply",
            "The fourth-nested first-summand curvature ceiling holds continuously from t=320 through t=1000.",
            "r_1''(t)<=600/t^2 for 320<=t<=1000",
            "Finite first-summand compact theorem; later rays and the full-kernel transfer remain separate.",
            {
                "blocks": len(blocks),
                "largest_scaled_upper": largest["scaled_curvature_upper"],
                "weakest_margin_lower": weakest_margin["curvature_margin_lower"],
            },
        ),
        CertificateRow(
            "co7sjb_05_ray_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the first-summand curvature theorem beyond t=1000 and then transfer it to the full kernel.",
            "r_1''(t)<=600/t^2 for t>=1000",
            "The common compact, finite-saddle, asymptotic, and full-kernel steps remain open here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate",
        "date": "2026-07-13",
        "status": "rigorous continuous fourth-nested first-summand curvature theorem on 320<=t<=1000",
        "proof_boundary": (
            "This artifact proves r_1''(t)<=600/t^2 only for the first Newman "
            "summand at lambda=-100 on 320<=t<=1000. It does not prove the "
            "remaining ray, full-kernel order seven, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source_contract(sources),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "point_order": POINT_ORDER,
            "remainder_order": REMAINDER_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "start_t": str(START_T),
            "end_t": str(END_T),
            "step_ratio": str(STEP_RATIO),
            "block_grid": str(BLOCK_GRID),
            "collar_radius": str(COLLAR_RADIUS),
        },
        "blocks": blocks,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 4,
            "open_rows": 1,
            "blocks": len(blocks),
            "source_blocks": source_counts,
            "all_blocks_passed": True,
            "continuous_curvature_theorems": 1,
            "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
            "weakest_curvature_margin_lower": weakest_margin["curvature_margin_lower"],
            "weakest_dimensionless_T_lower": weakest_t["common_coordinates"]["T"]["dimensionless_common_lower"],
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Shifted-Jet Compact Bridge",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous continuous fourth-nested first-summand curvature",
        "theorem on `320<=t<=1000`. This is not a proof of the remaining ray,",
        "full-kernel order seven, PF-infinity, RH, or `Lambda<=0`.",
        "",
        "## Method",
        "",
        "Rational geometric blocks use `(right-anchor)/anchor<=1/160`.",
        "At each anchor, exact-potential degree-30 panel Taylor quadrature",
        "encloses all eleven shifted point jets through order ten.",
        "",
        "A five-shift `H2..H21` collar supplies normalized stable jets through",
        "order eleven. With `rho_11=a^11 r^(11)/11!`, the curvature remainder is",
        "",
        "```text",
        "110*|rho_11|/a^2*((t-a)/a)^9.",
        "```",
        "",
        "Fine `0.01` tiles cover the endpoint-sensitive region; overlapping",
        "`0.1` tiles take over before the coarse remainder becomes restrictive.",
        "Both source caches are hashed in the JSON artifact.",
        "",
        "## Theorem",
        "",
        "```text",
        "r_1''(t)<=600/t^2 for every real 320<=t<=1000",
        f"blocks={summary['blocks']}",
        f"source blocks={summary['source_blocks']}",
        f"largest scaled upper={summary['largest_scaled_curvature_upper']}",
        f"weakest curvature margin={summary['weakest_curvature_margin_lower']}",
        f"weakest dimensionless T lower={summary['weakest_dimensionless_T_lower']}",
        "```",
        "",
        "## Remaining Range",
        "",
        "The next certificate begins at `t=1000`. The common compact ray,",
        "finite saddle ray, asymptotic ray, and full-kernel transfer remain",
        "separate proof tasks.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-cache", type=Path, default=DEFAULT_BLOCK_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-blocks", action="store_true")
    parser.add_argument("--max-blocks", type=int)
    args = parser.parse_args()

    flint.ctx.prec = PRECISION_BITS
    sources = {
        spec["label"]: load_source_cache(spec) for spec in CACHE_SPECS
    }
    tasks = deterministic_tasks(sources)
    blocks = build_block_cache(
        args.block_cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite_blocks,
        max_blocks=args.max_blocks,
    )
    print(f"order-seven shifted-jet bridge cache: {len(blocks)}/{len(tasks)}")
    if len(blocks) != len(tasks):
        return 0

    artifact = build_artifact(blocks, sources)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven shifted-jet t=320..1000 certificate: "
        f"{summary['blocks']} blocks, scaled upper "
        f"{summary['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
