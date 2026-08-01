#!/usr/bin/env python3
"""Certify order-twelve first-summand curvature on 2.001<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate as order11  # noqa: E402
from jensen_window_pf_compound_order12_nested_curvature_interval_core import (  # noqa: E402
    evaluate_dimensionless_ninth_curvature,
)


RESULTS = REPO_ROOT / "work/rh_compute/results"
DEFAULT_CACHE = RESULTS / "jensen_window_pf_compound_order12_nested_curvature_u2_u20_blocks.jsonl"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order12_nested_curvature_finite_ray_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order12_nested_curvature_finite_ray_certificate.md"
SOURCE_HIGH = RESULTS / "jensen_window_pf_compound_order12_high_cumulant_coarse_corridor.json"
SOURCE_ORDER11_FINITE = RESULTS / "jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate.json"
SOURCE_ORDER11_GLOBAL = RESULTS / "jensen_window_pf_compound_order11_first_summand_curvature_certificate.json"

MODE_START = order11.MODE_START
MODE_END = order11.MODE_END
RAY_WIDTH = order11.RAY_WIDTH
PRECISION_BITS = order11.PRECISION_BITS
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))
THEOREM = "v_1''(t)<=8000/t^2 for every saddle mode 2001/1000<=u<=20"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_contract() -> dict:
    high = load_json(SOURCE_HIGH)
    inherited = load_json(SOURCE_ORDER11_FINITE)
    global_order11 = load_json(SOURCE_ORDER11_GLOBAL)
    if high.get("exact", {}).get("exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=21,22, u>=2"
    ):
        raise RuntimeError("order-twelve high-cumulant source changed")
    if (
        inherited.get("status")
        != "rigorous order-eleven first-summand curvature theorem on 2001/1000<=u<=20"
        or inherited.get("finite_ray", {}).get("all_blocks_passed") is not True
    ):
        raise RuntimeError("inherited finite derivative-corridor source changed")
    if global_order11.get("theorem") != "y_1''(t)<=6000/t^2 for every real t>=1252":
        raise RuntimeError("global order-eleven curvature source changed")
    return {
        "inherited_H2_H20": inherited["theorem"],
        "new_H21_H22": high["exact"]["exact_corridor"],
        "order11_global_curvature": global_order11["theorem"],
        "high_cumulant_sha256": sha256(SOURCE_HIGH),
        "order11_finite_sha256": sha256(SOURCE_ORDER11_FINITE),
        "order11_global_sha256": sha256(SOURCE_ORDER11_GLOBAL),
    }


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order11.ray_tasks()


def ray_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    try:
        derivatives, diagnostics = order11.corridor_h_derivatives(
            left,
            right,
            extra_absolute_caps=((21, 1), (22, 1)),
        )
        result = evaluate_dimensionless_ninth_curvature(
            left,
            right,
            derivatives,
            diagnostics=diagnostics,
        )
    except Exception as exc:
        result = {"passed": False, "failure": "exception", "detail": repr(exc)}
    return {
        "kind": "order12_nested_curvature_dimensionless_finite_ray_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def load_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
    if not path.exists():
        return []
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            index = len(records)
            if index >= len(tasks):
                raise RuntimeError("finite-ray cache has extra rows")
            expected = tasks[index]
            if (
                record.get("kind")
                != "order12_nested_curvature_dimensionless_finite_ray_block"
                or record.get("index") != expected[0]
                or record.get("mode_left") != str(expected[1])
                or record.get("mode_right") != str(expected[2])
                or record.get("passed") is not True
            ):
                raise RuntimeError(f"invalid finite-ray cache row {line_number}")
            records.append(record)
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    runtime_seconds: float,
) -> list[dict]:
    records = load_cache(path, tasks)
    path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with (
        path.open("a", encoding="utf-8") as handle,
        ProcessPoolExecutor(max_workers=workers) as pool,
    ):
        while (
            len(records) < len(tasks)
            and time.monotonic() - started < runtime_seconds
        ):
            batch = tasks[len(records) : len(records) + workers]
            results = list(pool.map(ray_task, batch))
            for record in results:
                if record.get("passed") is not True:
                    raise RuntimeError(f"finite-ray block failed: {record}")
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                handle.flush()
                records.append(record)
            if len(records) % 1000 == 0 or len(records) == len(tasks):
                print(f"order-twelve finite ray: {len(records)}/{len(tasks)}")
    return records


def build_artifact(records: list[dict], cache_path: Path = DEFAULT_CACHE) -> dict:
    tasks = ray_tasks()
    if len(records) != len(tasks):
        raise RuntimeError(f"incomplete finite-ray cache: {len(records)}/{len(tasks)}")
    largest = max(records, key=lambda row: Fraction(row["scaled_curvature_upper"]))
    weakest = min(records, key=lambda row: Fraction(row["scaled_margin_lower"]))
    weakest_d9 = min(records, key=lambda row: Fraction(row["D9_lower"]))
    return {
        "kind": "jensen_window_pf_compound_order12_nested_curvature_finite_ray_certificate",
        "date": "2026-07-22",
        "status": "rigorous order-twelve first-summand curvature theorem on 2001/1000<=u<=20",
        "theorem": THEOREM,
        "proof_boundary": (
            "This artifact proves only the finite first-summand saddle ray. The lower and "
            "compact ranges, asymptotic ray, global composition, full-kernel transfer, "
            "PF-infinity, and RH remain separate."
        ),
        "source_contract": source_contract(),
        "finite_ray": {
            "mode_range": [str(MODE_START), str(MODE_END)],
            "block_width": str(RAY_WIDTH),
            "blocks": len(records),
            "all_blocks_passed": True,
            "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
            "largest_scaled_curvature_block": largest["index"],
            "smallest_margin_lower": weakest["scaled_margin_lower"],
            "weakest_D9_lower": weakest_d9["D9_lower"],
            "cache_sha256": sha256(cache_path),
        },
        "summary": {
            "finite_ray_theorems": 1,
            "blocks": len(records),
            "failed_blocks": 0,
            "open_lower_compact_ranges": 1,
            "open_asymptotic_ranges": 1,
            "global_first_summand_theorems": 0,
            "rh_claims": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order12_nested_curvature_finite_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_nested_curvature_finite_ray_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_ray"]
    lines = [
        "# Order-Twelve Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-22",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH.",
        "",
        "```text",
        artifact["theorem"],
        f"blocks={finite['blocks']}",
        f"largest scaled upper={finite['largest_scaled_curvature_upper']}",
        f"smallest margin={finite['smallest_margin_lower']}",
        f"weakest D9 lower={finite['weakest_D9_lower']}",
        "```",
        "",
        artifact["proof_boundary"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--runtime-seconds", type=float, default=7200)
    args = parser.parse_args()
    source_contract()
    tasks = ray_tasks()
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, min(4, args.workers)),
        runtime_seconds=max(1.0, args.runtime_seconds),
    )
    if len(records) != len(tasks):
        print(f"order-twelve finite-ray checkpoint: {len(records)}/{len(tasks)}")
        return 0
    artifact = build_artifact(records, args.cache)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-twelve finite-ray theorem: "
        f"{len(records)} blocks, largest scaled upper "
        f"{artifact['finite_ray']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
