#!/usr/bin/env python3
"""Build exact integer-grid H0-H7 jets for the order-ten compact bridge."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_jet_taylor_quadrature,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_point_h0_h7_integer_grid.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_point_h0_h7_cache.json"
)
DEFAULT_START_T = Fraction(5692)
DEFAULT_END_T = Fraction(38026)
DEFAULT_STEP_T = Fraction(1)
MODE_LEFT = Fraction(3, 2)
MODE_RIGHT = Fraction(201, 100)
WINDOW_Y = 15
TAYLOR_ORDER = 30
MAX_MOMENT = 7
DEFAULT_WORKERS = max(1, min(6, os.cpu_count() or 1))
PROFILE_SPECS = {
    "base": {
        "precision_bits": 768,
        "mode_bisections": 160,
        "panels": 64,
        "row_contract": "order10-point-h0-h7-p768-b160-n64-w15-t30-v1",
    },
    "medium": {
        "precision_bits": 896,
        "mode_bisections": 180,
        "panels": 80,
        "row_contract": "order10-point-h0-h7-p896-b180-n80-w15-t30-v1",
    },
    "far": {
        "precision_bits": 1024,
        "mode_bisections": 200,
        "panels": 96,
        "row_contract": "order10-point-h0-h7-p1024-b200-n96-w15-t30-v1",
    },
}
DEFAULT_PROFILE = "far"
# Backward-compatible aliases for the conservative default contract.
PRECISION_BITS = PROFILE_SPECS[DEFAULT_PROFILE]["precision_bits"]
MODE_BISECTIONS = PROFILE_SPECS[DEFAULT_PROFILE]["mode_bisections"]
PANELS = PROFILE_SPECS[DEFAULT_PROFILE]["panels"]
ROW_CONTRACT = PROFILE_SPECS[DEFAULT_PROFILE]["row_contract"]
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order10_compact_point_h0_h7_cache.py"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mode_bracket(
    target: Fraction,
    profile_name: str = DEFAULT_PROFILE,
) -> tuple[Fraction, Fraction]:
    """Enclose the unique saddle mode whose potential derivative is target."""
    left = MODE_LEFT
    right = MODE_RIGHT
    target_arb = arb_rational(target)
    left_value = potential_jet_arb(arb_rational(left), 1)[1]
    right_value = potential_jet_arb(arb_rational(right), 1)[1]
    if not bool(left_value < target_arb < right_value):
        raise RuntimeError(f"initial point mode bracket misses t={target}")
    profile = PROFILE_SPECS[profile_name]
    for _ in range(profile["mode_bisections"]):
        midpoint = (left + right) / 2
        value = potential_jet_arb(arb_rational(midpoint), 1)[1]
        if bool(value < target_arb):
            left = midpoint
        elif bool(value > target_arb):
            right = midpoint
        else:
            raise RuntimeError(f"inconclusive point mode bracket at t={target}")
    if not bool(
        potential_jet_arb(arb_rational(left), 1)[1]
        < target_arb
        < potential_jet_arb(arb_rational(right), 1)[1]
    ):
        raise RuntimeError(f"invalid point mode bracket at t={target}")
    return left, right


def deterministic_tasks(
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
    profile_name: str = DEFAULT_PROFILE,
) -> list[tuple[int, Fraction, str]]:
    if profile_name not in PROFILE_SPECS:
        raise ValueError(f"unknown point profile: {profile_name}")
    if not 0 < start_t <= end_t or step_t <= 0:
        raise ValueError("invalid point H cache range")
    quotient = (end_t - start_t) / step_t
    if quotient.denominator != 1:
        raise ValueError("point H cache range must align with its step")
    return [
        (index, start_t + index * step_t, profile_name)
        for index in range(quotient.numerator + 1)
    ]


def initialize_worker(profile_name: str = DEFAULT_PROFILE) -> None:
    flint.ctx.prec = PROFILE_SPECS[profile_name]["precision_bits"]


def point_task(task: tuple[int, Fraction, str]) -> dict:
    index, target, profile_name = task
    profile = PROFILE_SPECS[profile_name]
    flint.ctx.prec = profile["precision_bits"]
    mode_left, mode_right = mode_bracket(target, profile_name)
    result = integrate_h_jet_taylor_quadrature(
        mode_left,
        mode_right,
        profile["panels"],
        window_y=WINDOW_Y,
        taylor_order=TAYLOR_ORDER,
        max_moment=MAX_MOMENT,
    )
    if not result.get("passed"):
        return {
            "profile": profile_name,
            "contract_id": profile["row_contract"],
            "index": index,
            "target_t": str(target),
            "mode_left": str(mode_left),
            "mode_right": str(mode_right),
            "passed": False,
            "failure": result.get("failure"),
        }
    mode_ball = (
        (arb_rational(mode_left) + arb_rational(mode_right)) / 2
        + flint.arb(
            0,
            (arb_rational(mode_right) - arb_rational(mode_left)) / 2,
        )
    )
    target_ball = potential_jet_arb(mode_ball, 1)[1]
    if not bool(target_ball.contains(arb_rational(target))):
        return {
            "profile": profile_name,
            "contract_id": profile["row_contract"],
            "index": index,
            "target_t": str(target),
            "mode_left": str(mode_left),
            "mode_right": str(mode_right),
            "passed": False,
            "failure": "target-not-in-mode-bracket",
        }
    return {
        "kind": "order10_compact_point_h0_h7_jet",
        "profile": profile_name,
        "contract_id": profile["row_contract"],
        "index": index,
        "target_t": str(target),
        "mode_left": str(mode_left),
        "mode_right": str(mode_right),
        "target_t_ball": target_ball.str(50).replace("e", "E"),
        "passed": True,
        "h_derivatives": {
            str(order): compact.interval_text(result["h_derivatives"][order])
            for order in range(MAX_MOMENT + 1)
        },
        "maximum_panel_error_upper": result["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
    }


def load_cache(
    path: Path,
    expected: list[tuple[int, Fraction, str]],
) -> list[dict]:
    if not path.exists():
        return []
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid point H JSONL row {line_number}"
                ) from exc
    if len(records) > len(expected):
        raise RuntimeError("point H cache has too many rows")
    derivative_keys = {str(order) for order in range(MAX_MOMENT + 1)}
    for record, task in zip(records, expected):
        index, target, profile_name = task
        profile = PROFILE_SPECS[profile_name]
        if (
            record.get("kind") != "order10_compact_point_h0_h7_jet"
            or record.get("profile") != profile_name
            or record.get("contract_id") != profile["row_contract"]
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != derivative_keys
        ):
            raise RuntimeError(f"point H cache mismatch at row {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction, str]],
    *,
    workers: int,
    overwrite: bool,
    max_points: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, expected)
    stop = len(expected) if max_points is None else min(len(expected), max_points)
    remaining = expected[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker(remaining[0][2])
        iterator = map(point_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(remaining[0][2],),
        )
        iterator = executor.map(point_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"point H row {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 10 == 0:
                    handle.flush()
                if completed % 25 == 0:
                    elapsed = perf_counter() - started
                    rate = completed / elapsed if elapsed else 0.0
                    remaining_seconds = (
                        (len(remaining) - completed) / rate if rate else math.inf
                    )
                    print(
                        "order-ten compact point H0-H7 rows: "
                        f"{len(records)}/{stop} ({elapsed:.1f}s; "
                        f"ETA {remaining_seconds / 3600:.2f}h)",
                        flush=True,
                    )
            handle.flush()
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def write_manifest(
    path: Path,
    cache_path: Path,
    records: list[dict],
    *,
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
    profile_name: str,
) -> dict:
    profile = PROFILE_SPECS[profile_name]
    manifest = {
        "kind": "jensen_window_pf_compound_order10_compact_point_h0_h7_cache",
        "date": "2026-07-16",
        "status": (
            "rigorous reusable exact-point H0-H7 jet cache; "
            "computational input only"
        ),
        "proof_boundary": (
            "This cache encloses first-summand H Taylor jets at exact integer "
            "targets. It does not itself prove a stable-coordinate or "
            "curvature inequality."
        ),
        "parameters": {
            "start_t": str(start_t),
            "end_t": str(end_t),
            "step_t": str(step_t),
            "initial_mode_bracket": [str(MODE_LEFT), str(MODE_RIGHT)],
            "profile": profile_name,
            "precision_bits": profile["precision_bits"],
            "mode_bisections": profile["mode_bisections"],
            "panels": profile["panels"],
            "window_y": WINDOW_Y,
            "taylor_order": TAYLOR_ORDER,
            "max_moment": MAX_MOMENT,
            "row_contract": profile["row_contract"],
        },
        "cache": {
            "path": cache_path.resolve().relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "h_derivative_orders": [0, MAX_MOMENT],
        },
        "diagnostics": {
            "maximum_panel_error_upper": max(
                (record["maximum_panel_error_upper"] for record in records),
                key=float,
            ),
            "maximum_tail_moment_upper": max(
                (record["maximum_tail_moment_upper"] for record in records),
                key=float,
            ),
        },
        "generator": GENERATOR_PATH,
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
    parser.add_argument("--start-t", type=Fraction, default=DEFAULT_START_T)
    parser.add_argument("--end-t", type=Fraction, default=DEFAULT_END_T)
    parser.add_argument("--step-t", type=Fraction, default=DEFAULT_STEP_T)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--profile",
        choices=tuple(PROFILE_SPECS),
        default=DEFAULT_PROFILE,
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-points", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks(
        args.start_t,
        args.end_t,
        args.step_t,
        args.profile,
    )
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_points=args.max_points,
    )
    print(f"order-ten point H0-H7 cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(records) != len(tasks):
        raise RuntimeError("complete the point H cache before writing its manifest")
    manifest = write_manifest(
        args.manifest,
        args.cache,
        records,
        start_t=args.start_t,
        end_t=args.end_t,
        step_t=args.step_t,
        profile_name=args.profile,
    )
    print(
        "wrote order-ten point H0-H7 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
