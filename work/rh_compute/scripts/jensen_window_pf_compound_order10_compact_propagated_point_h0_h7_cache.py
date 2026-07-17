#!/usr/bin/env python3
"""Propagate sparse exact H jets to the compact integer point grid."""

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
import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as exact_source  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    hull,
)
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    _symmetric,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_upper_text,
    upper_absolute,
)


EXACT_CACHE = exact_source.DEFAULT_CACHE
EXACT_MANIFEST = exact_source.DEFAULT_MANIFEST
H_CACHE = h_source.DEFAULT_CACHE
H_MANIFEST = h_source.DEFAULT_MANIFEST
DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_integer_grid.jsonl"
)
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache.json"
)
TARGET_START = Fraction(5692)
TARGET_END = Fraction(38026)
TARGET_STEP = Fraction(1)
REMAINDER_MOMENT = 24
OUTPUT_MAX_MOMENT = 7
DEFAULT_WORKERS = max(1, min(4, os.cpu_count() or 1))
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache.py"
)


_EXACT_COEFFICIENTS: dict[Fraction, list] = {}
_EXACT_DIAGNOSTICS: dict[Fraction, dict] = {}
_H24_MAP: dict[Fraction, flint.arb] = {}
_SOURCE_CONTRACT_ID = ""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def expected_exact_parameters() -> dict:
    return {
        "start_t": str(exact_source.START_T),
        "end_t": str(exact_source.END_T),
        "step_t": str(exact_source.STEP_T),
        "max_moment": exact_source.MAX_MOMENT,
        "window_y": exact_source.WINDOW_Y,
        "taylor_order": exact_source.TAYLOR_ORDER,
        "profile_partition": {
            "medium": [
                str(exact_source.START_T),
                str(exact_source.MEDIUM_LAST_T),
            ],
            "far": [str(exact_source.FAR_FIRST_T), str(exact_source.END_T)],
        },
        "profiles": {
            profile_name: {
                **exact_source.profiles.PROFILE_SPECS[profile_name],
                "row_contract": exact_source.row_contract(profile_name),
            }
            for profile_name in ("medium", "far")
        },
    }


def validate_sources() -> dict:
    exact_manifest = load_json(EXACT_MANIFEST)
    exact_cache = exact_manifest.get("cache", {})
    if (
        exact_manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache"
        or exact_manifest.get("parameters") != expected_exact_parameters()
        or exact_cache.get("path") != relative(EXACT_CACHE)
        or exact_cache.get("sha256") != sha256(EXACT_CACHE)
        or exact_cache.get("row_count") != 4042
        or exact_cache.get("profile_row_counts")
        != {"medium": 3039, "far": 1003}
        or exact_cache.get("all_rows_passed") is not True
        or exact_cache.get("h_derivative_orders") != [0, 23]
    ):
        raise RuntimeError("invalid sparse exact H0-H23 source")
    h_manifest = load_json(H_MANIFEST)
    h_cache = h_manifest.get("cache", {})
    if (
        h_manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_h2_h24_unit_cache"
        or h_cache.get("path") != relative(H_CACHE)
        or h_cache.get("sha256") != sha256(H_CACHE)
        or h_cache.get("row_count") != 32336
        or h_cache.get("all_rows_passed") is not True
        or h_cache.get("h_derivative_orders") != [2, 24]
    ):
        raise RuntimeError("invalid compact H2-H24 source")
    sources = [
        {
            "kind": exact_manifest["kind"],
            "path": relative(EXACT_CACHE),
            "sha256": sha256(EXACT_CACHE),
            "manifest_path": relative(EXACT_MANIFEST),
            "manifest_sha256": sha256(EXACT_MANIFEST),
            "rows": 4042,
        },
        {
            "kind": h_manifest["kind"],
            "path": relative(H_CACHE),
            "sha256": sha256(H_CACHE),
            "manifest_path": relative(H_MANIFEST),
            "manifest_sha256": sha256(H_MANIFEST),
            "rows": 32336,
        },
    ]
    formula = {
        "target_domain": [str(TARGET_START), str(TARGET_END)],
        "sparse_origin": str(exact_source.START_T),
        "sparse_step": str(exact_source.STEP_T),
        "exact_orders": [0, exact_source.MAX_MOMENT],
        "output_orders": [0, OUTPUT_MAX_MOMENT],
        "remainder_order": REMAINDER_MOMENT,
        "identity": (
            "H^(j)(a+h)/j!=sum_{m=0}^{23-j} binom(j+m,j) "
            "H^(j+m)(a)/(j+m)! h^m + R_24/j!"
        ),
    }
    fingerprint = hashlib.sha256(
        json.dumps(
            {"sources": sources, "formula": formula},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
    ).hexdigest()
    return {
        "id": fingerprint,
        "precision_load_invariant": (
            f"set flint.ctx.prec={PRECISION_BITS} before parsing serialized "
            "Arb intervals"
        ),
        "formula": formula,
        "sources": sources,
    }


def sparse_base(target: Fraction) -> Fraction:
    quotient = (target - exact_source.START_T) // exact_source.STEP_T
    return exact_source.START_T + quotient * exact_source.STEP_T


def deterministic_tasks() -> list[tuple[int, Fraction]]:
    return [
        (index, TARGET_START + index)
        for index in range(int(TARGET_END - TARGET_START) + 1)
    ]


def load_exact_source(path: Path) -> tuple[dict, dict]:
    coefficients = {}
    diagnostics = {}
    expected = exact_source.deterministic_tasks()
    with path.open("r", encoding="utf-8") as handle:
        records = [json.loads(line) for line in handle if line.strip()]
    if len(records) != len(expected):
        raise RuntimeError("sparse exact source is incomplete")
    for record, task in zip(records, expected):
        index, target, profile_name = task
        derivatives = record.get("h_derivatives", {})
        if (
            record.get("kind") != "order10_compact_sparse_point_h0_h23_jet"
            or record.get("contract_id") != exact_source.row_contract(profile_name)
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("profile") != profile_name
            or record.get("passed") is not True
            or set(derivatives) != {str(order) for order in range(24)}
        ):
            raise RuntimeError(f"invalid sparse exact source row {index}")
        coefficients[target] = [
            compact.interval_from_text(derivatives[str(order)])
            / math.factorial(order)
            for order in range(24)
        ]
        diagnostics[target] = {
            "profile": profile_name,
            "mode_left": record["mode_left"],
            "mode_right": record["mode_right"],
            "maximum_panel_error_upper": record["maximum_panel_error_upper"],
            "maximum_tail_moment_upper": record["maximum_tail_moment_upper"],
            "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
        }
    return coefficients, diagnostics


def load_h24_source(path: Path) -> dict[Fraction, flint.arb]:
    result = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            right = Fraction(record["target_t_right"])
            if (
                record.get("kind") != "order10_compact_h2_h24_unit_tile"
                or record.get("contract_id") != h_source.ROW_CONTRACT
                or record.get("passed") is not True
                or right != left + 1
            ):
                raise RuntimeError(f"invalid H24 source row at {left}")
            result[left] = compact.interval_from_text(
                record["h_derivatives"][str(REMAINDER_MOMENT)]
            )
    if set(result) != {Fraction(value) for value in range(5692, 38028)}:
        raise RuntimeError("H24 source is not the exact unit grid")
    return result


def initialize_worker(contract_id: str) -> None:
    global _EXACT_COEFFICIENTS, _EXACT_DIAGNOSTICS, _H24_MAP, _SOURCE_CONTRACT_ID
    flint.ctx.prec = PRECISION_BITS
    _EXACT_COEFFICIENTS, _EXACT_DIAGNOSTICS = load_exact_source(EXACT_CACHE)
    _H24_MAP = load_h24_source(H_CACHE)
    _SOURCE_CONTRACT_ID = contract_id


def propagate_task(task: tuple[int, Fraction]) -> dict:
    index, target = task
    base = sparse_base(target)
    displacement = int(target - base)
    base_coefficients = _EXACT_COEFFICIENTS[base]
    if displacement:
        h24_abs = upper_absolute(
            hull(
                [
                    _H24_MAP[base + offset]
                    for offset in range(displacement)
                ]
            )
        )
    else:
        h24_abs = flint.arb(0)
    coefficients = {}
    remainders = {}
    for derivative in range(OUTPUT_MAX_MOMENT + 1):
        value = flint.arb(0)
        for moment in range(exact_source.MAX_MOMENT - derivative, -1, -1):
            value = value * displacement + (
                math.comb(derivative + moment, derivative)
                * base_coefficients[derivative + moment]
            )
        remainder = (
            h24_abs
            * displacement ** (REMAINDER_MOMENT - derivative)
            / math.factorial(REMAINDER_MOMENT - derivative)
            / math.factorial(derivative)
        )
        coefficients[str(derivative)] = compact.interval_text(
            value + _symmetric(remainder)
        )
        remainders[str(derivative)] = (
            "0" if displacement == 0 else arb_upper_text(remainder)
        )
    diagnostics = _EXACT_DIAGNOSTICS[base]
    return {
        "kind": "order10_compact_propagated_point_h0_h7_jet",
        "source_contract_id": _SOURCE_CONTRACT_ID,
        "index": index,
        "target_t": str(target),
        "sparse_base_t": str(base),
        "sparse_displacement": displacement,
        "base_profile": diagnostics["profile"],
        "base_mode_left": diagnostics["mode_left"],
        "base_mode_right": diagnostics["mode_right"],
        "h_taylor_coefficients": coefficients,
        "remainder_coefficient_upper": remainders,
        "maximum_panel_error_upper": diagnostics["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": diagnostics["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": diagnostics["minimum_tail_slope_lower"],
        "passed": True,
    }


def load_cache(
    path: Path,
    expected: list[tuple[int, Fraction]],
    contract_id: str,
) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8") as handle:
        records = [json.loads(line) for line in handle if line.strip()]
    if len(records) > len(expected):
        raise RuntimeError("propagated cache has too many rows")
    keys = {str(order) for order in range(8)}
    for record, task in zip(records, expected):
        index, target = task
        base = sparse_base(target)
        if (
            record.get("kind")
            != "order10_compact_propagated_point_h0_h7_jet"
            or record.get("source_contract_id") != contract_id
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("sparse_base_t") != str(base)
            or record.get("sparse_displacement") != int(target - base)
            or set(record.get("h_taylor_coefficients", {})) != keys
            or set(record.get("remainder_coefficient_upper", {})) != keys
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"propagated cache mismatch at row {index}")
    return records


def build_cache(
    path: Path,
    expected: list[tuple[int, Fraction]],
    contract_id: str,
    *,
    workers: int,
    overwrite: bool,
    max_points: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, expected, contract_id)
    stop = len(expected) if max_points is None else min(len(expected), max_points)
    remaining = expected[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker(contract_id)
        iterator = map(propagate_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(contract_id,),
        )
        iterator = executor.map(propagate_task, remaining, chunksize=16)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 50 == 0:
                    handle.flush()
                if completed % 500 == 0:
                    elapsed = perf_counter() - started
                    rate = completed / elapsed if elapsed else 0.0
                    eta = (len(remaining) - completed) / rate if rate else math.inf
                    print(
                        "order-ten propagated point H0-H7 rows: "
                        f"{len(records)}/{stop} ({elapsed:.1f}s; ETA {eta / 60:.1f}m)",
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
    contract: dict,
) -> dict:
    maxima = {
        str(order): max(
            (record["remainder_coefficient_upper"][str(order)] for record in records),
            key=float,
        )
        for order in range(8)
    }
    manifest = {
        "kind": "jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache",
        "date": "2026-07-16",
        "status": (
            "rigorous full integer-grid H0-H7 Taylor-coefficient cache "
            "propagated from sparse exact jets"
        ),
        "proof_boundary": (
            "This cache applies exact Taylor identities with H24 interval "
            "remainders. It is a point-jet input and does not itself prove "
            "a stable-coordinate curvature theorem."
        ),
        "parameters": {
            "target_start_t": str(TARGET_START),
            "target_end_t": str(TARGET_END),
            "target_step_t": str(TARGET_STEP),
            "output_orders": [0, OUTPUT_MAX_MOMENT],
            "remainder_order": REMAINDER_MOMENT,
        },
        "source_contract": contract,
        "cache": {
            "path": relative(cache_path),
            "sha256": sha256(cache_path),
            "row_count": len(records),
            "all_rows_passed": True,
            "coefficient_orders": [0, OUTPUT_MAX_MOMENT],
            "maximum_remainder_coefficient_upper": maxima,
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
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-points", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    contract = validate_sources()
    expected = deterministic_tasks()
    records = build_cache(
        args.cache,
        expected,
        contract["id"],
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_points=args.max_points,
    )
    print(f"order-ten propagated H0-H7 cache rows: {len(records)}/{len(expected)}")
    if args.cache_only:
        return 0
    if len(records) != len(expected):
        raise RuntimeError("complete the propagated cache before promotion")
    manifest = write_manifest(args.manifest, args.cache, records, contract)
    print(
        "wrote order-ten propagated H0-H7 cache manifest: "
        f"{manifest['cache']['row_count']} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
