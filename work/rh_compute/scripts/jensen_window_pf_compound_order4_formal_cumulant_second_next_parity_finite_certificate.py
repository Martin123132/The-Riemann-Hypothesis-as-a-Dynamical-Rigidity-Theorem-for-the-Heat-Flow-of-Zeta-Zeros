#!/usr/bin/env python3
"""Certify second-next cumulant coefficient bounds on 2<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
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

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    CompiledCoefficients,
    arb_interval,
    arb_rational,
    compile_coefficient_family,
    empty_extrema,
    serialize_extrema,
    sha256,
    half_power,
    integer_power,
    potential_jet_series,
    taylor_range,
    update_extrema,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_chunks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md"
)
SOURCE_SECOND_NEXT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.json"
)
MODE_START = Fraction(2)
MODE_END = Fraction(20)
BLOCK_WIDTH = Fraction(1, 200)
CHUNK_BLOCKS = 10
TAYLOR_DEGREE = 6
MAXIMUM_POTENTIAL_ORDER = 12
PRECISION_BITS = 192
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))
COEFFICIENT_BOUNDS = {
    2: (Fraction(0), Fraction(1, 20)),
    3: (Fraction(-1, 5), Fraction(0)),
    4: (Fraction(-3, 5), Fraction(0)),
    5: (Fraction(-1), Fraction(0)),
    6: (Fraction(-2), Fraction(0)),
    7: (Fraction(0), Fraction(9, 2)),
    8: (Fraction(0), Fraction(25, 4)),
}
NORMALIZED_JET_CAPS = {
    3: Fraction(6, 5),
    4: Fraction(3, 2),
    5: Fraction(2),
    6: Fraction(3),
    7: Fraction(9, 2),
    8: Fraction(7),
    9: Fraction(12),
    10: Fraction(21),
    11: Fraction(38),
    12: Fraction(71),
}
_WORKER_COEFFICIENTS: CompiledCoefficients | None = None


def compile_coefficients() -> CompiledCoefficients:
    return compile_coefficient_family(
        SOURCE_SECOND_NEXT,
        symbol_stop=13,
        expected_orders=set(COEFFICIENT_BOUNDS),
    )


def initialize_worker(coefficients: CompiledCoefficients) -> None:
    global _WORKER_COEFFICIENTS
    _WORKER_COEFFICIENTS = coefficients
    flint.ctx.prec = PRECISION_BITS
    flint.ctx.cap = TAYLOR_DEGREE + 2


def normalized_jet_series(mode: flint.arb_series) -> dict[int, flint.arb_series]:
    jet = potential_jet_series(mode, MAXIMUM_POTENTIAL_ORDER)
    q = flint.arb.pi() * (4 * mode).exp()
    curvature = jet[2]
    return {
        order: jet[order]
        * half_power(q, order - 2)
        * half_power(curvature, order).inv()
        for order in NORMALIZED_JET_CAPS
    }


def normalized_jet_taylor_range(
    left: Fraction, right: Fraction
) -> dict[int, flint.arb]:
    precision = TAYLOR_DEGREE + 2
    center = (left + right) / 2
    radius = (right - left) / 2
    point = normalized_jet_series(
        flint.arb_series([arb_rational(center), flint.arb(1)], prec=precision)
    )
    interval = normalized_jet_series(
        flint.arb_series([arb_interval(left, right), flint.arb(1)], prec=precision)
    )
    displacement = arb_interval(-radius, radius)
    remainder_power = integer_power(displacement, TAYLOR_DEGREE + 1)
    ranges = {}
    for order in NORMALIZED_JET_CAPS:
        coefficients = point[order].coeffs()
        value = flint.arb(0)
        for degree in range(TAYLOR_DEGREE, -1, -1):
            value = value * displacement + coefficients[degree]
        value += interval[order].coeffs()[TAYLOR_DEGREE + 1] * remainder_power
        ranges[order] = value
    return ranges


def chunk_task(task: tuple[int, int, int]) -> dict:
    if _WORKER_COEFFICIENTS is None:
        raise RuntimeError("second-next worker was not initialized")
    chunk_index, first_block, block_count = task
    extrema = {order: empty_extrema() for order in COEFFICIENT_BOUNDS}
    jet_extrema = {order: empty_extrema() for order in NORMALIZED_JET_CAPS}
    for offset in range(block_count):
        block_index = first_block + offset
        left = MODE_START + block_index * BLOCK_WIDTH
        right = left + BLOCK_WIDTH
        ranges = taylor_range(
            left,
            right,
            _WORKER_COEFFICIENTS,
            taylor_degree=TAYLOR_DEGREE,
            maximum_potential_order=MAXIMUM_POTENTIAL_ORDER,
        )
        jet_ranges = normalized_jet_taylor_range(left, right)
        for order, value in ranges.items():
            floor, ceiling = COEFFICIENT_BOUNDS[order]
            lower_margin = value - arb_rational(floor)
            upper_margin = arb_rational(ceiling) - value
            if not bool(lower_margin > 0 and upper_margin > 0):
                return {
                    "chunk_index": chunk_index,
                    "first_block": first_block,
                    "block_count": block_count,
                    "passed": False,
                    "failed_block": block_index,
                    "failed_order": order,
                    "value": value.str(40).replace("e", "E"),
                    "lower_margin": lower_margin.str(40).replace("e", "E"),
                    "upper_margin": upper_margin.str(40).replace("e", "E"),
                }
            update_extrema(
                extrema[order], value, lower_margin, upper_margin, block_index
            )
        for order, value in jet_ranges.items():
            lower_margin = value
            upper_margin = arb_rational(NORMALIZED_JET_CAPS[order]) - value
            if not bool(lower_margin > 0 and upper_margin > 0):
                return {
                    "chunk_index": chunk_index,
                    "first_block": first_block,
                    "block_count": block_count,
                    "passed": False,
                    "failed_block": block_index,
                    "failed_normalized_jet": order,
                    "value": value.str(40).replace("e", "E"),
                    "lower_margin": lower_margin.str(40).replace("e", "E"),
                    "upper_margin": upper_margin.str(40).replace("e", "E"),
                }
            update_extrema(
                jet_extrema[order], value, lower_margin, upper_margin, block_index
            )
    left = MODE_START + first_block * BLOCK_WIDTH
    right = left + block_count * BLOCK_WIDTH
    return {
        "chunk_index": chunk_index,
        "first_block": first_block,
        "block_count": block_count,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "block_width": str(BLOCK_WIDTH),
            "taylor_degree": TAYLOR_DEGREE,
            "maximum_potential_order": MAXIMUM_POTENTIAL_ORDER,
        },
        "coefficients": {
            str(order): serialize_extrema(values)
            for order, values in extrema.items()
        },
        "normalized_jets": {
            str(order): serialize_extrema(values)
            for order, values in jet_extrema.items()
        },
    }


def deterministic_tasks() -> list[tuple[int, int, int]]:
    total_fraction = (MODE_END - MODE_START) / BLOCK_WIDTH
    if total_fraction.denominator != 1:
        raise RuntimeError("second-next finite range is not aligned to the grid")
    total_blocks = total_fraction.numerator
    tasks = []
    first = 0
    chunk = 0
    while first < total_blocks:
        count = min(CHUNK_BLOCKS, total_blocks - first)
        tasks.append((chunk, first, count))
        first += count
        chunk += 1
    return tasks


def load_cache(path: Path, tasks: list[tuple[int, int, int]]) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("second-next finite cache has too many rows")
    for record, (chunk, first, count) in zip(records, tasks):
        parameters = record.get("parameters", {})
        if (
            record.get("chunk_index") != chunk
            or record.get("first_block") != first
            or record.get("block_count") != count
            or record.get("passed") is not True
            or parameters.get("taylor_degree") != TAYLOR_DEGREE
            or parameters.get("maximum_potential_order") != MAXIMUM_POTENTIAL_ORDER
            or set(record.get("normalized_jets", {}))
            != {str(order) for order in NORMALIZED_JET_CAPS}
        ):
            raise RuntimeError(f"second-next finite cache mismatch at chunk {chunk}")
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, int, int]],
    *,
    workers: int,
    overwrite: bool,
    max_chunks: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, tasks)
    stop = len(tasks) if max_chunks is None else min(max_chunks, len(tasks))
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records

    coefficients = compile_coefficients()
    path.parent.mkdir(parents=True, exist_ok=True)
    start_time = perf_counter()
    if workers == 1:
        initialize_worker(coefficients)
        iterator = map(chunk_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(coefficients,),
        )
        iterator = executor.map(chunk_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"second-next finite bound failed: {json.dumps(record, sort_keys=True)}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 20 == 0:
                    handle.flush()
                    print(
                        f"second-next finite chunks: {len(records)}/{stop} "
                        f"({perf_counter()-start_time:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def aggregate_section(records: list[dict], section: str, order: int) -> dict:
    rows = [record[section][str(order)] for record in records]
    keys = (
        ("minimum_value_lower", "minimum_value_block", "min"),
        ("maximum_value_upper", "maximum_value_block", "max"),
        ("minimum_lower_margin", "minimum_lower_margin_block", "min"),
        ("minimum_upper_margin", "minimum_upper_margin_block", "min"),
    )
    result = {}
    for value_key, block_key, direction in keys:
        selected = (
            min(rows, key=lambda row: Decimal(row[value_key]))
            if direction == "min"
            else max(rows, key=lambda row: Decimal(row[value_key]))
        )
        result[value_key] = selected[value_key]
        result[block_key] = selected[block_key]
    return result


def build_artifact(records: list[dict], cache_path: Path) -> dict:
    tasks = deterministic_tasks()
    if len(records) != len(tasks):
        raise RuntimeError("second-next finite artifact requires the complete cache")
    for previous, current in zip(records, records[1:]):
        if previous["mode_right"] != current["mode_left"]:
            raise RuntimeError("second-next finite chunks have a mode gap")
    if records[0]["mode_left"] != str(MODE_START) or records[-1]["mode_right"] != str(MODE_END):
        raise RuntimeError("second-next finite cache has the wrong endpoints")

    coefficients = {
        str(order): aggregate_section(records, "coefficients", order)
        for order in COEFFICIENT_BOUNDS
    }
    normalized_jets = {
        str(order): aggregate_section(records, "normalized_jets", order)
        for order in NORMALIZED_JET_CAPS
    }
    weakest = min(
        (
            (Decimal(row[margin]), order, margin)
            for order, row in coefficients.items()
            for margin in ("minimum_lower_margin", "minimum_upper_margin")
        ),
        key=lambda item: item[0],
    )
    total_blocks = sum(record["block_count"] for record in records)
    return {
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate",
        "date": "2026-07-13",
        "status": "rigorous finite interval theorem for the second-next coefficient functions",
        "proof_boundary": (
            "This artifact proves explicit signed bounds for all seven epsilon-nine/ten "
            "scaled coefficient functions on 2<=u<=20. It does not prove their "
            "u>=20 bounds, control the exact central remainder or tails, prove the "
            "exact cumulant corridors, curvature ray, order-four entry, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "mode_interval": [str(MODE_START), str(MODE_END)],
            "block_width": str(BLOCK_WIDTH),
            "block_count": total_blocks,
            "chunk_blocks": CHUNK_BLOCKS,
            "chunk_count": len(records),
            "taylor_degree": TAYLOR_DEGREE,
            "remainder_derivative_order": TAYLOR_DEGREE + 1,
            "maximum_potential_order": MAXIMUM_POTENTIAL_ORDER,
        },
        "method": {
            "centered_model": (
                "sum_(j=0)^6 D_r^(j)(c)*(u-c)^j/j! + "
                "D_r^(7)(I)*(u-c)^7/7!"
            ),
            "point_coefficients": "Arb formal series at the exact rational block center",
            "remainder_coefficients": "Arb natural interval extension of D_r^(7)/7! on the full block",
            "potential_jets": "exact nested-series evaluation through V^(12)",
        },
        "bounds": {
            str(order): [str(floor), str(ceiling)]
            for order, (floor, ceiling) in COEFFICIENT_BOUNDS.items()
        },
        "normalized_jet_caps": {
            str(order): str(cap) for order, cap in NORMALIZED_JET_CAPS.items()
        },
        "coefficients": coefficients,
        "normalized_jets": normalized_jets,
        "weakest_margin": {
            "value": str(weakest[0]),
            "order": weakest[1],
            "side": weakest[2],
        },
        "cache": {
            "path": cache_path.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(cache_path),
            "row_count": len(records),
        },
        "source_second_next": {
            "path": SOURCE_SECOND_NEXT.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(SOURCE_SECOND_NEXT),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.py"
        ),
        "remaining_target": (
            "Prove corresponding coefficient bounds on u>=20, then use the "
            "epsilon-ten subtraction in the exact central and two-tail theorem."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Second-Next Finite Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous finite interval theorem for the second-next coefficient",
        "functions. This is not a proof of the exact cumulant corridors,",
        "order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_chunks.jsonl",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.py",
        "```",
        "",
        "## Certified Bounds",
        "",
        f"A deterministic cover of `{artifact['parameters']['block_count']}` adjacent",
        f"blocks of width `{artifact['parameters']['block_width']}` proves:",
        "",
        "```text",
    ]
    for order, bounds in artifact["bounds"].items():
        lines.append(f"{bounds[0]} < D_{order}(u) < {bounds[1]}")
    lines.extend(
        [
            "```",
            "",
            "Here `D_r` multiplies `q^-4` for orders 2-4, `q^-3` for orders",
            "5-6, and `q^-2` for orders 7-8 in the corridor-scaled",
            "epsilon-nine/ten correction.",
            "",
            "| order | certified lower | certified upper | lower margin | upper margin |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for order, row in artifact["coefficients"].items():
        lines.append(
            f"| {order} | `{row['minimum_value_lower']}` | "
            f"`{row['maximum_value_upper']}` | `{row['minimum_lower_margin']}` | "
            f"`{row['minimum_upper_margin']}` |"
        )
    lines.extend(
        [
            "",
            "The same Taylor cells also prove the normalized potential-jet caps",
            "",
            "```text",
        ]
    )
    for order, cap in artifact["normalized_jet_caps"].items():
        lines.append(f"0 < L_{order}(u) < {cap}")
    lines.extend(["```", ""])
    lines.extend(
        [
            "",
            "## Taylor Gate",
            "",
            "Each block uses a centered sixth-order Arb Taylor model with exact",
            "rational-center derivatives and a full-block seventh-derivative",
            "remainder. The normalized potential geometry is evaluated through `V^(12)`,",
            "so all 30- and 42-monomial cancellations occur before range",
            "enclosure. This is not a point sample or floating-point fit.",
            "",
            f"The weakest outward-rounded margin is `{artifact['weakest_margin']['value']}`",
            f"for order `{artifact['weakest_margin']['order']}` on the",
            f"`{artifact['weakest_margin']['side']}` side.",
            "",
            "## Remaining Boundary",
            "",
            "This closes the explicit epsilon-nine/ten coefficient layer only on",
            "`2<=u<=20`. Its asymptotic companion and the exact central and two-tail",
            "density theorem remain open.",
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/formal_core.md",
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    parser.add_argument("--max-chunks", type=int)
    parser.add_argument("--cache-only", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    tasks = deterministic_tasks()
    records = build_cache(
        args.cache,
        tasks,
        workers=args.workers,
        overwrite=args.overwrite_cache,
        max_chunks=args.max_chunks,
    )
    print(f"second-next finite cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only or len(records) != len(tasks):
        return 0
    artifact = build_artifact(records, args.cache)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four second-next finite bounds: "
        f"{artifact['parameters']['block_count']} Taylor blocks, "
        f"{len(artifact['coefficients'])} signed coefficients, "
        f"weakest margin {artifact['weakest_margin']['value']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
