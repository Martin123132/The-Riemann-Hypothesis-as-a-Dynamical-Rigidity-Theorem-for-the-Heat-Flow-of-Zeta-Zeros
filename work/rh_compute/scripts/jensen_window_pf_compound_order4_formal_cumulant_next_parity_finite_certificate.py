#!/usr/bin/env python3
"""Certify next-parity cumulant coefficient bounds on 2<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
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
import sympy as sp  # noqa: E402

from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_interval,
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    exact_lower,
    exact_upper,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_chunks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.md"
)
SOURCE_PARITY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.json"
)
MODE_START = Fraction(2)
MODE_END = Fraction(20)
BLOCK_WIDTH = Fraction(1, 100)
CHUNK_BLOCKS = 10
TAYLOR_DEGREE = 6
PRECISION_BITS = 192
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))
COEFFICIENT_BOUNDS = {
    2: (Fraction(-1, 20), Fraction(0)),
    3: (Fraction(-3, 20), Fraction(0)),
    4: (Fraction(-1, 2), Fraction(0)),
    5: (Fraction(0), Fraction(2)),
    6: (Fraction(0), Fraction(16, 5)),
    7: (Fraction(0), Fraction(5)),
    8: (Fraction(0), Fraction(61, 10)),
}


CompiledTerm = tuple[tuple[int, ...], int, int]
CompiledCoefficients = dict[int, list[CompiledTerm]]
_WORKER_COEFFICIENTS: CompiledCoefficients | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def compile_coefficient_family(
    path: Path,
    *,
    symbol_stop: int,
    expected_orders: set[int],
) -> CompiledCoefficients:
    artifact = load_json(path)
    symbols = sp.symbols(f"L_3:{symbol_stop}")
    locals_map = {str(symbol): symbol for symbol in symbols}
    compiled = {}
    for order_text, row in artifact["coefficient_rows"].items():
        expression = sp.sympify(row["scaled_coefficient"], locals=locals_map)
        terms = []
        for powers, coefficient in sp.Poly(expression, *symbols).terms():
            rational = sp.Rational(coefficient)
            terms.append((powers, int(rational.p), int(rational.q)))
        compiled[int(order_text)] = terms
    if set(compiled) != expected_orders:
        raise RuntimeError("coefficient source has the wrong orders")
    return compiled


def compile_coefficients(path: Path = SOURCE_PARITY) -> CompiledCoefficients:
    return compile_coefficient_family(
        path,
        symbol_stop=11,
        expected_orders=set(COEFFICIENT_BOUNDS),
    )


def initialize_worker(coefficients: CompiledCoefficients) -> None:
    global _WORKER_COEFFICIENTS
    _WORKER_COEFFICIENTS = coefficients
    flint.ctx.prec = PRECISION_BITS
    flint.ctx.cap = TAYLOR_DEGREE + 2


def half_power(value: object, doubled_exponent: int) -> object:
    integer_part = value ** (doubled_exponent // 2)
    return integer_part if doubled_exponent % 2 == 0 else integer_part * value.sqrt()


def nested_series_multiply(left: list[object], right: list[object], order: int) -> list[object]:
    return [
        sum(left[index] * right[degree - index] for index in range(degree + 1))
        for degree in range(order + 1)
    ]


def nested_series_exponential(values: list[object], order: int) -> list[object]:
    result = [flint.arb(0) for _ in range(order + 1)]
    result[0] = values[0].exp()
    for degree in range(1, order + 1):
        result[degree] = (
            sum(
                flint.arb(index) * values[index] * result[degree - index]
                for index in range(1, degree + 1)
            )
            * (flint.arb(1) / degree)
        )
    return result


def nested_series_logarithm(values: list[object], order: int) -> list[object]:
    """Logarithm for an auxiliary series with arb_series coefficients."""
    result = [flint.arb(0) for _ in range(order + 1)]
    result[0] = values[0].log()
    inverse_constant = values[0].inv()
    for degree in range(1, order + 1):
        prior = sum(
            flint.arb(index) * result[index] * values[degree - index]
            for index in range(1, degree)
        )
        result[degree] = (
            (flint.arb(degree) * values[degree] - prior)
            * inverse_constant
            * (flint.arb(1) / degree)
        )
    return result


def potential_jet_series(mode: flint.arb_series, order: int = 10) -> list[flint.arb_series]:
    """Evaluate x-jets of the exact potential with an outer Arb series mode."""
    mode_series = [
        mode * (flint.arb(1) / ((2**degree) * math.factorial(degree)))
        for degree in range(order + 1)
    ]
    q_series = [
        flint.arb.pi() * value
        for value in nested_series_exponential(
            [4 * value for value in mode_series], order
        )
    ]
    denominator = [2 * value for value in q_series]
    denominator[0] -= 3
    log_denominator = nested_series_logarithm(denominator, order)
    log_mode = nested_series_logarithm(mode_series, order)
    potential = [
        100 * value
        for value in nested_series_multiply(mode_series, mode_series, order)
    ]
    for degree in range(order + 1):
        potential[degree] += q_series[degree] - 5 * mode_series[degree]
        potential[degree] -= log_denominator[degree]
        potential[degree] -= log_mode[degree]
    return [
        potential[degree] * math.factorial(degree)
        for degree in range(order + 1)
    ]


def coefficient_series(
    mode: flint.arb_series,
    coefficients: CompiledCoefficients,
    *,
    maximum_potential_order: int = 10,
) -> dict[int, flint.arb_series]:
    jet = potential_jet_series(mode, maximum_potential_order)
    q = flint.arb.pi() * (4 * mode).exp()
    curvature = jet[2]
    normalized = [
        jet[order]
        * half_power(q, order - 2)
        * half_power(curvature, order).inv()
        for order in range(3, maximum_potential_order + 1)
    ]
    values = {}
    for order, terms in coefficients.items():
        result = flint.arb_series([], prec=mode.prec)
        for powers, numerator, denominator in terms:
            term = flint.arb(numerator) / denominator
            for value, power in zip(normalized, powers):
                if power:
                    term *= value**power
            result += term
        values[order] = result
    return values


def integer_power(value: flint.arb, exponent: int) -> flint.arb:
    """Avoid the python-flint 0.8 indeterminate power path for zero-containing balls."""
    result = flint.arb(1)
    for _ in range(exponent):
        result *= value
    return result


def taylor_range(
    left: Fraction,
    right: Fraction,
    coefficients: CompiledCoefficients,
    *,
    taylor_degree: int = TAYLOR_DEGREE,
    maximum_potential_order: int = 10,
) -> dict[int, flint.arb]:
    precision = taylor_degree + 2
    center = (left + right) / 2
    radius = (right - left) / 2
    point_mode = flint.arb_series(
        [arb_rational(center), flint.arb(1)], prec=precision
    )
    interval_mode = flint.arb_series(
        [arb_interval(left, right), flint.arb(1)], prec=precision
    )
    point_series = coefficient_series(
        point_mode,
        coefficients,
        maximum_potential_order=maximum_potential_order,
    )
    interval_series = coefficient_series(
        interval_mode,
        coefficients,
        maximum_potential_order=maximum_potential_order,
    )
    displacement = arb_interval(-radius, radius)
    remainder_power = integer_power(displacement, taylor_degree + 1)
    ranges = {}
    for order in coefficients:
        point_coefficients = point_series[order].coeffs()
        value = flint.arb(0)
        for degree in range(taylor_degree, -1, -1):
            value = value * displacement + point_coefficients[degree]
        value += interval_series[order].coeffs()[taylor_degree + 1] * remainder_power
        ranges[order] = value
    return ranges


def empty_extrema() -> dict:
    return {
        "minimum_value_lower": None,
        "minimum_value_block": None,
        "maximum_value_upper": None,
        "maximum_value_block": None,
        "minimum_lower_margin": None,
        "minimum_lower_margin_block": None,
        "minimum_upper_margin": None,
        "minimum_upper_margin_block": None,
    }


def update_extrema(
    extrema: dict,
    value: flint.arb,
    lower_margin: flint.arb,
    upper_margin: flint.arb,
    block_index: int,
) -> None:
    candidates = (
        ("minimum_value_lower", "minimum_value_block", exact_lower(value), "min"),
        ("maximum_value_upper", "maximum_value_block", exact_upper(value), "max"),
        (
            "minimum_lower_margin",
            "minimum_lower_margin_block",
            exact_lower(lower_margin),
            "min",
        ),
        (
            "minimum_upper_margin",
            "minimum_upper_margin_block",
            exact_lower(upper_margin),
            "min",
        ),
    )
    for value_key, block_key, candidate, direction in candidates:
        current = extrema[value_key]
        if current is None or (
            direction == "min" and candidate.lower() < current.lower()
        ) or (direction == "max" and candidate.upper() > current.upper()):
            extrema[value_key] = candidate
            extrema[block_key] = block_index


def serialize_extrema(extrema: dict) -> dict:
    return {
        key: (
            arb_lower_text(value)
            if key.startswith("minimum") and not key.endswith("block")
            else arb_upper_text(value)
            if key == "maximum_value_upper"
            else value
        )
        for key, value in extrema.items()
    }


def chunk_task(task: tuple[int, int, int]) -> dict:
    if _WORKER_COEFFICIENTS is None:
        raise RuntimeError("next-parity worker was not initialized")
    chunk_index, first_block, block_count = task
    extrema = {order: empty_extrema() for order in COEFFICIENT_BOUNDS}
    for offset in range(block_count):
        block_index = first_block + offset
        left = MODE_START + block_index * BLOCK_WIDTH
        right = left + BLOCK_WIDTH
        ranges = taylor_range(left, right, _WORKER_COEFFICIENTS)
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
        },
        "coefficients": {
            str(order): serialize_extrema(values)
            for order, values in extrema.items()
        },
    }


def deterministic_tasks() -> list[tuple[int, int, int]]:
    total_fraction = (MODE_END - MODE_START) / BLOCK_WIDTH
    if total_fraction.denominator != 1:
        raise RuntimeError("next-parity finite range is not aligned to the grid")
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
        raise RuntimeError("next-parity finite cache has too many rows")
    for record, (chunk, first, count) in zip(records, tasks):
        if (
            record.get("chunk_index") != chunk
            or record.get("first_block") != first
            or record.get("block_count") != count
            or record.get("passed") is not True
            or record.get("parameters", {}).get("taylor_degree") != TAYLOR_DEGREE
        ):
            raise RuntimeError(f"next-parity finite cache mismatch at chunk {chunk}")
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
                        f"next-parity finite bound failed: {json.dumps(record, sort_keys=True)}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 20 == 0:
                    handle.flush()
                    print(
                        f"next-parity finite chunks: {len(records)}/{stop} "
                        f"({perf_counter()-start_time:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def aggregate_coefficient(records: list[dict], order: int) -> dict:
    rows = [record["coefficients"][str(order)] for record in records]
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
        raise RuntimeError("next-parity finite artifact requires the complete cache")
    for previous, current in zip(records, records[1:]):
        if previous["mode_right"] != current["mode_left"]:
            raise RuntimeError("next-parity finite chunks have a mode gap")
    if records[0]["mode_left"] != str(MODE_START) or records[-1]["mode_right"] != str(MODE_END):
        raise RuntimeError("next-parity finite cache has the wrong endpoints")

    coefficients = {
        str(order): aggregate_coefficient(records, order)
        for order in COEFFICIENT_BOUNDS
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
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate",
        "date": "2026-07-13",
        "status": "rigorous finite interval theorem for the next-parity coefficient functions",
        "proof_boundary": (
            "This artifact proves explicit signed bounds for all seven first omitted "
            "scaled formal coefficient functions on 2<=u<=20. It does not prove "
            "their u>=20 bounds, control terms beyond epsilon order eight, bound the "
            "exact density or tails, prove the exact curvature ray, order-four entry, "
            "PF-infinity, RH, or Lambda<=0."
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
        },
        "method": {
            "centered_model": (
                "sum_(j=0)^6 C_r^(j)(c)*(u-c)^j/j! + "
                "C_r^(7)(I)*(u-c)^7/7!"
            ),
            "point_coefficients": "Arb formal series at the exact rational block center",
            "remainder_coefficients": "Arb natural interval extension of C_r^(7)/7! on the full block",
            "potential_jets": "exact symbolic-series evaluation through V^(10)",
        },
        "bounds": {
            str(order): [str(floor), str(ceiling)]
            for order, (floor, ceiling) in COEFFICIENT_BOUNDS.items()
        },
        "coefficients": coefficients,
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
        "source_next_parity": {
            "path": SOURCE_PARITY.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(SOURCE_PARITY),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.py"
        ),
        "remaining_target": (
            "Prove corresponding coefficient bounds on u>=20, then bound the "
            "beyond-epsilon-eight central remainder and both exact-density tails."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Next-Parity Finite Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous finite interval theorem for the next-parity coefficient",
        "functions. This is not a proof of the exact cumulant ray, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.json",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_chunks.jsonl",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.py",
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
        lines.append(f"{bounds[0]} < C_{order}(u) < {bounds[1]}")
    lines.extend(
        [
            "```",
            "",
            "Here `C_r` is the exact rational polynomial multiplying the first",
            "omitted corridor-scaled power: `q^-3` for orders 2-4, `q^-2` for",
            "orders 5-6, and `q^-1` for orders 7-8.",
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
            "## Taylor Gate",
            "",
            "Direct interval substitution loses the cancellations in the 15-22-term",
            "coefficient polynomials. Each block instead uses a centered sixth-order",
            "Arb Taylor model. Exact rational-center coefficients enclose derivatives",
            "through order six, while a full-block natural interval extension encloses",
            "the seventh derivative remainder. Thus the range enclosure is rigorous,",
            "not a point sample or floating-point fit.",
            "",
            f"The weakest outward-rounded margin is `{artifact['weakest_margin']['value']}`",
            f"for order `{artifact['weakest_margin']['order']}` on the",
            f"`{artifact['weakest_margin']['side']}` side.",
            "",
            "## Remaining Boundary",
            "",
            "This closes only the explicit next-parity coefficient layer on the finite",
            "interval. The asymptotic coefficient ray, all terms beyond epsilon eight,",
            "and the exact-density central and two-tail estimates remain open.",
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.md",
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
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
    print(f"next-parity finite cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only or len(records) != len(tasks):
        return 0
    artifact = build_artifact(records, args.cache)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four next-parity finite bounds: "
        f"{artifact['parameters']['block_count']} Taylor blocks, "
        f"{len(artifact['coefficients'])} signed coefficients, "
        f"weakest margin {artifact['weakest_margin']['value']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
