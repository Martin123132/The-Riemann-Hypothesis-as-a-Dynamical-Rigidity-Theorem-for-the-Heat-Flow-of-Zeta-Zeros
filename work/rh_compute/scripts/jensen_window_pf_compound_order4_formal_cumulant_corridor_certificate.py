#!/usr/bin/env python3
"""Certify the epsilon^6 formal cumulant corridors on 2<=u<=20."""

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

from jensen_window_pf_compound_order4_gaussian_cumulant_ray_target import (  # noqa: E402
    CORRIDOR_CAPS,
    CORRIDOR_K2,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
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
    "jensen_window_pf_compound_order4_formal_cumulant_corridor_chunks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md"
)
SOURCE_TARGET = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json"
)
MODE_START = Fraction(2)
MODE_END = Fraction(20)
BLOCK_WIDTH = Fraction(1, 100_000)
CHUNK_BLOCKS = 1_000
PRECISION_BITS = 192
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


CompiledTerm = tuple[tuple[int, ...], int, int]
CompiledSeries = dict[int, list[tuple[int, list[CompiledTerm]]]]
_WORKER_SERIES: CompiledSeries | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def compile_formal_series(path: Path = SOURCE_TARGET) -> CompiledSeries:
    artifact = load_json(path)
    symbols = sp.symbols("L_3:9")
    locals_map = {str(symbol): symbol for symbol in symbols}
    compiled: CompiledSeries = {}
    for order_text, row in artifact["exact"]["cumulants"].items():
        series = []
        for term in row["terms"]:
            expression = sp.sympify(term["coefficient"], locals=locals_map)
            polynomial = sp.Poly(expression, *symbols)
            monomials: list[CompiledTerm] = []
            for powers, coefficient in polynomial.terms():
                rational = sp.Rational(coefficient)
                monomials.append(
                    (powers, int(rational.p), int(rational.q))
                )
            series.append((int(term["epsilon_power"]), monomials))
        compiled[int(order_text)] = series
    return compiled


def initialize_worker(series: CompiledSeries) -> None:
    global _WORKER_SERIES
    _WORKER_SERIES = series
    flint.ctx.prec = PRECISION_BITS


def evaluate_polynomial(
    terms: list[CompiledTerm], values: list[flint.arb]
) -> flint.arb:
    result = flint.arb(0)
    for powers, numerator, denominator in terms:
        term = flint.arb(numerator) / denominator
        for value, power in zip(values, powers):
            if power:
                term *= value**power
        result += term
    return result


def formal_scaled_quantities(
    left: Fraction, right: Fraction, series: CompiledSeries
) -> dict[str, flint.arb]:
    mode = arb_interval(left, right)
    jet = potential_jet_arb(mode, 8)
    q = flint.arb.pi() * (4 * mode).exp()
    epsilon = 1 / q.sqrt()
    curvature = jet[2]
    normalized = [
        jet[order]
        / curvature ** (flint.arb(order) / 2)
        / epsilon ** (order - 2)
        for order in range(3, 9)
    ]

    quantities: dict[str, flint.arb] = {}
    kappa2_scaled = flint.arb(0)
    for epsilon_power, terms in series[2]:
        kappa2_scaled += epsilon ** (epsilon_power - 2) * evaluate_polynomial(
            terms, normalized
        )
    quantities["q_kappa2_minus_1"] = kappa2_scaled

    for order in range(3, 9):
        scaled = flint.arb(0)
        for epsilon_power, terms in series[order]:
            scaled += epsilon ** (epsilon_power - (order - 2)) * evaluate_polynomial(
                terms, normalized
            )
        scaled *= arb_rational(
            Fraction((-1) ** order, math.factorial(order - 2))
        )
        quantities[f"scaled_kappa_{order}"] = scaled
    return quantities


def corridor_bounds() -> dict[str, tuple[Fraction, Fraction]]:
    bounds = {"q_kappa2_minus_1": CORRIDOR_K2}
    bounds.update(
        {
            f"scaled_kappa_{order}": (Fraction(1), cap)
            for order, cap in CORRIDOR_CAPS.items()
        }
    )
    return bounds


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
    if _WORKER_SERIES is None:
        raise RuntimeError("formal-series worker was not initialized")
    chunk_index, first_block, block_count = task
    bounds = corridor_bounds()
    extrema = {name: empty_extrema() for name in bounds}
    for offset in range(block_count):
        block_index = first_block + offset
        left = MODE_START + block_index * BLOCK_WIDTH
        right = left + BLOCK_WIDTH
        quantities = formal_scaled_quantities(left, right, _WORKER_SERIES)
        for name, value in quantities.items():
            floor, ceiling = bounds[name]
            lower_margin = value - arb_rational(floor)
            upper_margin = arb_rational(ceiling) - value
            if not bool(lower_margin > 0 and upper_margin > 0):
                return {
                    "chunk_index": chunk_index,
                    "first_block": first_block,
                    "block_count": block_count,
                    "passed": False,
                    "failed_block": block_index,
                    "failed_quantity": name,
                    "value": value.str(40).replace("e", "E"),
                    "lower_margin": lower_margin.str(40).replace("e", "E"),
                    "upper_margin": upper_margin.str(40).replace("e", "E"),
                }
            update_extrema(
                extrema[name], value, lower_margin, upper_margin, block_index
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
        },
        "quantities": {
            name: serialize_extrema(values) for name, values in extrema.items()
        },
    }


def deterministic_tasks() -> list[tuple[int, int, int]]:
    total_fraction = (MODE_END - MODE_START) / BLOCK_WIDTH
    if total_fraction.denominator != 1:
        raise RuntimeError("formal corridor range is not aligned to the block grid")
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
        raise RuntimeError("formal corridor cache has too many rows")
    for record, (chunk, first, count) in zip(records, tasks):
        if (
            record.get("chunk_index") != chunk
            or record.get("first_block") != first
            or record.get("block_count") != count
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"formal corridor cache mismatch at chunk {chunk}")
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

    compiled = compile_formal_series()
    path.parent.mkdir(parents=True, exist_ok=True)
    start_time = perf_counter()
    if workers == 1:
        initialize_worker(compiled)
        iterator = map(chunk_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(compiled,),
        )
        iterator = executor.map(chunk_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"formal corridor failed: {json.dumps(record, sort_keys=True)}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 100 == 0:
                    handle.flush()
                    print(
                        f"formal-corridor chunks: {len(records)}/{stop} "
                        f"({perf_counter()-start_time:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def cache_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def aggregate_quantity(records: list[dict], name: str) -> dict:
    rows = [record["quantities"][name] for record in records]
    keys = (
        ("minimum_value_lower", "minimum_value_block", "min"),
        ("maximum_value_upper", "maximum_value_block", "max"),
        ("minimum_lower_margin", "minimum_lower_margin_block", "min"),
        ("minimum_upper_margin", "minimum_upper_margin_block", "min"),
    )
    result = {}
    for value_key, block_key, direction in keys:
        selected = min(rows, key=lambda row: Decimal(row[value_key])) if direction == "min" else max(
            rows, key=lambda row: Decimal(row[value_key])
        )
        result[value_key] = selected[value_key]
        result[block_key] = selected[block_key]
    return result


def build_artifact(records: list[dict], cache_path: Path) -> dict:
    tasks = deterministic_tasks()
    if len(records) != len(tasks):
        raise RuntimeError("formal corridor artifact requires the complete chunk cache")
    for previous, current in zip(records, records[1:]):
        if previous["mode_right"] != current["mode_left"]:
            raise RuntimeError("formal corridor chunks have a mode gap")
    if records[0]["mode_left"] != str(MODE_START) or records[-1]["mode_right"] != str(MODE_END):
        raise RuntimeError("formal corridor cache has the wrong endpoints")

    quantities = {
        name: aggregate_quantity(records, name) for name in corridor_bounds()
    }
    weakest = min(
        (
            (Decimal(row[margin]), name, margin)
            for name, row in quantities.items()
            for margin in ("minimum_lower_margin", "minimum_upper_margin")
        ),
        key=lambda item: item[0],
    )
    total_blocks = sum(record["block_count"] for record in records)
    return {
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate",
        "date": "2026-07-13",
        "status": "rigorous interval certificate for the epsilon^6 formal cumulant model",
        "proof_boundary": (
            "This artifact proves that the exact epsilon^6 formal cumulant polynomial "
            "lies strictly inside the candidate corridors on 2<=u<=20. It does not "
            "bound the difference between the formal and exact cumulants, prove the "
            "u>=20 formal ray, prove the exact u>=2 curvature ray, order-four entry, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "mode_interval": [str(MODE_START), str(MODE_END)],
            "block_width": str(BLOCK_WIDTH),
            "block_count": total_blocks,
            "chunk_blocks": CHUNK_BLOCKS,
            "chunk_count": len(records),
        },
        "corridors": {
            name: [str(floor), str(ceiling)]
            for name, (floor, ceiling) in corridor_bounds().items()
        },
        "quantities": quantities,
        "weakest_margin": {
            "value": str(weakest[0]),
            "quantity": weakest[1],
            "side": weakest[2],
        },
        "cache": {
            "path": cache_path.relative_to(REPO_ROOT).as_posix(),
            "sha256": cache_sha256(cache_path),
            "row_count": len(records),
        },
        "source_formal_target": {
            "path": SOURCE_TARGET.relative_to(REPO_ROOT).as_posix(),
            "sha256": cache_sha256(SOURCE_TARGET),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.py"
        ),
        "remaining_target": (
            "Prove the formal corridors on u>=20 and bound the exact-minus-formal "
            "central and tail remainders inside the certified margins."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Formal Cumulant Corridor Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous interval certificate for the epsilon-six formal cumulant",
        "model. This is not a proof of the exact cumulant corridors, the complete",
        "curvature ray, order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.json",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_corridor_chunks.jsonl",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.py",
        "```",
        "",
        "## Certified Formal Range",
        "",
        f"`{artifact['parameters']['block_count']}` adjacent Arb blocks of width",
        f"`{artifact['parameters']['block_width']}` cover `2<=u<=20`. On every block,",
        "the exact epsilon-six formal cumulant polynomial lies strictly inside all",
        "seven candidate corridors.",
        "",
        "| quantity | value lower | value upper | lower margin | upper margin |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, row in artifact["quantities"].items():
        lines.append(
            f"| `{name}` | `{row['minimum_value_lower']}` | "
            f"`{row['maximum_value_upper']}` | `{row['minimum_lower_margin']}` | "
            f"`{row['minimum_upper_margin']}` |"
        )
    lines.extend(
        [
            "",
            f"The weakest outward-rounded margin is `{artifact['weakest_margin']['value']}`",
            f"for `{artifact['weakest_margin']['quantity']}` on the",
            f"`{artifact['weakest_margin']['side']}` side.",
            "",
            "## Proof Boundary",
            "",
            "This certificate concerns the exact truncated formal polynomial, not the",
            "exact standardized density. Two upgrades remain: prove the formal corridor",
            "for `u>=20`, then prove that the central expansion remainder and both tails",
            "fit inside the displayed finite and asymptotic margins.",
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
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
    print(f"formal-corridor cache rows: {len(records)}/{len(tasks)}")
    if args.cache_only or len(records) != len(tasks):
        return 0
    artifact = build_artifact(records, args.cache)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "certified order-four formal cumulant corridors: "
        f"{artifact['parameters']['block_count']} blocks, "
        f"weakest margin {artifact['weakest_margin']['value']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
