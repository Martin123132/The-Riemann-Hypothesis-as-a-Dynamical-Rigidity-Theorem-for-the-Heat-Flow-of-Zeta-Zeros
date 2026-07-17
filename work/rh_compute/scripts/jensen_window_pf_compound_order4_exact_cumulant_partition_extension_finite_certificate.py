#!/usr/bin/env python3
"""Certify the finite epsilon-11 through epsilon-14 partition extension."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
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

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    CompiledCoefficients,
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    exact_lower,
    exact_upper,
    sha256,
    taylor_range,
)


DEFAULT_PARTITION_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_chunks.jsonl"
)
DEFAULT_JET_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_jet17_chunks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.md"
)
SOURCE_SECOND_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json"
)
MODE_START = Fraction(2)
MODE_END = Fraction(20)
SHIFTED_START = Fraction(39, 20)
SHIFTED_END = Fraction(401, 20)
BLOCK_WIDTH = Fraction(1, 300)
CHUNK_BLOCKS = 10
TAYLOR_DEGREE = 6
JET_TAYLOR_DEGREE = 8
MAXIMUM_POTENTIAL_ORDER = 16
PRECISION_BITS = 192
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))
PARTITION_BOUNDS = {
    11: Fraction(2),
    12: Fraction(2),
    13: Fraction(21, 10),
    14: Fraction(12, 5),
}
NEW_JET_CAPS = {
    13: Fraction(200),
    14: Fraction(400),
    15: Fraction(800),
    16: Fraction(1600),
}
SHIFTED_JET17_CAP = Fraction(4000)


@dataclass(frozen=True)
class FiniteRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


_WORKER_COEFFICIENTS: CompiledCoefficients | None = None
_WORKER_METADATA: dict[int, tuple[str, int, int | None]] | None = None
_WORKER_JET17: CompiledCoefficients | None = None


def formal_partition_extension() -> dict:
    y, z = sp.symbols("y z")
    symbols = tuple(sp.symbols("L_3:17"))
    maximum_degree = 14
    perturbation = [sp.Integer(0) for _ in range(maximum_degree + 1)]
    for order, symbol in zip(range(3, 17), symbols):
        perturbation[order - 2] = symbol * y**order / sp.factorial(order)
    exponential = [sp.Integer(0) for _ in range(maximum_degree + 1)]
    exponential[0] = sp.Integer(1)
    for degree in range(1, maximum_degree + 1):
        exponential[degree] = sp.expand(
            -sum(
                index * perturbation[index] * exponential[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )

    tilted_moments = [sp.Integer(1)]
    maximum_y_degree = max(sp.Poly(value, y).degree() for value in exponential)
    for _degree in range(maximum_y_degree):
        previous = tilted_moments[-1]
        tilted_moments.append(sp.expand(z * previous + sp.diff(previous, z)))

    def tilted_expectation(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient * tilted_moments[monomial[0]]
                for monomial, coefficient in sp.Poly(polynomial, y).terms()
            )
        )

    partitions = {
        degree: tilted_expectation(exponential[degree])
        for degree in PARTITION_BOUNDS
    }
    compiled: CompiledCoefficients = {}
    metadata: dict[int, tuple[str, int, int | None]] = {}
    partition_rows = {}
    for degree, expression in partitions.items():
        z_rows = {}
        for (z_power,), coefficient in sp.Poly(expression, z).terms():
            key = degree * 100 + z_power
            terms = []
            for powers, scalar in sp.Poly(coefficient, *symbols).terms():
                terms.append(
                    (
                        tuple(int(power) for power in powers),
                        int(sp.numer(scalar)),
                        int(sp.denom(scalar)),
                    )
                )
            compiled[key] = terms
            metadata[key] = ("partition", degree, z_power)
            z_rows[str(z_power)] = {
                "formula": sp.sstr(coefficient),
                "terms": len(terms),
            }
        partition_rows[str(degree)] = {
            "formula": sp.sstr(expression),
            "terms": len(sp.Poly(expression, *symbols, z).terms()),
            "z_degree": sp.Poly(expression, z).degree(),
            "z_coefficients": z_rows,
        }
    for order in NEW_JET_CAPS:
        key = 10_000 + order
        powers = tuple(1 if index == order - 3 else 0 for index in range(14))
        compiled[key] = [(powers, 1, 1)]
        metadata[key] = ("jet", order, None)

    finite_caps = {
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
        **NEW_JET_CAPS,
    }
    weight_norms = {}
    for degree in PARTITION_BOUNDS:
        total = sp.Rational(0)
        for powers, coefficient in sp.Poly(
            exponential[degree], *symbols, y
        ).terms():
            term = abs(sp.Rational(coefficient))
            for order, power in zip(range(3, 17), powers[:-1]):
                cap = finite_caps[order]
                term *= sp.Rational(cap.numerator, cap.denominator) ** power
            total += term
        weight_norms[str(degree)] = {
            "coefficient_norm": str(sp.factor(total)),
            "y_degree": sp.Poly(exponential[degree], y).degree(),
            "terms": len(sp.Poly(exponential[degree], *symbols, y).terms()),
        }
    return {
        "compiled": compiled,
        "metadata": metadata,
        "partition_rows": partition_rows,
        "weight_norms": weight_norms,
        "symbol_count": len(symbols),
        "partition_scalar_functions": sum(
            len(row["z_coefficients"]) for row in partition_rows.values()
        ),
        "partition_monomials": sum(
            row["terms"] for row in partition_rows.values()
        ),
        "maximum_y_degree": maximum_y_degree,
    }


def initialize_worker(
    coefficients: CompiledCoefficients,
    metadata: dict[int, tuple[str, int, int | None]],
) -> None:
    global _WORKER_COEFFICIENTS, _WORKER_METADATA, _WORKER_JET17
    _WORKER_COEFFICIENTS = coefficients
    _WORKER_METADATA = metadata
    powers = tuple(1 if index == 14 else 0 for index in range(15))
    _WORKER_JET17 = {17: [(powers, 1, 1)]}
    flint.ctx.prec = PRECISION_BITS
    flint.ctx.cap = JET_TAYLOR_DEGREE + 2


def deterministic_tasks(
    start: Fraction, end: Fraction
) -> list[tuple[int, int, int]]:
    total = (end - start) / BLOCK_WIDTH
    if total.denominator != 1:
        raise RuntimeError("finite partition range is not aligned to the grid")
    tasks = []
    first = 0
    chunk = 0
    while first < total.numerator:
        count = min(CHUNK_BLOCKS, total.numerator - first)
        tasks.append((chunk, first, count))
        first += count
        chunk += 1
    return tasks


def partition_chunk_task(task: tuple[int, int, int]) -> dict:
    if _WORKER_COEFFICIENTS is None or _WORKER_METADATA is None:
        raise RuntimeError("partition-extension worker was not initialized")
    chunk_index, first_block, block_count = task
    maxima = {degree: (flint.arb(0), None) for degree in PARTITION_BOUNDS}
    jet_maxima = {order: (flint.arb(0), None) for order in NEW_JET_CAPS}
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
        norms = {degree: flint.arb(0) for degree in PARTITION_BOUNDS}
        for key, value in ranges.items():
            kind, order, _power = _WORKER_METADATA[key]
            if kind == "partition":
                norms[order] += exact_upper(abs(value))
            else:
                if not bool(value > 0 and value < arb_rational(NEW_JET_CAPS[order])):
                    return {
                        "passed": False,
                        "failure": "normalized-jet-cap",
                        "failed_block": block_index,
                        "failed_order": order,
                        "value": value.str(40).replace("e", "E"),
                    }
                if value.upper() > jet_maxima[order][0].upper():
                    jet_maxima[order] = (exact_upper(value), block_index)
        for degree, norm in norms.items():
            if not bool(norm < arb_rational(PARTITION_BOUNDS[degree])):
                return {
                    "passed": False,
                    "failure": "partition-norm",
                    "failed_block": block_index,
                    "failed_degree": degree,
                    "value": norm.str(40).replace("e", "E"),
                }
            if norm.upper() > maxima[degree][0].upper():
                maxima[degree] = (exact_upper(norm), block_index)
    return {
        "kind": "partition",
        "chunk_index": chunk_index,
        "first_block": first_block,
        "block_count": block_count,
        "mode_left": str(MODE_START + first_block * BLOCK_WIDTH),
        "mode_right": str(
            MODE_START + (first_block + block_count) * BLOCK_WIDTH
        ),
        "passed": True,
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "block_width": str(BLOCK_WIDTH),
            "taylor_degree": TAYLOR_DEGREE,
            "maximum_potential_order": MAXIMUM_POTENTIAL_ORDER,
        },
        "partition_norms": {
            str(degree): {
                "maximum_upper": arb_upper_text(value),
                "maximum_block": block,
            }
            for degree, (value, block) in maxima.items()
        },
        "normalized_jets": {
            str(order): {
                "maximum_upper": arb_upper_text(value),
                "maximum_block": block,
            }
            for order, (value, block) in jet_maxima.items()
        },
    }


def jet_chunk_task(task: tuple[int, int, int]) -> dict:
    if _WORKER_JET17 is None:
        raise RuntimeError("partition-extension jet worker was not initialized")
    chunk_index, first_block, block_count = task
    maximum = (flint.arb(0), None)
    for offset in range(block_count):
        block_index = first_block + offset
        left = SHIFTED_START + block_index * BLOCK_WIDTH
        right = left + BLOCK_WIDTH
        value = taylor_range(
            left,
            right,
            _WORKER_JET17,
            taylor_degree=JET_TAYLOR_DEGREE,
            maximum_potential_order=17,
        )[17]
        if not bool(value > 0 and value < arb_rational(SHIFTED_JET17_CAP)):
            return {
                "passed": False,
                "failure": "shifted-jet17-cap",
                "failed_block": block_index,
                "value": value.str(40).replace("e", "E"),
            }
        if value.upper() > maximum[0].upper():
            maximum = (exact_upper(value), block_index)
    return {
        "kind": "shifted_jet17",
        "chunk_index": chunk_index,
        "first_block": first_block,
        "block_count": block_count,
        "mode_left": str(SHIFTED_START + first_block * BLOCK_WIDTH),
        "mode_right": str(
            SHIFTED_START + (first_block + block_count) * BLOCK_WIDTH
        ),
        "passed": True,
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "block_width": str(BLOCK_WIDTH),
            "taylor_degree": JET_TAYLOR_DEGREE,
            "maximum_potential_order": 17,
        },
        "normalized_jet17": {
            "maximum_upper": arb_upper_text(maximum[0]),
            "maximum_block": maximum[1],
        },
    }


def load_cache(path: Path, tasks: list[tuple[int, int, int]], kind: str) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError(f"{kind} cache has too many rows")
    for record, (chunk, first, count) in zip(records, tasks):
        if (
            record.get("kind") != kind
            or record.get("chunk_index") != chunk
            or record.get("first_block") != first
            or record.get("block_count") != count
            or record.get("passed") is not True
            or record.get("parameters", {}).get("block_width") != str(BLOCK_WIDTH)
        ):
            raise RuntimeError(f"{kind} cache mismatch at chunk {chunk}")
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, int, int]],
    kind: str,
    worker,
    coefficients: CompiledCoefficients,
    metadata: dict[int, tuple[str, int, int | None]],
    *,
    workers: int,
    overwrite: bool,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, tasks, kind)
    remaining = tasks[len(records) :]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    if workers == 1:
        initialize_worker(coefficients, metadata)
        iterator = map(worker, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(coefficients, metadata),
        )
        iterator = executor.map(worker, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"{kind} finite bound failed: {json.dumps(record, sort_keys=True)}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 20 == 0:
                    handle.flush()
                    print(
                        f"partition-extension {kind} chunks: "
                        f"{len(records)}/{len(tasks)} ({perf_counter()-started:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def aggregate_maximum(records: list[dict], section: str, key: str) -> dict:
    selected = max(
        records,
        key=lambda record: Decimal(record[section][key]["maximum_upper"]),
    )
    row = selected[section][key]
    return {
        "maximum_upper": row["maximum_upper"],
        "maximum_block": row["maximum_block"],
    }


def validate_cover(
    records: list[dict], start: Fraction, end: Fraction, tasks: list[tuple[int, int, int]]
) -> None:
    if len(records) != len(tasks):
        raise RuntimeError("finite partition artifact requires a complete cache")
    if records[0]["mode_left"] != str(start) or records[-1]["mode_right"] != str(end):
        raise RuntimeError("finite partition cache has the wrong endpoints")
    for previous, current in zip(records, records[1:]):
        if previous["mode_right"] != current["mode_left"]:
            raise RuntimeError("finite partition cache has a gap")


def build_artifact(
    algebra: dict,
    partition_records: list[dict],
    jet_records: list[dict],
    partition_cache: Path,
    jet_cache: Path,
) -> dict:
    partition_tasks = deterministic_tasks(MODE_START, MODE_END)
    jet_tasks = deterministic_tasks(SHIFTED_START, SHIFTED_END)
    validate_cover(partition_records, MODE_START, MODE_END, partition_tasks)
    validate_cover(jet_records, SHIFTED_START, SHIFTED_END, jet_tasks)
    second_finite = json.loads(SOURCE_SECOND_FINITE.read_text(encoding="utf-8"))
    if second_finite.get("parameters", {}).get("block_count") != 3600:
        raise RuntimeError("second-next finite source is not closed")

    partition_norms = {
        str(degree): aggregate_maximum(
            partition_records, "partition_norms", str(degree)
        )
        for degree in PARTITION_BOUNDS
    }
    normalized_jets = {
        str(order): aggregate_maximum(
            partition_records, "normalized_jets", str(order)
        )
        for order in NEW_JET_CAPS
    }
    shifted_jet17 = aggregate_maximum(
        jet_records, "normalized_jet17", "maximum_upper"
    ) if False else max(
        (record["normalized_jet17"] for record in jet_records),
        key=lambda row: Decimal(row["maximum_upper"]),
    )
    rows = [
        FiniteRow(
            id="co4ecpefc_01_graded_partition_extension",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="The exact graded exponential and Gaussian tilted-moment recurrences extend the partition through epsilon fourteen.",
            formula="P^[14]=P^[10]+sum_(n=11)^14 epsilon^n Z_n",
            proof_boundary="Exact formal partition algebra only.",
            diagnostics={
                "scalar_functions": algebra["partition_scalar_functions"],
                "partition_monomials": algebra["partition_monomials"],
                "maximum_y_degree": algebra["maximum_y_degree"],
            },
        ),
        FiniteRow(
            id="co4ecpefc_02_partition_taylor_cover",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="Centered Arb Taylor enclosures preserve the partition cancellations on every finite mode block.",
            formula="sum_k |[z^k]Z_n(u,z)|<B_n, n=11,...,14",
            proof_boundary="Finite interval 2<=u<=20 and |z|<=1 only.",
            diagnostics=partition_norms,
        ),
        FiniteRow(
            id="co4ecpefc_03_order_eleven_bound",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="The epsilon-eleven partition coefficient has unit-disk coefficient norm below two.",
            formula="||Z_11||_1<2",
            proof_boundary="Formal partition coefficient only.",
            diagnostics=partition_norms["11"],
        ),
        FiniteRow(
            id="co4ecpefc_04_orders_twelve_to_fourteen",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="The next three cancellation-bearing partition coefficients obey fixed unit-disk norm caps.",
            formula="||Z_12||_1<2; ||Z_13||_1<21/10; ||Z_14||_1<12/5",
            proof_boundary="Formal partition coefficients only.",
            diagnostics={key: partition_norms[key] for key in ("12", "13", "14")},
        ),
        FiniteRow(
            id="co4ecpefc_05_new_mode_jet_caps",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="The same finite Taylor cover bounds the normalized potential jets through order sixteen.",
            formula="L_13<200, L_14<400, L_15<800, L_16<1600",
            proof_boundary="Mode jets on 2<=u<=20 only.",
            diagnostics=normalized_jets,
        ),
        FiniteRow(
            id="co4ecpefc_06_shifted_jet17_cap",
            role="exact_interval_theorem",
            readiness="ready_to_apply",
            claim="A companion collar cover bounds the normalized seventeenth potential derivative.",
            formula="0<L_17(v)<4000 for 39/20<=v<=401/20",
            proof_boundary="Shifted normalized jet only; central normalization transfer is separate.",
            diagnostics=shifted_jet17,
        ),
        FiniteRow(
            id="co4ecpefc_07_formal_tail_inputs",
            role="exact_inequality_input",
            readiness="ready_to_apply",
            claim="Exact coefficient norms record the four added formal-density polynomials for tail accounting.",
            formula="deg_y E_n=3n, n=11,...,14",
            proof_boundary="Formal weight coefficient norms only.",
            diagnostics=algebra["weight_norms"],
        ),
        FiniteRow(
            id="co4ecpefc_08_central_composition_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Compose these cancellation-bearing coefficients with the epsilon-fifteen Bell remainder and exact potential Taylor remainder.",
            formula="central exact-minus-P^[10] <1/(500000*q^3)",
            proof_boundary="Open composition target; no exact cumulant corridor is asserted here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate",
        "date": "2026-07-13",
        "status": "rigorous finite partition extension and high-jet theorem",
        "proof_boundary": (
            "This artifact proves the cancellation-bearing epsilon-eleven through "
            "epsilon-fourteen partition coefficient bounds and normalized potential-jet "
            "caps needed for the finite central residual. It does not yet compose the "
            "Bell and potential remainders, prove the exact cumulant corridors, "
            "curvature ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "mode_interval": [str(MODE_START), str(MODE_END)],
            "shifted_interval": [str(SHIFTED_START), str(SHIFTED_END)],
            "block_width": str(BLOCK_WIDTH),
            "partition_blocks": sum(row["block_count"] for row in partition_records),
            "shifted_jet_blocks": sum(row["block_count"] for row in jet_records),
            "partition_taylor_degree": TAYLOR_DEGREE,
            "shifted_jet_taylor_degree": JET_TAYLOR_DEGREE,
        },
        "partition_bounds": {
            str(degree): str(bound) for degree, bound in PARTITION_BOUNDS.items()
        },
        "partition_norms": partition_norms,
        "normalized_jet_caps": {
            str(order): str(bound) for order, bound in NEW_JET_CAPS.items()
        },
        "normalized_jets": normalized_jets,
        "shifted_jet17_cap": str(SHIFTED_JET17_CAP),
        "shifted_jet17": shifted_jet17,
        "partition_rows": algebra["partition_rows"],
        "formal_weight_norms": algebra["weight_norms"],
        "cache": {
            "partition_path": partition_cache.relative_to(REPO_ROOT).as_posix(),
            "partition_sha256": sha256(partition_cache),
            "partition_records": len(partition_records),
            "jet17_path": jet_cache.relative_to(REPO_ROOT).as_posix(),
            "jet17_sha256": sha256(jet_cache),
            "jet17_records": len(jet_records),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 7,
            "open_analytic_rows": 1,
            "partition_orders": 4,
            "partition_scalar_functions": algebra["partition_scalar_functions"],
            "partition_blocks": sum(row["block_count"] for row in partition_records),
            "shifted_jet_blocks": sum(row["block_count"] for row in jet_records),
            "new_normalized_jet_caps": 5,
            "finite_partition_extension_closed": True,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            SOURCE_SECOND_FINITE.relative_to(REPO_ROOT).as_posix(): sha256(
                SOURCE_SECOND_FINITE
            )
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.py"
        ),
        "remaining_target": (
            "Use the four partition coefficient caps, the finite normalized-jet "
            "caps, and a cancellation-preserving epsilon-fifteen remainder to close "
            "the sole central unit-disk partition residual."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Partition-Extension Finite Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous finite partition extension and high-jet theorem.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.py",
        "```",
        "",
        "## Partition Extension",
        "",
        "Exact graded exponential and Gaussian tilted-moment recurrences extend the",
        "formal partition through epsilon fourteen. Centered Arb Taylor enclosures",
        "preserve the cancellations in all 78 scalar coefficient functions and prove",
        "on `2<=u<=20`, uniformly for `|z|<=1`,",
        "",
        "```text",
        "||Z_11||_1<2,",
        "||Z_12||_1<2,",
        "||Z_13||_1<21/10,",
        "||Z_14||_1<12/5.",
        "```",
        "",
        "## High Jets",
        "",
        "The finite covers also prove",
        "",
        "```text",
        "L_13<200, L_14<400, L_15<800, L_16<1600,",
        "L_17(v)<4000 for 39/20<=v<=401/20.",
        "```",
        "",
        f"The partition cover contains `{artifact['parameters']['partition_blocks']}` blocks;",
        f"the shifted collar cover contains `{artifact['parameters']['shifted_jet_blocks']}` blocks.",
        "Both caches are deterministic and source-hashed in the JSON artifact.",
        "",
        "## Remaining Boundary",
        "",
        "The remaining finite step is scalar composition: combine these four small",
        "partition coefficients with the epsilon-fifteen Bell remainder, the added",
        "formal tails, and the exact potential Taylor remainder. No exact cumulant",
        "corridor is promoted by this finite input theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.md",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-cache", type=Path, default=DEFAULT_PARTITION_CACHE)
    parser.add_argument("--jet-cache", type=Path, default=DEFAULT_JET_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    algebra = formal_partition_extension()
    coefficients = algebra.pop("compiled")
    metadata = algebra.pop("metadata")
    partition_tasks = deterministic_tasks(MODE_START, MODE_END)
    jet_tasks = deterministic_tasks(SHIFTED_START, SHIFTED_END)
    partition_records = build_cache(
        args.partition_cache,
        partition_tasks,
        "partition",
        partition_chunk_task,
        coefficients,
        metadata,
        workers=max(1, args.workers),
        overwrite=args.overwrite_cache,
    )
    jet_records = build_cache(
        args.jet_cache,
        jet_tasks,
        "shifted_jet17",
        jet_chunk_task,
        coefficients,
        metadata,
        workers=max(1, args.workers),
        overwrite=args.overwrite_cache,
    )
    artifact = build_artifact(
        algebra,
        partition_records,
        jet_records,
        args.partition_cache,
        args.jet_cache,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four finite partition extension: "
        "8 rows, 7 exact rows, 4 partition orders, 78 scalar functions, "
        "5400 partition blocks, 5430 shifted-jet blocks, 5 new jet caps"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
