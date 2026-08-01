#!/usr/bin/env python3
"""Certify the complete Q208 bottom edge by fixed-time phase cells."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import time
from typing import Sequence

import psutil

import jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate as bridge


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q208_bottom_phase_cell_certificate"
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

SCHEMA = "newman_q208_bottom_phase_cell_v1"
DATE = "2026-07-25"
TIME_VALUE = Fraction(1, 1040)
X_LOWER = Fraction(0)
X_UPPER = Fraction(246)
INITIAL_X_STEP = Fraction(1, 2)
TAYLOR_X_ORDER = 24
MAX_OSCILLATORY_MOMENT = TAYLOR_X_ORDER + 1
MAX_ABSOLUTE_MOMENT = TAYLOR_X_ORDER + 2
MAX_SUBDIVISION_DEPTH = 12
ENDPOINT_DIGITS = 120
BASELINE_SECONDS = 6
BASELINE_LIMIT = 60.0
CPU_PARK_THRESHOLD = 75.0
CPU_CONSECUTIVE_PARK_SAMPLES = 2
DEFAULT_RUNTIME_LIMIT_SECONDS = 10_800.0

compact = bridge.compact


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_range(
    lower: Fraction,
    upper: Fraction,
    step: Fraction,
) -> list[tuple[Fraction, Fraction]]:
    rows: list[tuple[Fraction, Fraction]] = []
    current = lower
    while current < upper:
        following = min(current + step, upper)
        rows.append((current, following))
        current = following
    return rows


def panel_tasks() -> list[dict]:
    return [
        {
            "sequence": sequence,
            "x_low": str(x_low),
            "x_high": str(x_high),
        }
        for sequence, (x_low, x_high) in enumerate(
            fraction_range(X_LOWER, X_UPPER, INITIAL_X_STEP),
            start=1,
        )
    ]


def source_hashes() -> dict:
    return {
        "bridge_builder": file_hash(Path(bridge.__file__).resolve()),
        "six_term_tail": file_hash(bridge.TAIL_SOURCE),
        "outer_tail": file_hash(bridge.OUTER_SOURCE),
    }


def contract_payload() -> dict:
    return {
        "schema": SCHEMA,
        "time": str(TIME_VALUE),
        "x_domain": [str(X_LOWER), str(X_UPPER)],
        "initial_x_step": str(INITIAL_X_STEP),
        "taylor_x_order": TAYLOR_X_ORDER,
        "max_subdivision_depth": MAX_SUBDIVISION_DEPTH,
        "precision_bits": bridge.PRECISION_BITS,
        "retained_terms": bridge.RETAINED_TERMS,
        "first_omitted": bridge.FIRST_OMITTED,
        "integration_cutoff": str(bridge.INTEGRATION_CUTOFF),
        "proxy": "F_t(x)=16*(1+x^4)*H_t(x)",
        "proxy_derivative": (
            "F_t'(x)=64*x^3*H_t(x)+16*(1+x^4)*H_t'(x)"
        ),
        "source_sha256": source_hashes(),
    }


def contract_hash() -> str:
    return sha256(canonical_json(contract_payload()).encode()).hexdigest()


def serialize_interval(value) -> str:
    if not value.is_finite():
        raise RuntimeError(f"nonfinite interval: {value}")
    return value.str(ENDPOINT_DIGITS, more=True)


def fraction_to_arb(value: Fraction):
    return compact.arb(value.numerator) / value.denominator


def dyadic_fraction(value) -> Fraction:
    mantissa, exponent = value.mid().man_exp()
    mantissa = int(mantissa)
    exponent = int(exponent)
    if exponent >= 0:
        return Fraction(mantissa * (1 << exponent))
    return Fraction(mantissa, 1 << (-exponent))


def interval_midpoint_fraction(value) -> Fraction:
    point = dyadic_fraction(value)
    if not value.contains(fraction_to_arb(point)):
        raise RuntimeError("dyadic midpoint escaped its interval")
    return point


def interval_intersection(left, right):
    try:
        return left.intersection(right)
    except ValueError as exc:
        raise RuntimeError(
            f"adjacent phase cells do not intersect: {left} and {right}"
        ) from exc


def interval_contains_zero(value) -> bool:
    return value.contains(compact.arb(0))


def exact_cross(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def exact_dot(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def open_ray_crossing(vertices: Sequence[tuple[Fraction, Fraction]]) -> int:
    total = 0
    for left, right in zip(vertices, vertices[1:]):
        determinant = exact_cross(left, right)
        if (
            left == (0, 0)
            or right == (0, 0)
            or (determinant == 0 and exact_dot(left, right) <= 0)
        ):
            raise RuntimeError("witness edge meets the origin")
        if left[1] <= 0 < right[1] and determinant > 0:
            total += 1
        elif right[1] <= 0 < left[1] and determinant < 0:
            total -= 1
    return total


class FixedTimePanelCertifier:
    def __init__(self, time_value: Fraction = TIME_VALUE) -> None:
        compact.flint.ctx.prec = bridge.PRECISION_BITS
        self.base = bridge.SharedPanelCertifier()
        self.time_value = time_value
        self.time_ball = compact.acb(
            compact.fraction_decimal(self.time_value)
        )
        self.models_built = 0
        self.oscillatory_integrals = 0

    def oscillatory_integral(
        self,
        moment: int,
        frequency: Fraction,
    ):
        frequency_ball = compact.acb(
            compact.fraction_decimal(frequency)
        )

        def integrand(u, _analytic: bool):
            return (
                u**moment
                * (self.time_ball * u * u).exp()
                * self.base.retained_kernel(u)
                * (compact.acb(0, 1) * frequency_ball * u).exp()
            )

        value = compact.acb.integral(
            integrand,
            self.base.zero,
            self.base.cutoff,
            **self.base.options,
        )
        self.oscillatory_integrals += 1
        return value + self.base.cutoff_tail_complex(moment)

    def remainder(self, derivative_order: int, radius: Fraction):
        radius_ball = compact.arb(compact.fraction_decimal(radius))
        moment = derivative_order + TAYLOR_X_ORDER + 1
        return (
            compact.arb_power(radius_ball, TAYLOR_X_ORDER + 1)
            / math.factorial(TAYLOR_X_ORDER + 1)
            * self.base.absolute_moment(moment, self.time_value).upper()
        ).upper()

    def build_model(
        self,
        x_low: Fraction,
        x_high: Fraction,
    ) -> dict:
        center = (x_low + x_high) / 2
        radius = (x_high - x_low) / 2
        integrals = [
            self.oscillatory_integral(moment, center)
            for moment in range(MAX_OSCILLATORY_MOMENT + 1)
        ]
        coefficients: dict[int, list] = {}
        for derivative_order in (0, 1):
            coefficients[derivative_order] = [
                (
                    compact.acb(0, 1) ** (derivative_order + order)
                    * integrals[derivative_order + order]
                ).real
                / math.factorial(order)
                for order in range(TAYLOR_X_ORDER + 1)
            ]
        remainders = {
            derivative_order: self.remainder(derivative_order, radius)
            for derivative_order in (0, 1)
        }
        self.models_built += 1
        return {
            "x_low": x_low,
            "x_high": x_high,
            "x_center": center,
            "x_radius": radius,
            "coefficients": coefficients,
            "remainders": remainders,
            "minimum_integral_accuracy_bits": min(
                integral.rel_accuracy_bits() for integral in integrals
            ),
        }

    def transform_box(
        self,
        model: dict,
        derivative_order: int,
        x_low: Fraction,
        x_high: Fraction,
    ):
        xi = bridge.interval_ball(
            x_low - model["x_center"],
            x_high - model["x_center"],
        )
        value = compact.arb(0)
        for coefficient in reversed(model["coefficients"][derivative_order]):
            value = value * xi + coefficient
        return value + compact.arb(
            0,
            model["remainders"][derivative_order].str(
                ENDPOINT_DIGITS,
                more=True,
            ),
        )

    def evaluate_box(
        self,
        model: dict,
        x_low: Fraction,
        x_high: Fraction,
        depth: int,
    ) -> dict:
        h_retained = self.transform_box(model, 0, x_low, x_high)
        h_prime_retained = self.transform_box(
            model,
            1,
            x_low,
            x_high,
        )
        x_box = bridge.interval_ball(x_low, x_high)
        f_retained = 16 * (1 + x_box**4) * h_retained
        f_prime_retained = (
            64 * x_box**3 * h_retained
            + 16 * (1 + x_box**4) * h_prime_retained
        )

        x_ceiling = compact.arb(compact.fraction_decimal(x_high))
        tail_f = (
            16 * (1 + x_ceiling**4) * self.base.e0
        ).upper()
        tail_f_prime = (
            64 * x_ceiling**3 * self.base.e0
            + 16 * (1 + x_ceiling**4) * self.base.e1
        ).upper()
        full_f = f_retained + compact.arb(
            0,
            tail_f.str(ENDPOINT_DIGITS, more=True),
        )
        full_f_prime = f_prime_retained + compact.arb(
            0,
            tail_f_prime.str(ENDPOINT_DIGITS, more=True),
        )
        f_margin = full_f.abs_lower()
        f_prime_margin = full_f_prime.abs_lower()
        certified = bool(f_margin > 0 or f_prime_margin > 0)
        if f_margin > 0 and (
            f_prime_margin <= 0 or f_margin > f_prime_margin
        ):
            branch = "value"
        elif f_prime_margin > 0:
            branch = "derivative"
        else:
            branch = None
        norm_lower = (f_margin**2 + f_prime_margin**2).sqrt()
        return {
            "depth": depth,
            "x_low": str(x_low),
            "x_high": str(x_high),
            "h_retained": serialize_interval(h_retained),
            "h_prime_retained": serialize_interval(h_prime_retained),
            "f_retained": serialize_interval(f_retained),
            "f_prime_retained": serialize_interval(f_prime_retained),
            "tail_f_upper": serialize_interval(tail_f),
            "tail_f_prime_upper": serialize_interval(tail_f_prime),
            "full_f": serialize_interval(full_f),
            "full_f_prime": serialize_interval(full_f_prime),
            "closest_point_norm_lower": serialize_interval(norm_lower),
            "certified": certified,
            "branch": branch,
        }

    def certify_panel(self, task: dict) -> dict:
        x_low = Fraction(task["x_low"])
        x_high = Fraction(task["x_high"])
        started = time.monotonic()
        model = self.build_model(x_low, x_high)
        queue: list[tuple[Fraction, Fraction, int]] = [
            (x_low, x_high, 0)
        ]
        cursor = 0
        leaves: list[dict] = []
        unresolved: list[dict] = []
        subdivisions = 0
        while cursor < len(queue):
            left, right, depth = queue[cursor]
            cursor += 1
            evaluation = self.evaluate_box(model, left, right, depth)
            if evaluation["certified"]:
                leaves.append(evaluation)
            elif depth < MAX_SUBDIVISION_DEPTH:
                middle = (left + right) / 2
                queue.extend(
                    [
                        (left, middle, depth + 1),
                        (middle, right, depth + 1),
                    ]
                )
                subdivisions += 1
            else:
                unresolved.append(evaluation)
        leaves.sort(key=lambda row: Fraction(row["x_low"]))
        unresolved.sort(key=lambda row: Fraction(row["x_low"]))
        return {
            "status": "certified" if not unresolved else "unresolved",
            "evaluated_boxes": cursor,
            "certified_leaf_boxes": len(leaves),
            "subdivisions": subdivisions,
            "unresolved_boxes": len(unresolved),
            "maximum_depth": max(
                [
                    row["depth"]
                    for row in leaves + unresolved
                ],
                default=0,
            ),
            "branch_counts": dict(
                Counter(row["branch"] for row in leaves)
            ),
            "model": {
                "x_center": str(model["x_center"]),
                "x_radius": str(model["x_radius"]),
                "oscillatory_moments": MAX_OSCILLATORY_MOMENT + 1,
                "minimum_integral_accuracy_bits": (
                    model["minimum_integral_accuracy_bits"]
                ),
                "h_remainder": serialize_interval(
                    model["remainders"][0]
                ),
                "h_prime_remainder": serialize_interval(
                    model["remainders"][1]
                ),
            },
            "certified_leaves": leaves,
            "unresolved_leaves": unresolved,
            "elapsed_seconds": time.monotonic() - started,
        }


def cache_record(
    task: dict,
    result: dict,
    previous_hash: str,
) -> dict:
    payload = {
        "schema": SCHEMA,
        "contract_sha256": contract_hash(),
        "sequence": task["sequence"],
        "task": task,
        "previous_hash": previous_hash,
        "result": result,
    }
    payload["record_sha256"] = sha256(
        canonical_json(payload).encode()
    ).hexdigest()
    return payload


def validate_cache_record(
    record: dict,
    task: dict,
    previous_hash: str,
) -> None:
    if record.get("schema") != SCHEMA:
        raise RuntimeError("cache schema mismatch")
    if record.get("contract_sha256") != contract_hash():
        raise RuntimeError("cache contract hash mismatch")
    if record.get("sequence") != task["sequence"]:
        raise RuntimeError("cache sequence mismatch")
    if record.get("task") != task:
        raise RuntimeError("cache task mismatch")
    if record.get("previous_hash") != previous_hash:
        raise RuntimeError("cache hash chain is broken")
    claimed = record.get("record_sha256")
    replay = dict(record)
    replay.pop("record_sha256", None)
    actual = sha256(canonical_json(replay).encode()).hexdigest()
    if claimed != actual:
        raise RuntimeError("cache record hash mismatch")


def load_cache(path: Path, tasks: list[dict]) -> list[dict]:
    if not path.exists():
        return []
    records: list[dict] = []
    previous_hash = "GENESIS"
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) > len(tasks):
        raise RuntimeError("cache has more records than deterministic tasks")
    for index, line in enumerate(lines):
        if not line.strip():
            raise RuntimeError("cache contains a blank row")
        record = json.loads(line)
        validate_cache_record(record, tasks[index], previous_hash)
        records.append(record)
        previous_hash = record["record_sha256"]
    return records


def append_cache(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(canonical_json(record) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def parse_cell(row: dict) -> tuple:
    return compact.arb(row["full_f"]), compact.arb(row["full_f_prime"])


def point_in_cell(
    point: tuple[Fraction, Fraction],
    cell: tuple,
) -> bool:
    return cell[0].contains(fraction_to_arb(point[0])) and cell[
        1
    ].contains(fraction_to_arb(point[1]))


def build_open_chain(records: list[dict]) -> dict:
    leaves = [
        leaf
        for record in records
        for leaf in record["result"]["certified_leaves"]
    ]
    leaves.sort(key=lambda row: Fraction(row["x_low"]))
    if not leaves:
        raise RuntimeError("cannot build an empty phase chain")
    if Fraction(leaves[0]["x_low"]) != X_LOWER:
        raise RuntimeError("phase chain misses x=0")
    if Fraction(leaves[-1]["x_high"]) != X_UPPER:
        raise RuntimeError("phase chain misses x=246")
    for left, right in zip(leaves, leaves[1:]):
        if Fraction(left["x_high"]) != Fraction(right["x_low"]):
            raise RuntimeError("phase-cell x cover is not contiguous")

    cells = [parse_cell(leaf) for leaf in leaves]
    for index, cell in enumerate(cells):
        if interval_contains_zero(cell[0]) and interval_contains_zero(
            cell[1]
        ):
            raise RuntimeError(f"phase cell {index} contains the origin")

    first_f, first_f_prime = cells[0]
    if not interval_contains_zero(first_f_prime):
        raise RuntimeError("first cell does not contain the exact axis derivative")
    first_x = interval_midpoint_fraction(first_f)
    if first_x <= 0:
        if first_f.upper() <= 0:
            raise RuntimeError("axis-corner value is not positive")
        first_x = dyadic_fraction(first_f.upper() / 2)
    first = (first_x, Fraction(0))
    if not point_in_cell(first, cells[0]):
        raise RuntimeError("axis-corner witness misses the first cell")

    witnesses: list[tuple[Fraction, Fraction]] = [first]
    for index in range(len(cells) - 1):
        intersection_f = interval_intersection(
            cells[index][0],
            cells[index + 1][0],
        )
        intersection_f_prime = interval_intersection(
            cells[index][1],
            cells[index + 1][1],
        )
        witness = (
            interval_midpoint_fraction(intersection_f),
            interval_midpoint_fraction(intersection_f_prime),
        )
        if not point_in_cell(witness, cells[index]) or not point_in_cell(
            witness,
            cells[index + 1],
        ):
            raise RuntimeError(
                f"join witness {index + 1} misses an adjacent cell"
            )
        witnesses.append(witness)

    last = (
        interval_midpoint_fraction(cells[-1][0]),
        interval_midpoint_fraction(cells[-1][1]),
    )
    if not point_in_cell(last, cells[-1]):
        raise RuntimeError("right-corner witness misses the last cell")
    witnesses.append(last)

    return {
        "cell_count": len(cells),
        "witness_count": len(witnesses),
        "x_cover": [str(X_LOWER), str(X_UPPER)],
        "axis_witness_positive_real": bool(
            witnesses[0][0] > 0 and witnesses[0][1] == 0
        ),
        "open_positive_ray_crossing_count": open_ray_crossing(witnesses),
        "witnesses": [
            [str(coordinate) for coordinate in point]
            for point in witnesses
        ],
    }


def summarize(
    records: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    leaves = [
        leaf
        for record in records
        for leaf in record["result"].get("certified_leaves", [])
    ]
    unresolved = sum(
        record["result"].get("unresolved_boxes", 0)
        for record in records
    )
    complete = len(records) == len(tasks) and unresolved == 0
    chain = build_open_chain(records) if complete else None
    margins = [
        compact.arb(leaf["closest_point_norm_lower"])
        for leaf in leaves
    ]
    branch_counts = Counter(leaf["branch"] for leaf in leaves)
    return {
        "tasks_total": len(tasks),
        "tasks_completed": len(records),
        "certified_leaf_cells": len(leaves),
        "unresolved_boxes": unresolved,
        "branch_counts": dict(branch_counts),
        "maximum_subdivision_depth": max(
            (
                record["result"].get("maximum_depth", 0)
                for record in records
            ),
            default=0,
        ),
        "minimum_closest_point_norm_lower": (
            serialize_interval(min(margins)) if margins else None
        ),
        "complete_bottom_edge_certificate": complete,
        "stop_reason": stop_reason,
        "open_phase_chain": chain,
        "resource": resource,
    }


def write_artifact(
    out: Path,
    note: Path,
    cache: Path,
    records: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    summary = summarize(records, tasks, stop_reason, resource)
    complete = summary["complete_bottom_edge_certificate"]
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "rigorous complete Q208 bottom-edge phase-cell certificate"
            if complete
            else "resumable partial Q208 bottom-edge phase-cell computation"
        ),
        "proof_boundary": (
            "The complete status certifies only the fixed-time bottom path "
            "t=1/1040 on 0<=x<=246 for the proxy "
            "F=16*(1+x^4)*H and constructs an exact rational open phase "
            "chain. The top edge, transformed right-edge joins, cyclic "
            "winding, Q208, every later shell, Lambda<=0, RH, and the Clay "
            "prize remain open."
        ),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "builder_sha256": file_hash(Path(__file__).resolve()),
        "cache": str(cache.relative_to(REPO_ROOT)).replace("\\", "/"),
        "cache_sha256": file_hash(cache) if cache.exists() else None,
        "cache_last_record_sha256": (
            records[-1]["record_sha256"] if records else "GENESIS"
        ),
        "summary": summary,
        "records": records,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    note.write_text(render_note(artifact), encoding="utf-8")
    return artifact


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    chain = summary["open_phase_chain"]
    lines = [
        "# Newman Q208 Bottom Phase-Cell Certificate",
        "",
        f"Date: {DATE}",
        "",
        (
            "Status: rigorous complete fixed-time bottom certificate."
            if summary["complete_bottom_edge_certificate"]
            else "Status: resumable partial fixed-time computation."
        ),
        "This is not a proof of Q208, `Lambda<=0`, RH, or the Clay prize.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"work/rh_compute/results/{STEM}.jsonl",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Contract",
        "",
        "The certified proxy is",
        "",
        "```text",
        "F_t(x)=16*(1+x^4)*H_t(x),",
        "F_t'(x)=64*x^3*H_t(x)+16*(1+x^4)*H_t'(x).",
        "```",
        "",
        "Its triangular map from `(H,H')` has positive determinant and a",
        "canonical homotopy to the identity. It therefore has the same",
        "contact set and closed-boundary winding as the ordinary first jet,",
        "while remaining nonsingular at `x=0` and asymptotic to `J=16x^4H`.",
        "",
        "The retained transform uses six ordinary theta terms at the single",
        "time `t=1/1040`, Taylor order 24 in `x`, and the independently",
        "proved raw `E0(7),E1(7)` omitted-tail bars.",
        "",
        "## Progress",
        "",
        "```text",
        f"panels={summary['tasks_completed']}/{summary['tasks_total']}",
        f"certified cells={summary['certified_leaf_cells']}",
        f"unresolved={summary['unresolved_boxes']}",
        f"branches={summary['branch_counts']}",
        (
            "minimum closest-point norm lower="
            f"{summary['minimum_closest_point_norm_lower']}"
        ),
        f"maximum subdivision depth={summary['maximum_subdivision_depth']}",
        f"stop reason={summary['stop_reason']}",
        "```",
        "",
    ]
    if chain is not None:
        lines.extend(
            [
                "## Exact Open Phase Chain",
                "",
                "Every cell encloses the whole path segment and excludes the",
                "origin. Each adjacent pair has an exact dyadic rational",
                "intersection witness. The resulting open polygon has",
                "",
                "```text",
                f"cells={chain['cell_count']}",
                f"witnesses={chain['witness_count']}",
                (
                    "positive-ray crossing count="
                    f"{chain['open_positive_ray_crossing_count']}"
                ),
                "```",
                "",
                "This crossing count is an open-edge contribution, not a",
                "winding number. The right-corner witness must still be",
                "matched to the transformed right-edge cells, and the top",
                "edge must be traversed in reverse before the cyclic integer",
                "can be evaluated.",
                "",
            ]
        )
    lines.extend(
        [
            "## Scope",
            "",
            "A complete bottom certificate does not by itself certify Q208.",
            "The top phase chain, right-edge coordinate transfer, corner",
            "intersections, and one exact closed polygon winding remain open.",
            "",
        ]
    )
    return "\n".join(lines)


def baseline_samples() -> list[float]:
    return [
        float(psutil.cpu_percent(interval=1.0))
        for _ in range(BASELINE_SECONDS)
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument(
        "--max-panels",
        type=int,
        default=-1,
        help="Maximum new panels; -1 means unlimited and 0 initializes only.",
    )
    parser.add_argument(
        "--runtime-limit-seconds",
        type=float,
        default=DEFAULT_RUNTIME_LIMIT_SECONDS,
    )
    args = parser.parse_args()

    args.cache.parent.mkdir(parents=True, exist_ok=True)
    args.cache.touch(exist_ok=True)
    tasks = panel_tasks()
    records = load_cache(args.cache, tasks)
    priority_lowered = compact.request_below_normal_priority()
    baseline = baseline_samples()
    baseline_mean = sum(baseline) / len(baseline)
    resource = {
        "mode": "day",
        "worker_count": 1,
        "below_normal_priority_applied": priority_lowered,
        "baseline_samples": baseline,
        "baseline_mean": baseline_mean,
        "runtime_cpu_samples": [],
        "elapsed_seconds": 0.0,
    }
    if baseline_mean > BASELINE_LIMIT:
        write_artifact(
            args.out,
            args.note,
            args.cache,
            records,
            tasks,
            "deferred_baseline_cpu",
            resource,
        )
        print(
            "deferred Q208 bottom phase-cell certificate: "
            f"baseline mean {baseline_mean:.2f}% exceeds "
            f"{BASELINE_LIMIT:.2f}%"
        )
        return 0

    started = time.monotonic()
    certifier = FixedTimePanelCertifier()
    previous_hash = (
        records[-1]["record_sha256"] if records else "GENESIS"
    )
    high_cpu_streak = 0
    stop_reason = "complete"
    newly_completed = 0
    for task in tasks[len(records) :]:
        if (
            args.max_panels >= 0
            and newly_completed >= args.max_panels
        ):
            stop_reason = "max_panels"
            break
        if time.monotonic() - started >= args.runtime_limit_seconds:
            stop_reason = "runtime_limit"
            break
        result = certifier.certify_panel(task)
        record = cache_record(task, result, previous_hash)
        append_cache(args.cache, record)
        records.append(record)
        previous_hash = record["record_sha256"]
        newly_completed += 1
        cpu_sample = float(psutil.cpu_percent(interval=0.25))
        resource["runtime_cpu_samples"].append(cpu_sample)
        print(
            f"Q208 bottom {task['sequence']}/{len(tasks)} "
            f"x={task['x_low']}..{task['x_high']} "
            f"status={result['status']} "
            f"leaves={result['certified_leaf_boxes']} "
            f"depth={result['maximum_depth']} "
            f"in {result['elapsed_seconds']:.3f}s "
            f"cpu={cpu_sample:.1f}%",
            flush=True,
        )
        if result["status"] != "certified":
            stop_reason = "unresolved_panel"
            break
        if cpu_sample > CPU_PARK_THRESHOLD:
            high_cpu_streak += 1
        else:
            high_cpu_streak = 0
        if high_cpu_streak >= CPU_CONSECUTIVE_PARK_SAMPLES:
            stop_reason = "sustained_cpu_limit"
            break

    if len(records) < len(tasks) and stop_reason == "complete":
        stop_reason = "partial_existing_cache"
    resource["elapsed_seconds"] = time.monotonic() - started
    artifact = write_artifact(
        args.out,
        args.note,
        args.cache,
        records,
        tasks,
        stop_reason,
        resource,
    )
    summary = artifact["summary"]
    print(
        "stored Q208 bottom phase-cell certificate: "
        f"{summary['tasks_completed']}/{summary['tasks_total']} panels, "
        f"{summary['certified_leaf_cells']} cells, "
        f"{summary['unresolved_boxes']} unresolved, "
        f"complete={summary['complete_bottom_edge_certificate']}, "
        f"stop={summary['stop_reason']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
