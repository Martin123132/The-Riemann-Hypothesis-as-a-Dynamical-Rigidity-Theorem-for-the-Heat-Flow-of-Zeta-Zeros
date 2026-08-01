#!/usr/bin/env python3
"""Certify the finite Newman bridge with shared six-term Taylor panels."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import sys
import time

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_compact_transversality_interval_certificate as compact  # noqa: E402
import jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate as outer  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_forward_six_term_"
    "finite_bridge_interval_certificate"
)
RESULT_DIR = REPO_ROOT / "work" / "rh_compute" / "results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
TAIL_STEM = "jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate"
TAIL_SOURCE = RESULT_DIR / f"{TAIL_STEM}.json"
TAIL_BUILDER = SCRIPT_DIR / f"{TAIL_STEM}.py"
TAIL_CHECKER = SCRIPT_DIR / f"check_{TAIL_STEM}.py"
MODULAR_SOURCE = (
    RESULT_DIR / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
OUTER_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
COMPACT_SOURCE = compact.DEFAULT_OUT
Q31_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate.json"
)
DIAGONAL_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
WINDING_SOURCE = (
    RESULT_DIR / "jensen_window_pf_newman_first_jet_winding_gate.json"
)
BOUNDARY_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)

SCHEMA = "newman_theta_forward_six_term_finite_bridge_interval_v1"
DATE = "2026-07-25"
PRECISION_BITS = 352
RETAINED_TERMS = 6
FIRST_OMITTED = 7
TIME_GLOBAL_LOWER = Fraction(1, 1035)
TIME_GLOBAL_UPPER = Fraction(1, 5)
TIME_GLOBAL_CENTER = (TIME_GLOBAL_LOWER + TIME_GLOBAL_UPPER) / 2
TIME_GLOBAL_RADIUS = (TIME_GLOBAL_UPPER - TIME_GLOBAL_LOWER) / 2
LOW_STRIP_TIME_UPPER = Fraction(1, 155)
LOW_STRIP_X_LOWER = Fraction(38)
LOW_STRIP_X_UPPER = Fraction(69)
RIGHT_STRIP_X_LOWER = Fraction(69)
RIGHT_STRIP_X_UPPER = Fraction(245)
INITIAL_TIME_STEP = Fraction(1, 100)
INITIAL_X_STEP = Fraction(1, 2)
TAYLOR_TIME_ORDER = 30
TAYLOR_X_ORDER = 24
MAX_SUBDIVISION_DEPTH = 8
INTEGRATION_CUTOFF = Fraction(11, 5)
INTEGRATION_OPTIONS = {
    "deg_limit": 140,
    "eval_limit": 1_500_000,
    "depth_limit": 80,
    "use_heap": True,
    "abs_tol": "2^-290",
    "rel_tol": "2^-280",
}
ENDPOINT_DIGITS = 75
BASELINE_SECONDS = 6
BASELINE_ONE_WORKER_LIMIT = 55.0
CPU_PARK_THRESHOLD = 85.0
CPU_CONSECUTIVE_PARK_SAMPLES = 2
MAX_OSCILLATORY_MOMENT = (
    1 + 2 * TAYLOR_TIME_ORDER + TAYLOR_X_ORDER
)
MAX_ABSOLUTE_MOMENT = max(
    1 + 2 * (TAYLOR_TIME_ORDER + 1),
    1 + 2 * TAYLOR_TIME_ORDER + TAYLOR_X_ORDER + 1,
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_text(value: Fraction) -> str:
    return str(value)


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


def interval_ball(lower: Fraction, upper: Fraction) -> compact.arb:
    center = (lower + upper) / 2
    radius = (upper - lower) / 2
    return compact.interval_ball(center, radius)


def panel_tasks() -> list[dict]:
    tasks: list[dict] = []
    for lower, upper in fraction_range(
        LOW_STRIP_X_LOWER,
        LOW_STRIP_X_UPPER,
        INITIAL_X_STEP,
    ):
        tasks.append(
            {
                "region": "low_time_left_strip",
                "x_low": fraction_text(lower),
                "x_high": fraction_text(upper),
            }
        )
    for lower, upper in fraction_range(
        RIGHT_STRIP_X_LOWER,
        RIGHT_STRIP_X_UPPER,
        INITIAL_X_STEP,
    ):
        tasks.append(
            {
                "region": "full_time_right_strip",
                "x_low": fraction_text(lower),
                "x_high": fraction_text(upper),
            }
        )
    return tasks


def time_cells(region: str) -> list[tuple[Fraction, Fraction]]:
    upper = (
        LOW_STRIP_TIME_UPPER
        if region == "low_time_left_strip"
        else TIME_GLOBAL_UPPER
    )
    return fraction_range(
        TIME_GLOBAL_LOWER,
        upper,
        INITIAL_TIME_STEP,
    )


def task_id(task: dict) -> str:
    return f"{task['region']}:x{task['x_low']}..{task['x_high']}"


def serialize_positive(value: compact.arb) -> dict:
    if not value.is_finite() or value.lower() <= 0:
        raise RuntimeError(f"nonpositive or nonfinite enclosure: {value}")
    return {
        "enclosure": value.str(ENDPOINT_DIGITS, more=True),
        "upper": value.upper().str(ENDPOINT_DIGITS),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def integration_options() -> dict:
    options = dict(INTEGRATION_OPTIONS)
    for key in ("abs_tol", "rel_tol"):
        options[key] = compact.arb(2) ** (-int(options[key][3:]))
    return options


class SharedPanelCertifier:
    def __init__(self) -> None:
        compact.flint.ctx.prec = PRECISION_BITS
        self.pi = compact.acb.pi()
        self.zero = compact.acb(0)
        self.cutoff = compact.acb(
            compact.fraction_decimal(INTEGRATION_CUTOFF)
        )
        self.options = integration_options()
        self.polynomials = outer.polynomial_rows()
        self.cutoff_tails = self.build_cutoff_tails()
        self.absolute_moment_cache: dict[
            tuple[int, Fraction], compact.arb
        ] = {}
        self.models_built = 0
        self.oscillatory_integrals = 0
        self.absolute_integrals = 0
        tail = json.loads(TAIL_SOURCE.read_text(encoding="utf-8"))
        witnesses = {
            int(row["first_omitted"]): row
            for row in tail.get("witnesses", [])
        }
        if FIRST_OMITTED not in witnesses:
            raise RuntimeError("six-term tail witness is missing")
        witness = witnesses[FIRST_OMITTED]
        self.e0 = compact.arb(witness["E0"]["enclosure"])
        self.e1 = compact.arb(witness["E1"]["enclosure"])
        if self.e0.lower() <= 0 or self.e1.lower() <= 0:
            raise RuntimeError("six-term direct tail bars are not positive")

    def build_cutoff_tails(self) -> list[compact.arb]:
        compact.flint.ctx.prec = PRECISION_BITS + 48
        values: list[compact.arb] = []
        for moment in range(MAX_ABSOLUTE_MOMENT + 1):
            value, _ = outer.forward_outer_bound(
                0,
                moment,
                0,
                0,
                self.polynomials,
            )
            values.append(
                compact.arb(value.upper().str(ENDPOINT_DIGITS + 35))
            )
        compact.flint.ctx.prec = PRECISION_BITS
        return [
            compact.arb(value.str(ENDPOINT_DIGITS + 25, more=True))
            for value in values
        ]

    def phi(self, u: compact.acb, index: int) -> compact.acb:
        n2 = index * index
        phase = self.pi * n2 * (4 * u).exp()
        return (
            self.pi
            * n2
            * (5 * u).exp()
            * (2 * phase - 3)
            * (-phase).exp()
        )

    def retained_kernel(self, u: compact.acb) -> compact.acb:
        return sum(
            (
                self.phi(u, index)
                for index in range(1, RETAINED_TERMS + 1)
            ),
            compact.acb(0),
        )

    def cutoff_tail_complex(self, moment: int) -> compact.acb:
        radius = self.cutoff_tails[moment].upper()
        ball = compact.arb(0, radius.str(ENDPOINT_DIGITS + 20))
        return compact.acb(ball, ball)

    def oscillatory_integral(
        self,
        moment: int,
        frequency: Fraction,
    ) -> compact.acb:
        time_ball = compact.acb(
            compact.fraction_decimal(TIME_GLOBAL_CENTER)
        )
        frequency_ball = compact.acb(
            compact.fraction_decimal(frequency)
        )

        def integrand(
            u: compact.acb,
            _analytic: bool,
        ) -> compact.acb:
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_kernel(u)
                * (compact.acb(0, 1) * frequency_ball * u).exp()
            )

        value = compact.acb.integral(
            integrand,
            self.zero,
            self.cutoff,
            **self.options,
        )
        self.oscillatory_integrals += 1
        return value + self.cutoff_tail_complex(moment)

    def absolute_moment(
        self,
        moment: int,
        time_value: Fraction,
    ) -> compact.arb:
        key = (moment, time_value)
        cached = self.absolute_moment_cache.get(key)
        if cached is not None:
            return cached
        time_ball = compact.acb(
            compact.fraction_decimal(time_value)
        )

        def integrand(
            u: compact.acb,
            _analytic: bool,
        ) -> compact.acb:
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_kernel(u)
            )

        value = (
            compact.acb.integral(
                integrand,
                self.zero,
                self.cutoff,
                **self.options,
            ).real
            + self.cutoff_tails[moment].upper()
        )
        if value.lower() <= 0:
            raise RuntimeError(
                f"absolute moment is not positive at order {moment}"
            )
        self.absolute_integrals += 1
        self.absolute_moment_cache[key] = value
        return value

    def remainder(self, derivative_order: int, x_radius: Fraction) -> dict:
        time_radius = compact.arb(
            compact.fraction_decimal(TIME_GLOBAL_RADIUS)
        )
        x_radius_ball = compact.arb(
            compact.fraction_decimal(x_radius)
        )
        time_moment = (
            derivative_order + 2 * (TAYLOR_TIME_ORDER + 1)
        )
        time_remainder = (
            compact.arb_power(
                time_radius,
                TAYLOR_TIME_ORDER + 1,
            )
            / math.factorial(TAYLOR_TIME_ORDER + 1)
            * self.absolute_moment(
                time_moment,
                TIME_GLOBAL_UPPER,
            ).upper()
        )
        x_remainder = compact.arb(0)
        time_power = compact.arb(1)
        for time_order in range(TAYLOR_TIME_ORDER + 1):
            moment = (
                derivative_order
                + 2 * time_order
                + TAYLOR_X_ORDER
                + 1
            )
            x_remainder += (
                time_power
                / math.factorial(time_order)
                * compact.arb_power(
                    x_radius_ball,
                    TAYLOR_X_ORDER + 1,
                )
                / math.factorial(TAYLOR_X_ORDER + 1)
                * self.absolute_moment(
                    moment,
                    TIME_GLOBAL_CENTER,
                ).upper()
            )
            time_power *= time_radius
        total = (time_remainder + x_remainder).upper()
        return {
            "time": time_remainder,
            "frequency": x_remainder,
            "total": total,
        }

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
        coefficients: dict[int, list[list[compact.arb]]] = {}
        for derivative_order in (0, 1):
            rows: list[list[compact.arb]] = []
            for time_order in range(TAYLOR_TIME_ORDER + 1):
                row: list[compact.arb] = []
                for x_order in range(TAYLOR_X_ORDER + 1):
                    moment = (
                        derivative_order
                        + 2 * time_order
                        + x_order
                    )
                    coefficient = (
                        compact.acb(0, 1)
                        ** (derivative_order + x_order)
                        * integrals[moment]
                    ).real / (
                        math.factorial(time_order)
                        * math.factorial(x_order)
                    )
                    row.append(coefficient)
                rows.append(row)
            coefficients[derivative_order] = rows
        remainders = {
            derivative_order: self.remainder(
                derivative_order,
                radius,
            )
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
                value.rel_accuracy_bits() for value in integrals
            ),
        }

    def transform_box(
        self,
        model: dict,
        derivative_order: int,
        time_low: Fraction,
        time_high: Fraction,
        x_low: Fraction,
        x_high: Fraction,
    ) -> compact.arb:
        tau = interval_ball(
            time_low - TIME_GLOBAL_CENTER,
            time_high - TIME_GLOBAL_CENTER,
        )
        xi = interval_ball(
            x_low - model["x_center"],
            x_high - model["x_center"],
        )
        rows = model["coefficients"][derivative_order]
        value = compact.arb(0)
        for time_order in range(TAYLOR_TIME_ORDER, -1, -1):
            inner = compact.arb(0)
            for x_order in range(TAYLOR_X_ORDER, -1, -1):
                inner = inner * xi + rows[time_order][x_order]
            value = value * tau + inner
        radius = model["remainders"][derivative_order]["total"]
        return value + compact.arb(
            0,
            radius.str(ENDPOINT_DIGITS + 20),
        )

    def evaluate_box(
        self,
        model: dict,
        time_low: Fraction,
        time_high: Fraction,
        x_low: Fraction,
        x_high: Fraction,
        depth: int,
    ) -> dict:
        h_box = self.transform_box(
            model,
            0,
            time_low,
            time_high,
            x_low,
            x_high,
        )
        h_prime_box = self.transform_box(
            model,
            1,
            time_low,
            time_high,
            x_low,
            x_high,
        )
        x_box = interval_ball(x_low, x_high)
        j_box = 16 * x_box**4 * h_box
        j_prime_box = (
            64 * x_box**3 * h_box
            + 16 * x_box**4 * h_prime_box
        )

        x_ceiling = compact.arb(
            compact.fraction_decimal(x_high)
        )
        tail_value = (16 * x_ceiling**4 * self.e0).upper()
        tail_derivative = (
            64 * x_ceiling**3 * self.e0
            + 16 * x_ceiling**4 * self.e1
        ).upper()
        value_lower = j_box.abs_lower()
        derivative_lower = j_prime_box.abs_lower()
        value_ratio = value_lower / tail_value
        derivative_ratio = derivative_lower / tail_derivative
        value_pass = value_lower > tail_value
        derivative_pass = derivative_lower > tail_derivative
        if value_pass and (
            not derivative_pass or value_ratio > derivative_ratio
        ):
            branch = "value"
            ratio = value_ratio
        elif derivative_pass:
            branch = "derivative"
            ratio = derivative_ratio
        else:
            branch = None
            ratio = max(value_ratio, derivative_ratio)
        return {
            "depth": depth,
            "t_low": str(time_low),
            "t_high": str(time_high),
            "x_low": str(x_low),
            "x_high": str(x_high),
            "j_retained_lower": j_box.lower().str(ENDPOINT_DIGITS),
            "j_retained_upper": j_box.upper().str(ENDPOINT_DIGITS),
            "j_retained_prime_lower": j_prime_box.lower().str(
                ENDPOINT_DIGITS
            ),
            "j_retained_prime_upper": j_prime_box.upper().str(
                ENDPOINT_DIGITS
            ),
            "tail_value_upper": tail_value.str(ENDPOINT_DIGITS),
            "tail_derivative_upper": tail_derivative.str(
                ENDPOINT_DIGITS
            ),
            "value_ratio_lower": value_ratio.lower().str(
                ENDPOINT_DIGITS
            ),
            "derivative_ratio_lower": derivative_ratio.lower().str(
                ENDPOINT_DIGITS
            ),
            "certified": bool(value_pass or derivative_pass),
            "branch": branch,
            "certified_ratio_lower": ratio.lower().str(
                ENDPOINT_DIGITS
            ),
            "full_value_negative": bool(
                j_box.upper() < -tail_value
            ),
            "full_value_positive": bool(
                j_box.lower() > tail_value
            ),
        }


def split_bounds(
    bounds: tuple[Fraction, Fraction, Fraction, Fraction],
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    time_low, time_high, x_low, x_high = bounds
    normalized_time = (time_high - time_low) / INITIAL_TIME_STEP
    normalized_x = (x_high - x_low) / INITIAL_X_STEP
    if normalized_x >= normalized_time:
        middle = (x_low + x_high) / 2
        return [
            (time_low, time_high, x_low, middle),
            (time_low, time_high, middle, x_high),
        ]
    middle = (time_low + time_high) / 2
    return [
        (time_low, middle, x_low, x_high),
        (middle, time_high, x_low, x_high),
    ]


def certify_panel(
    certifier: SharedPanelCertifier,
    task: dict,
) -> dict:
    x_low = Fraction(task["x_low"])
    x_high = Fraction(task["x_high"])
    started = time.monotonic()
    model = certifier.build_model(x_low, x_high)
    queue = [
        ((time_low, time_high, x_low, x_high), 0)
        for time_low, time_high in time_cells(task["region"])
    ]
    cursor = 0
    certified_leaves: list[dict] = []
    unresolved: list[dict] = []
    subdivisions = 0
    evaluations = 0
    while cursor < len(queue):
        bounds, depth = queue[cursor]
        cursor += 1
        evaluation = certifier.evaluate_box(
            model,
            *bounds,
            depth,
        )
        evaluations += 1
        if evaluation["certified"]:
            certified_leaves.append(evaluation)
        elif depth < MAX_SUBDIVISION_DEPTH:
            queue.extend(
                (child, depth + 1)
                for child in split_bounds(bounds)
            )
            subdivisions += 1
        else:
            unresolved.append(evaluation)

    ratios = [
        compact.arb(row["certified_ratio_lower"])
        for row in certified_leaves
    ]
    minimum = (
        min(ratios, key=lambda value: float(value.lower()))
        if ratios
        else None
    )
    return {
        "status": "certified" if not unresolved else "unresolved",
        "initial_time_cells": len(time_cells(task["region"])),
        "evaluated_boxes": evaluations,
        "certified_leaf_boxes": len(certified_leaves),
        "subdivisions": subdivisions,
        "unresolved_boxes": len(unresolved),
        "maximum_depth": max(
            (
                row["depth"]
                for row in certified_leaves + unresolved
            ),
            default=0,
        ),
        "minimum_certified_ratio_lower": (
            minimum.lower().str(ENDPOINT_DIGITS)
            if minimum is not None
            else None
        ),
        "branch_counts": {
            "value": sum(
                row["branch"] == "value"
                for row in certified_leaves
            ),
            "derivative": sum(
                row["branch"] == "derivative"
                for row in certified_leaves
            ),
        },
        "model": {
            "x_center": str(model["x_center"]),
            "x_radius": str(model["x_radius"]),
            "oscillatory_moments": MAX_OSCILLATORY_MOMENT + 1,
            "minimum_integral_accuracy_bits": (
                model["minimum_integral_accuracy_bits"]
            ),
            "h_remainder": {
                key: serialize_positive(value)
                for key, value in model["remainders"][0].items()
            },
            "h_prime_remainder": {
                key: serialize_positive(value)
                for key, value in model["remainders"][1].items()
            },
        },
        "certified_leaves": certified_leaves,
        "unresolved_leaves": unresolved,
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }


def config_payload() -> dict:
    return {
        "schema": SCHEMA,
        "implementation_sha256": file_hash(Path(__file__)),
        "precision_bits": PRECISION_BITS,
        "retained_terms": RETAINED_TERMS,
        "first_omitted": FIRST_OMITTED,
        "time_global": [
            str(TIME_GLOBAL_LOWER),
            str(TIME_GLOBAL_UPPER),
        ],
        "low_strip": {
            "time": [
                str(TIME_GLOBAL_LOWER),
                str(LOW_STRIP_TIME_UPPER),
            ],
            "x": [
                str(LOW_STRIP_X_LOWER),
                str(LOW_STRIP_X_UPPER),
            ],
        },
        "right_strip": {
            "time": [
                str(TIME_GLOBAL_LOWER),
                str(TIME_GLOBAL_UPPER),
            ],
            "x": [
                str(RIGHT_STRIP_X_LOWER),
                str(RIGHT_STRIP_X_UPPER),
            ],
        },
        "initial_time_step": str(INITIAL_TIME_STEP),
        "initial_x_step": str(INITIAL_X_STEP),
        "taylor_time_order": TAYLOR_TIME_ORDER,
        "taylor_x_order": TAYLOR_X_ORDER,
        "max_subdivision_depth": MAX_SUBDIVISION_DEPTH,
        "integration_cutoff": str(INTEGRATION_CUTOFF),
        "integration_options": INTEGRATION_OPTIONS,
        "max_oscillatory_moment": MAX_OSCILLATORY_MOMENT,
        "max_absolute_moment": MAX_ABSOLUTE_MOMENT,
        "tail_source_sha256": file_hash(TAIL_SOURCE),
        "outer_source_sha256": file_hash(OUTER_SOURCE),
        "task_ids": [task_id(task) for task in panel_tasks()],
    }


def config_hash() -> str:
    return sha256(
        canonical_json(config_payload()).encode("utf-8")
    ).hexdigest()


def load_cache(path: Path) -> tuple[list[dict], dict[str, dict]]:
    if not path.exists():
        return [], {}
    records: list[dict] = []
    by_id: dict[str, dict] = {}
    previous = "0" * 64
    expected_config = config_hash()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.endswith("\n"):
                raise RuntimeError(
                    f"cache line {line_number} lacks a terminator"
                )
            record = json.loads(line)
            stored_hash = record.pop("row_sha256", None)
            actual_hash = sha256(
                canonical_json(record).encode("utf-8")
            ).hexdigest()
            record["row_sha256"] = stored_hash
            if stored_hash != actual_hash:
                raise RuntimeError(
                    f"cache hash mismatch at line {line_number}"
                )
            if record.get("sequence") != line_number:
                raise RuntimeError(
                    f"cache sequence mismatch at line {line_number}"
                )
            if record.get("previous_row_sha256") != previous:
                raise RuntimeError(
                    f"cache chain mismatch at line {line_number}"
                )
            if record.get("config_sha256") != expected_config:
                raise RuntimeError(
                    f"cache config mismatch at line {line_number}"
                )
            identifier = task_id(record["task"])
            if identifier in by_id:
                raise RuntimeError(
                    f"duplicate cache task: {identifier}"
                )
            records.append(record)
            by_id[identifier] = record
            previous = stored_hash
    return records, by_id


def append_record(
    path: Path,
    records: list[dict],
    task: dict,
    result: dict,
) -> dict:
    record = {
        "schema": SCHEMA,
        "sequence": len(records) + 1,
        "config_sha256": config_hash(),
        "task": task,
        "result": result,
        "previous_row_sha256": (
            records[-1]["row_sha256"] if records else "0" * 64
        ),
    }
    record["row_sha256"] = sha256(
        canonical_json(record).encode("utf-8")
    ).hexdigest()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(canonical_json(record) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    records.append(record)
    return record


def baseline_samples() -> list[float]:
    return [
        float(psutil.cpu_percent(interval=1.0))
        for _ in range(BASELINE_SECONDS)
    ]


def run_tasks(
    cache: Path,
    max_new_panels: int | None,
    max_seconds: float | None,
    stop_file: Path | None,
    skip_baseline: bool,
) -> tuple[
    list[dict],
    dict[str, dict],
    int,
    list[float],
    bool,
    list[float],
    bool,
    dict,
]:
    records, by_id = load_cache(cache)
    pending = [
        task
        for task in panel_tasks()
        if task_id(task) not in by_id
    ]
    baseline = [] if skip_baseline or not pending else baseline_samples()
    baseline_deferred = bool(
        baseline
        and sum(baseline) / len(baseline)
        > BASELINE_ONE_WORKER_LIMIT
    )
    if baseline_deferred:
        return (
            records,
            by_id,
            0,
            [],
            False,
            baseline,
            True,
            {},
        )

    certifier: SharedPanelCertifier | None = None
    completed = 0
    cpu_samples: list[float] = []
    consecutive_high = 0
    resource_parked = False
    started = time.monotonic()
    for task in pending:
        if stop_file is not None and stop_file.exists():
            break
        if max_new_panels is not None and completed >= max_new_panels:
            break
        if (
            max_seconds is not None
            and completed > 0
            and time.monotonic() - started >= max_seconds
        ):
            break
        if certifier is None:
            certifier = SharedPanelCertifier()
        panel_started = time.monotonic()
        result = certify_panel(certifier, task)
        record = append_record(
            cache,
            records,
            task,
            result,
        )
        by_id[task_id(task)] = record
        completed += 1
        print(
            f"completed {record['sequence']}/{len(panel_tasks())} "
            f"{task_id(task)} status={result['status']} "
            f"leaves={result['certified_leaf_boxes']} "
            f"unresolved={result['unresolved_boxes']} "
            f"in {time.monotonic() - panel_started:.3f}s",
            flush=True,
        )
        sample = float(psutil.cpu_percent(interval=0.2))
        cpu_samples.append(sample)
        if sample > CPU_PARK_THRESHOLD:
            consecutive_high += 1
        else:
            consecutive_high = 0
        if consecutive_high >= CPU_CONSECUTIVE_PARK_SAMPLES:
            resource_parked = True
            break
    counters = (
        {
            "models_built": certifier.models_built,
            "oscillatory_integrals": certifier.oscillatory_integrals,
            "absolute_integrals": certifier.absolute_integrals,
        }
        if certifier is not None
        else {}
    )
    return (
        records,
        by_id,
        completed,
        cpu_samples,
        resource_parked,
        baseline,
        False,
        counters,
    )


def source_audit() -> dict:
    paths = {
        "six_term_tail": TAIL_SOURCE,
        "modular_partition": MODULAR_SOURCE,
        "outer_tail": OUTER_SOURCE,
        "compact_core": COMPACT_SOURCE,
        "q31": Q31_SOURCE,
        "diagonal_exhaustion": DIAGONAL_SOURCE,
        "winding": WINDING_SOURCE,
        "boundary": BOUNDARY_SOURCE,
    }
    sources = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in paths.items()
    }
    text = {
        key: json.dumps(source, sort_keys=True)
        for key, source in sources.items()
    }
    markers = {
        "six_term_tail": (
            "38<=x<=245",
            "10^-11*exp(-pi*x/8)",
        ),
        "modular_partition": (
            "Phi(u)=sum_(n>=1)phi_n(u)=Phi(-u)",
            "phi_n(u)>0",
        ),
        "outer_tail": ("u>11/5",),
        "compact_core": ("|x|<=38",),
        "q31": (
            "Q_31=[1/155,1/4]x[0,69]",
            '"theorem_ready": true',
        ),
        "diagonal_exhaustion": (
            "R_j=38+j",
            "delta_j=1/(5j)",
        ),
        "winding": (
            "Q_j=[1/(5j),1/4]x[0,38+j]",
            "wind(Z(partial Q_j),0)=0",
        ),
        "boundary": ("Lambda<=1/5",),
    }
    for key, required in markers.items():
        for marker in required:
            if marker not in text[key]:
                raise RuntimeError(
                    f"{key} source marker missing: {marker}"
                )
    return {
        key: {
            "kind": source["kind"],
            "path": str(paths[key].relative_to(REPO_ROOT)).replace(
                "\\",
                "/",
            ),
            "sha256": file_hash(paths[key]),
        }
        for key, source in sources.items()
    } | {
        "tail_builder_sha256": file_hash(TAIL_BUILDER),
        "tail_checker_sha256": file_hash(TAIL_CHECKER),
    }


def build_artifact(
    cache: Path,
    records: list[dict],
    by_id: dict[str, dict],
    cpu_samples: list[float],
    resource_parked: bool,
    baseline: list[float],
    baseline_deferred: bool,
    counters: dict,
) -> dict:
    tasks = panel_tasks()
    missing = [
        task_id(task)
        for task in tasks
        if task_id(task) not in by_id
    ]
    results = [record["result"] for record in records]
    unresolved = [
        result
        for result in results
        if result["status"] != "certified"
    ]
    complete = not missing
    theorem_ready = complete and not unresolved
    certified_leaves = [
        leaf
        for result in results
        for leaf in result["certified_leaves"]
    ]
    ratios = [
        compact.arb(leaf["certified_ratio_lower"])
        for leaf in certified_leaves
    ]
    minimum = (
        min(ratios, key=lambda value: float(value.lower()))
        if ratios
        else None
    )
    finite_readiness = (
        "proved" if theorem_ready else "not_ready_to_apply"
    )
    rows = [
        GateRow(
            id="ntfsfbic_01_six_term_forward_partition",
            role="exact_input",
            readiness="proved",
            claim=(
                "The retained transform is the first six ordinary positive-"
                "half-line theta summands."
            ),
            formula=(
                "S_(6,t)^F=sum_(n=1)^6 integral_0^infinity "
                "exp(tu^2)phi_n(u)cos(xu)du"
            ),
            proof_boundary=(
                "No termwise evenness or integration-by-parts cancellation "
                "is assumed."
            ),
        ),
        GateRow(
            id="ntfsfbic_02_direct_six_term_tail",
            role="rigorous_input",
            readiness="proved",
            claim=(
                "The direct omitted value and first-moment tails are below "
                "10^-11 of the gamma scale on the complete bridge."
            ),
            formula=(
                "|J-J_6^F|,|J'-(J_6^F)'|"
                "<10^-11*exp(-pi*x/8)"
            ),
            proof_boundary="Uses the independently checked six-term tail gate.",
        ),
        GateRow(
            id="ntfsfbic_03_shared_taylor_panel",
            role="exact_inequality",
            readiness="proved",
            claim=(
                "One directed bivariate Taylor model per half-unit x panel "
                "covers every time cell with explicit absolute remainders."
            ),
            formula=(
                "A=30, B=24; "
                "R_t<=T^(A+1)M_(d+2A+2)/(A+1)!; "
                "R_x<=sum_a T^a X^(B+1)M_(d+2a+B+1)/(a!(B+1)!)"
            ),
            proof_boundary=(
                "All moments through order 86 include an analytic "
                "u>11/5 tail."
            ),
        ),
        GateRow(
            id="ntfsfbic_04_finite_bridge_cover",
            role="directed_rounding_certificate",
            readiness=finite_readiness,
            claim=(
                "Every leaf in the low-time left strip and full-time right "
                "strip satisfies the retained value-or-derivative "
                "disjunction."
            ),
            formula=(
                "[1/1035,1/155]x[38,69] union "
                "[1/1035,1/5]x[69,245]"
            ),
            proof_boundary=(
                "Claim activates only when all 414 panel records are present "
                "and no unresolved leaf remains."
            ),
        ),
        GateRow(
            id="ntfsfbic_05_q207_no_contact",
            role="finite_diagonal_theorem",
            readiness=finite_readiness,
            claim=(
                "The compact core, Q31, finite bridge, and t>1/5 theorem "
                "exclude every first-jet contact in Q207."
            ),
            formula="Q_207=[1/1035,1/4]x[0,245]",
            proof_boundary="This is a finite diagonal theorem, not cofinality.",
        ),
        GateRow(
            id="ntfsfbic_06_q207_winding",
            role="finite_diagonal_theorem",
            readiness=finite_readiness,
            claim=(
                "Boundary nonvanishing and the contact-index theorem give "
                "zero first-jet winding on Q207."
            ),
            formula="wind((H+iH_x)(partial Q_207),0)=0",
            proof_boundary="Finite winding does not prove Lambda<=0.",
        ),
        GateRow(
            id="ntfsfbic_07_cofinal_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "Replace the finite x<=245 cover by retained adaptive "
                "first-jet separation on every cofinal transition cell."
            ),
            formula=(
                "x>=245; N(x)=ceil((1+x)^(3/4)); "
                "max(|J_N|,|J_N'|)>quarter-gamma tail"
            ),
            proof_boundary=(
                "No cofinal retained theorem, strict Laguerre positivity, "
                "Lambda<=0, RH, or Clay-prize proof is supplied."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "complete rigorous finite Q207 six-term interval theorem"
            if theorem_ready
            else "partial resumable six-term finite-bridge interval certificate"
        ),
        "config": config_payload(),
        "config_sha256": config_hash(),
        "cache": {
            "path": str(cache.relative_to(REPO_ROOT)).replace(
                "\\",
                "/",
            ),
            "records": len(records),
            "expected_records": len(tasks),
            "last_row_sha256": (
                records[-1]["row_sha256"] if records else None
            ),
            "complete": complete,
            "missing_count": len(missing),
            "missing_task_ids": missing,
            "resource_parked_last_run": resource_parked,
            "cpu_samples_last_run": cpu_samples,
            "baseline_samples_last_run": baseline,
            "baseline_deferred_last_run": baseline_deferred,
            "run_counters_last_run": counters,
        },
        "summary": {
            "theorem_ready": theorem_ready,
            "certified_panels": sum(
                result["status"] == "certified"
                for result in results
            ),
            "unresolved_panels": len(unresolved),
            "initial_time_cells": sum(
                result["initial_time_cells"]
                for result in results
            ),
            "evaluated_boxes": sum(
                result["evaluated_boxes"]
                for result in results
            ),
            "certified_leaf_boxes": len(certified_leaves),
            "subdivisions": sum(
                result["subdivisions"]
                for result in results
            ),
            "branch_counts": {
                "value": sum(
                    result["branch_counts"]["value"]
                    for result in results
                ),
                "derivative": sum(
                    result["branch_counts"]["derivative"]
                    for result in results
                ),
            },
            "minimum_certified_ratio_lower": (
                minimum.lower().str(ENDPOINT_DIGITS)
                if minimum is not None
                else None
            ),
        },
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": (
            "Even when complete, this artifact proves only the finite "
            "Q207 rectangle through x=245. It does not prove retained "
            "adaptive first-jet separation on the cofinal transition cells, "
            "strict Laguerre positivity, Lambda<=0 or RH, and does not "
            "supply a Clay-prize conclusion."
        ),
    }


def render_note(artifact: dict) -> str:
    cache = artifact["cache"]
    summary = artifact["summary"]
    lines = [
        "# Newman Theta Forward Six-Term Finite-Bridge Interval Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: "
        + (
            "complete rigorous finite Q207 theorem."
            if summary["theorem_ready"]
            else "partial resumable finite-bridge certificate."
        ),
        "This is not a cofinal theorem and not a proof of `Lambda<=0` or RH.",
        "",
        "## Domain",
        "",
        "```text",
        "low-time left strip: [1/1035,1/155] x [38,69]",
        "full-time right strip: [1/1035,1/5] x [69,245]",
        "retained ordinary theta terms: n=1,...,6",
        "```",
        "",
        "The completed Q31 theorem supplies the complementary",
        "`[1/155,1/4] x [0,69]` region, while the compact theorem",
        "supplies `x<=38`. Together with the published `Lambda<=1/5`",
        "boundary input, a complete cache promotes the finite rectangle",
        "",
        "```text",
        "Q_207=[1/1035,1/4] x [0,245].",
        "```",
        "",
        "## Shared Models",
        "",
        "Each half-unit x panel uses one 352-bit model centered at",
        "`t=(1/1035+1/5)/2`, shared by every time cell:",
        "",
        "```text",
        "time Taylor order: 30",
        "frequency Taylor order: 24",
        "oscillatory moments per panel: 86",
        "maximum absolute moment: 86",
        "cutoff: u=11/5 with analytic forward tail",
        "```",
        "",
        "The direct six-term omitted-tail theorem is used without",
        "integration by parts:",
        "",
        "```text",
        "|J-J_6^F|<=16*x^4*E0(7),",
        "|J'-(J_6^F)'|<=64*x^3*E0(7)+16*x^4*E1(7).",
        "```",
        "",
        "## Progress",
        "",
        f"Panels: `{cache['records']}/{cache['expected_records']}`.",
        f"Certified panels: `{summary['certified_panels']}`.",
        f"Unresolved panels: `{summary['unresolved_panels']}`.",
        f"Initial time cells represented: `{summary['initial_time_cells']}`.",
        f"Certified leaves: `{summary['certified_leaf_boxes']}`.",
        f"Subdivisions: `{summary['subdivisions']}`.",
        "Minimum certified ratio: "
        f"`{summary['minimum_certified_ratio_lower']}`.",
        "",
        "## Proof Boundary",
        "",
        "Completion proves only the finite Q207 diagonal. The surviving",
        "Newman obligation is retained adaptive first-jet separation on",
        "every cofinal transition cell for `x>=245`.",
        "",
    ]
    return "\n".join(lines)


def write_artifact(
    artifact: dict,
    out: Path,
    note: Path,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--max-new-panels", type=int)
    parser.add_argument("--max-seconds", type=float)
    parser.add_argument("--stop-file", type=Path)
    parser.add_argument("--skip-baseline", action="store_true")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()

    priority_lowered = compact.request_below_normal_priority()
    (
        records,
        by_id,
        completed,
        cpu_samples,
        resource_parked,
        baseline,
        baseline_deferred,
        counters,
    ) = run_tasks(
        args.cache,
        args.max_new_panels,
        args.max_seconds,
        args.stop_file,
        args.skip_baseline,
    )
    artifact = build_artifact(
        args.cache,
        records,
        by_id,
        cpu_samples,
        resource_parked,
        baseline,
        baseline_deferred,
        counters,
    )
    artifact["below_normal_priority_applied"] = priority_lowered
    write_artifact(artifact, args.out, args.note)
    print(
        f"finite-bridge progress: new={completed} "
        f"cache={artifact['cache']['records']}/"
        f"{artifact['cache']['expected_records']} "
        f"unresolved={artifact['summary']['unresolved_panels']} "
        f"baseline_deferred={baseline_deferred} "
        f"resource_parked={resource_parked}",
        flush=True,
    )
    if args.require_complete and not artifact["summary"]["theorem_ready"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
