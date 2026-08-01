#!/usr/bin/env python3
"""Certify the Q207-to-Q208 bottom collar by heat-jet transport."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import time

import psutil
import sympy as sp

import jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate as bridge
from jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate import (
    forward_tail,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_certificate"
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
REFERENCE_STEM = (
    "jensen_window_pf_newman_q208_bottom_phase_cell_certificate"
)
REFERENCE_CACHE = RESULT_DIR / f"{REFERENCE_STEM}.jsonl"
REFERENCE_RESULT = RESULT_DIR / f"{REFERENCE_STEM}.json"
SUCCESSOR_RESULT = (
    RESULT_DIR
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)

SCHEMA = "newman_q207_q208_adiabatic_bottom_collar_v1"
DATE = "2026-07-25"
TIME_LOW = Fraction(1, 1040)
TIME_HIGH = Fraction(1, 1035)
TIME_CENTER = (TIME_LOW + TIME_HIGH) / 2
TIME_RADIUS = (TIME_HIGH - TIME_LOW) / 2
TIME_STEP = TIME_HIGH - TIME_LOW
X_LOWER = Fraction(0)
X_UPPER = Fraction(245)
X_STEP = Fraction(1, 2)
TIME_TAYLOR_ORDER = 2
X_TAYLOR_ORDER = 24
DERIVATIVE_ORDERS = (2, 3)
MAX_OSCILLATORY_MOMENT = (
    max(DERIVATIVE_ORDERS)
    + 2 * TIME_TAYLOR_ORDER
    + X_TAYLOR_ORDER
)
ENDPOINT_DIGITS = 105
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


def fraction_to_arb(value: Fraction):
    return compact.arb(value.numerator) / value.denominator


def serialize_interval(value) -> str:
    if not value.is_finite():
        raise RuntimeError(f"nonfinite interval: {value}")
    return value.str(ENDPOINT_DIGITS, more=True)


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
            fraction_range(X_LOWER, X_UPPER, X_STEP),
            start=1,
        )
    ]


def load_reference_leaves() -> list[dict]:
    records = [
        json.loads(line)
        for line in REFERENCE_CACHE.read_text(
            encoding="utf-8"
        ).splitlines()
        if line.strip()
    ]
    leaves = [
        leaf
        for record in records
        for leaf in record["result"]["certified_leaves"]
        if Fraction(leaf["x_high"]) <= X_UPPER
    ]
    leaves.sort(key=lambda row: Fraction(row["x_low"]))
    tasks = panel_tasks()
    if len(leaves) != len(tasks):
        raise RuntimeError("reference phase-cell count mismatch")
    for task, leaf in zip(tasks, leaves, strict=True):
        if (
            task["x_low"] != leaf["x_low"]
            or task["x_high"] != leaf["x_high"]
        ):
            raise RuntimeError("reference phase-cell partition mismatch")
    return leaves


def direct_tail(moment: int):
    return (
        5
        * forward_tail(
            moment,
            kernel_power=2,
            beta=sp.Rational(9, 4),
            first=bridge.FIRST_OMITTED,
        )
    ).upper()


def source_hashes() -> dict:
    return {
        "bridge_builder": file_hash(Path(bridge.__file__).resolve()),
        "six_term_tail": file_hash(bridge.TAIL_SOURCE),
        "outer_tail": file_hash(bridge.OUTER_SOURCE),
        "reference_cache": file_hash(REFERENCE_CACHE),
        "reference_result": file_hash(REFERENCE_RESULT),
        "successor_lemma": file_hash(SUCCESSOR_RESULT),
    }


def contract_payload() -> dict:
    return {
        "schema": SCHEMA,
        "time_collar": [str(TIME_LOW), str(TIME_HIGH)],
        "time_step": str(TIME_STEP),
        "x_domain": [str(X_LOWER), str(X_UPPER)],
        "x_step": str(X_STEP),
        "time_taylor_order": TIME_TAYLOR_ORDER,
        "x_taylor_order": X_TAYLOR_ORDER,
        "derivative_orders": list(DERIVATIVE_ORDERS),
        "precision_bits": bridge.PRECISION_BITS,
        "retained_terms": bridge.RETAINED_TERMS,
        "first_omitted": bridge.FIRST_OMITTED,
        "serialization_policy": (
            "form every scalar transport bound from reparsed outward "
            "serialized interval enclosures"
        ),
        "proxy": "F=16*(1+x^4)*H",
        "transport": (
            "F_t=-16*(1+x^4)*H_xx; "
            "F_xt=-64*x^3*H_xx-16*(1+x^4)*H_xxx"
        ),
        "strict_gate": "delta*sup||(F_t,F_xt)||<dist(0,C_k)",
        "source_sha256": source_hashes(),
    }


def contract_hash() -> str:
    return sha256(canonical_json(contract_payload()).encode()).hexdigest()


class CollarPanelCertifier:
    def __init__(self) -> None:
        compact.flint.ctx.prec = bridge.PRECISION_BITS
        self.base = bridge.SharedPanelCertifier()
        self.time_ball = compact.acb(
            compact.fraction_decimal(TIME_CENTER)
        )
        self.tail = {
            derivative_order: direct_tail(derivative_order)
            for derivative_order in DERIVATIVE_ORDERS
        }
        if any(value.lower() <= 0 for value in self.tail.values()):
            raise RuntimeError("direct derivative tail is not positive")
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

    def remainder(
        self,
        derivative_order: int,
        x_radius: Fraction,
    ) -> dict:
        time_radius = fraction_to_arb(TIME_RADIUS)
        x_radius_ball = fraction_to_arb(x_radius)
        time_moment = (
            derivative_order + 2 * (TIME_TAYLOR_ORDER + 1)
        )
        time_remainder = (
            compact.arb_power(
                time_radius,
                TIME_TAYLOR_ORDER + 1,
            )
            / math.factorial(TIME_TAYLOR_ORDER + 1)
            * self.base.absolute_moment(
                time_moment,
                TIME_HIGH,
            ).upper()
        )
        x_remainder = compact.arb(0)
        time_power = compact.arb(1)
        for time_order in range(TIME_TAYLOR_ORDER + 1):
            moment = (
                derivative_order
                + 2 * time_order
                + X_TAYLOR_ORDER
                + 1
            )
            x_remainder += (
                time_power
                / math.factorial(time_order)
                * compact.arb_power(
                    x_radius_ball,
                    X_TAYLOR_ORDER + 1,
                )
                / math.factorial(X_TAYLOR_ORDER + 1)
                * self.base.absolute_moment(
                    moment,
                    TIME_CENTER,
                ).upper()
            )
            time_power *= time_radius
        return {
            "time": time_remainder.upper(),
            "frequency": x_remainder.upper(),
            "total": (time_remainder + x_remainder).upper(),
        }

    def build_model(
        self,
        x_low: Fraction,
        x_high: Fraction,
    ) -> dict:
        x_center = (x_low + x_high) / 2
        x_radius = (x_high - x_low) / 2
        integrals = [
            self.oscillatory_integral(moment, x_center)
            for moment in range(MAX_OSCILLATORY_MOMENT + 1)
        ]
        coefficients: dict[int, list[list]] = {}
        for derivative_order in DERIVATIVE_ORDERS:
            rows: list[list] = []
            for time_order in range(TIME_TAYLOR_ORDER + 1):
                row: list = []
                for x_order in range(X_TAYLOR_ORDER + 1):
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
                x_radius,
            )
            for derivative_order in DERIVATIVE_ORDERS
        }
        self.models_built += 1
        return {
            "x_center": x_center,
            "x_radius": x_radius,
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
        tau = bridge.interval_ball(
            TIME_LOW - TIME_CENTER,
            TIME_HIGH - TIME_CENTER,
        )
        xi = bridge.interval_ball(
            x_low - model["x_center"],
            x_high - model["x_center"],
        )
        rows = model["coefficients"][derivative_order]
        value = compact.arb(0)
        for time_order in range(TIME_TAYLOR_ORDER, -1, -1):
            inner = compact.arb(0)
            for x_order in range(X_TAYLOR_ORDER, -1, -1):
                inner = (
                    inner * xi + rows[time_order][x_order]
                )
            value = value * tau + inner
        return value + compact.arb(
            0,
            model["remainders"][derivative_order][
                "total"
            ].str(ENDPOINT_DIGITS, more=True),
        )

    def certify_panel(self, task: dict, reference: dict) -> dict:
        x_low = Fraction(task["x_low"])
        x_high = Fraction(task["x_high"])
        started = time.monotonic()
        model = self.build_model(x_low, x_high)
        h_xx_retained = self.transform_box(
            model,
            2,
            x_low,
            x_high,
        )
        h_xxx_retained = self.transform_box(
            model,
            3,
            x_low,
            x_high,
        )
        h_xx_raw = h_xx_retained + compact.arb(
            0,
            self.tail[2].str(ENDPOINT_DIGITS, more=True),
        )
        h_xxx_raw = h_xxx_retained + compact.arb(
            0,
            self.tail[3].str(ENDPOINT_DIGITS, more=True),
        )
        h_xx_text = serialize_interval(h_xx_raw)
        h_xxx_text = serialize_interval(h_xxx_raw)
        # Ratios are formed from the outward serialized enclosures so an
        # independent checker can reconstruct every inequality exactly.
        h_xx = compact.arb(h_xx_text)
        h_xxx = compact.arb(h_xxx_text)
        x_box = bridge.interval_ball(x_low, x_high)
        f_t_raw = -16 * (1 + x_box**4) * h_xx
        f_xt_raw = (
            -64 * x_box**3 * h_xx
            - 16 * (1 + x_box**4) * h_xxx
        )
        f_t_text = serialize_interval(f_t_raw)
        f_xt_text = serialize_interval(f_xt_raw)
        f_t = compact.arb(f_t_text)
        f_xt = compact.arb(f_xt_text)
        f_t_upper = f_t.abs_upper().upper()
        f_xt_upper = f_xt.abs_upper().upper()
        transport_norm_upper = (
            f_t_upper**2 + f_xt_upper**2
        ).sqrt().upper()
        transport_upper = (
            fraction_to_arb(TIME_STEP) * transport_norm_upper
        ).upper()

        reference_f = compact.arb(reference["full_f"])
        reference_f_prime = compact.arb(
            reference["full_f_prime"]
        )
        reference_distance = (
            reference_f.abs_lower() ** 2
            + reference_f_prime.abs_lower() ** 2
        ).sqrt().lower()
        if reference_distance <= 0:
            raise RuntimeError("reference phase cell meets the origin")
        ratio = (transport_upper / reference_distance).upper()
        certified = bool(ratio < 1)
        return {
            "status": "certified" if certified else "failed",
            "x_low": str(x_low),
            "x_high": str(x_high),
            "reference_phase_cell": {
                "full_f": reference["full_f"],
                "full_f_prime": reference["full_f_prime"],
                "distance_lower": serialize_interval(
                    reference_distance
                ),
            },
            "retained": {
                "h_xx": serialize_interval(h_xx_retained),
                "h_xxx": serialize_interval(h_xxx_retained),
            },
            "direct_tail": {
                "E2": serialize_interval(self.tail[2]),
                "E3": serialize_interval(self.tail[3]),
            },
            "full": {
                "h_xx": h_xx_text,
                "h_xxx": h_xxx_text,
            },
            "transport": {
                "f_t": f_t_text,
                "f_xt": f_xt_text,
                "f_t_abs_upper": serialize_interval(f_t_upper),
                "f_xt_abs_upper": serialize_interval(f_xt_upper),
                "jet_time_derivative_norm_upper": serialize_interval(
                    transport_norm_upper
                ),
                "time_step": str(TIME_STEP),
                "displacement_upper": serialize_interval(
                    transport_upper
                ),
                "displacement_to_cell_distance_ratio_upper": (
                    serialize_interval(ratio)
                ),
            },
            "model": {
                "x_center": str(model["x_center"]),
                "x_radius": str(model["x_radius"]),
                "time_center": str(TIME_CENTER),
                "time_radius": str(TIME_RADIUS),
                "oscillatory_moments": (
                    MAX_OSCILLATORY_MOMENT + 1
                ),
                "minimum_integral_accuracy_bits": (
                    model["minimum_integral_accuracy_bits"]
                ),
                "h_xx_remainder": {
                    key: serialize_interval(value)
                    for key, value in model["remainders"][2].items()
                },
                "h_xxx_remainder": {
                    key: serialize_interval(value)
                    for key, value in model["remainders"][3].items()
                },
            },
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
        raise RuntimeError("cache has too many rows")
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


def summarize(
    records: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    ratios = [
        compact.arb(
            record["result"]["transport"][
                "displacement_to_cell_distance_ratio_upper"
            ]
        )
        for record in records
    ]
    complete = (
        len(records) == len(tasks)
        and all(
            record["result"]["status"] == "certified"
            for record in records
        )
    )
    worst_index = (
        max(range(len(ratios)), key=ratios.__getitem__)
        if ratios
        else None
    )
    return {
        "tasks_total": len(tasks),
        "tasks_completed": len(records),
        "certified_panels": sum(
            record["result"]["status"] == "certified"
            for record in records
        ),
        "failed_panels": sum(
            record["result"]["status"] != "certified"
            for record in records
        ),
        "complete_collar_certificate": complete,
        "maximum_transport_ratio_upper": (
            serialize_interval(ratios[worst_index])
            if worst_index is not None
            else None
        ),
        "maximum_transport_ratio_x": (
            [
                records[worst_index]["task"]["x_low"],
                records[worst_index]["task"]["x_high"],
            ]
            if worst_index is not None
            else None
        ),
        "stop_reason": stop_reason,
        "resource": resource,
    }


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    return "\n".join(
        [
            "# Newman Q207-Q208 Adiabatic Bottom-Collar Certificate",
            "",
            f"Date: {DATE}",
            "",
            (
                "Status: rigorous complete bottom-collar certificate."
                if summary["complete_collar_certificate"]
                else "Status: resumable partial bottom-collar computation."
            ),
            "This is a finite calibration, not an all-j successor theorem",
            "and not a proof of `Lambda<=0` or RH.",
            "",
            "## Domain",
            "",
            "```text",
            "time collar: [1/1040,1/1035]",
            "x range: [0,245]",
            "panels: half-unit",
            "reference: Q208 bottom phase cells at t=1/1040",
            "```",
            "",
            "For `F=16*(1+x^4)*H`, the exact heat equation gives",
            "",
            "```text",
            "F_t=-16*(1+x^4)*H_xx,",
            "F_xt=-64*x^3*H_xx-16*(1+x^4)*H_xxx.",
            "```",
            "",
            "Each panel proves",
            "",
            "```text",
            "delta*sup||(F_t,F_xt)|| < dist(0,C_k),",
            "```",
            "",
            "where `C_k` is the stored full Q208 bottom phase cell.",
            "Reverse triangle therefore excludes the origin throughout",
            "the complete collar over that panel.",
            "",
            "## Progress",
            "",
            "```text",
            (
                f"panels={summary['tasks_completed']}/"
                f"{summary['tasks_total']}"
            ),
            f"certified={summary['certified_panels']}",
            f"failed={summary['failed_panels']}",
            (
                "maximum transport ratio upper="
                f"{summary['maximum_transport_ratio_upper']}"
            ),
            (
                "worst panel="
                f"{summary['maximum_transport_ratio_x']}"
            ),
            f"stop reason={summary['stop_reason']}",
            "```",
            "",
            "## Consequence",
            "",
            (
                "When complete, this proves the actual Q207-to-Q208 "
                "bottom collar contact-free and validates the adiabatic "
                "successor interface on one finite transition."
                if summary["complete_collar_certificate"]
                else "No collar theorem is promoted from a partial cache."
            ),
            "",
            "The already certified derivative-positive strip",
            "`[1/1040,1/5]x[245,246]`, together with the Q207 base, then",
            "gives an independent finite no-contact decomposition of the",
            "Q208 low rectangle.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def write_artifact(
    out: Path,
    note: Path,
    cache: Path,
    records: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    summary = summarize(
        records,
        tasks,
        stop_reason,
        resource,
    )
    complete = summary["complete_collar_certificate"]
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "rigorous complete Q207-to-Q208 adiabatic bottom-collar certificate"
            if complete
            else "resumable partial Q207-to-Q208 bottom-collar computation"
        ),
        "proof_boundary": (
            "The complete status certifies only "
            "[1/1040,1/1035]x[0,245] by transport from the stored "
            "Q208 bottom phase cells. It calibrates one finite successor. "
            "No Q208-to-Q209 transport bound, no uniform all-j collar "
            "estimate, no all-j right-strip cone, no cofinal theorem, "
            "no Lambda<=0, and no RH proof is supplied."
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
    out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")
    return artifact


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
        help="Maximum new panels; -1 means unlimited.",
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
    references = load_reference_leaves()
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
            "deferred Q207-Q208 collar certificate: "
            f"baseline mean {baseline_mean:.2f}% exceeds "
            f"{BASELINE_LIMIT:.2f}%"
        )
        return 0

    started = time.monotonic()
    certifier = CollarPanelCertifier()
    previous_hash = (
        records[-1]["record_sha256"] if records else "GENESIS"
    )
    newly_completed = 0
    high_cpu_streak = 0
    stop_reason = "complete"
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
        reference = references[task["sequence"] - 1]
        result = certifier.certify_panel(task, reference)
        record = cache_record(task, result, previous_hash)
        append_cache(args.cache, record)
        records.append(record)
        previous_hash = record["record_sha256"]
        newly_completed += 1
        cpu_sample = float(psutil.cpu_percent(interval=0.25))
        resource["runtime_cpu_samples"].append(cpu_sample)
        ratio = result["transport"][
            "displacement_to_cell_distance_ratio_upper"
        ]
        if (
            task["sequence"] % 10 == 0
            or task["sequence"] == len(tasks)
            or result["status"] != "certified"
        ):
            print(
                f"Q207-Q208 collar {task['sequence']}/{len(tasks)} "
                f"x={task['x_low']}..{task['x_high']} "
                f"status={result['status']} ratio={ratio} "
                f"in {result['elapsed_seconds']:.3f}s "
                f"cpu={cpu_sample:.1f}%",
                flush=True,
            )
        if result["status"] != "certified":
            stop_reason = "failed_panel"
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
        "Q207-Q208 collar artifact: "
        f"{summary['tasks_completed']}/{summary['tasks_total']} panels, "
        f"certified={summary['certified_panels']}, "
        f"failed={summary['failed_panels']}, "
        f"complete={summary['complete_collar_certificate']}, "
        f"stop={summary['stop_reason']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
