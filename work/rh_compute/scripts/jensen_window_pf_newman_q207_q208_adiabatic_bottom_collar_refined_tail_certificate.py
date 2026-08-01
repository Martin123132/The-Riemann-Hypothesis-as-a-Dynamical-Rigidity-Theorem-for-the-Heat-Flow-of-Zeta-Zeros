#!/usr/bin/env python3
"""Refine the Q207-Q208 collar tail after the coarse certified prefix."""

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

import jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_certificate as coarse


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_refined_tail_certificate"
)
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
COARSE_CACHE = RESULT_DIR / f"{coarse.STEM}.jsonl"
COARSE_RESULT = RESULT_DIR / f"{coarse.STEM}.json"

SCHEMA = "newman_q207_q208_adiabatic_refined_tail_v1"
DATE = "2026-07-25"
PREFIX_END = Fraction(379, 2)
X_UPPER = Fraction(245)
REFINED_X_STEP = Fraction(1, 4)
REFINED_TIME_TAYLOR_ORDER = 6
REFINED_X_TAYLOR_ORDER = coarse.X_TAYLOR_ORDER
MAX_OSCILLATORY_MOMENT = (
    max(coarse.DERIVATIVE_ORDERS)
    + 2 * REFINED_TIME_TAYLOR_ORDER
    + REFINED_X_TAYLOR_ORDER
)
BASELINE_SECONDS = 6
BASELINE_LIMIT = 60.0
CPU_PARK_THRESHOLD = 75.0
CPU_CONSECUTIVE_PARK_SAMPLES = 2
DEFAULT_RUNTIME_LIMIT_SECONDS = 10_800.0
compact = coarse.compact


class RefinedCollarPanelCertifier(coarse.CollarPanelCertifier):
    """Use a higher time order on the cancellation-sensitive tail."""

    def remainder(
        self,
        derivative_order: int,
        x_radius: Fraction,
    ) -> dict:
        time_radius = coarse.fraction_to_arb(coarse.TIME_RADIUS)
        x_radius_ball = coarse.fraction_to_arb(x_radius)
        time_moment = (
            derivative_order
            + 2 * (REFINED_TIME_TAYLOR_ORDER + 1)
        )
        time_remainder = (
            compact.arb_power(
                time_radius,
                REFINED_TIME_TAYLOR_ORDER + 1,
            )
            / math.factorial(REFINED_TIME_TAYLOR_ORDER + 1)
            * self.base.absolute_moment(
                time_moment,
                coarse.TIME_HIGH,
            ).upper()
        )
        x_remainder = compact.arb(0)
        time_power = compact.arb(1)
        for time_order in range(REFINED_TIME_TAYLOR_ORDER + 1):
            moment = (
                derivative_order
                + 2 * time_order
                + REFINED_X_TAYLOR_ORDER
                + 1
            )
            x_remainder += (
                time_power
                / math.factorial(time_order)
                * compact.arb_power(
                    x_radius_ball,
                    REFINED_X_TAYLOR_ORDER + 1,
                )
                / math.factorial(REFINED_X_TAYLOR_ORDER + 1)
                * self.base.absolute_moment(
                    moment,
                    coarse.TIME_CENTER,
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
        for derivative_order in coarse.DERIVATIVE_ORDERS:
            rows: list[list] = []
            for time_order in range(
                REFINED_TIME_TAYLOR_ORDER + 1
            ):
                row: list = []
                for x_order in range(
                    REFINED_X_TAYLOR_ORDER + 1
                ):
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
            for derivative_order in coarse.DERIVATIVE_ORDERS
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
        tau = coarse.bridge.interval_ball(
            coarse.TIME_LOW - coarse.TIME_CENTER,
            coarse.TIME_HIGH - coarse.TIME_CENTER,
        )
        xi = coarse.bridge.interval_ball(
            x_low - model["x_center"],
            x_high - model["x_center"],
        )
        rows = model["coefficients"][derivative_order]
        value = compact.arb(0)
        for time_order in range(
            REFINED_TIME_TAYLOR_ORDER,
            -1,
            -1,
        ):
            inner = compact.arb(0)
            for x_order in range(
                REFINED_X_TAYLOR_ORDER,
                -1,
                -1,
            ):
                inner = inner * xi + rows[time_order][x_order]
            value = value * tau + inner
        return value + compact.arb(
            0,
            model["remainders"][derivative_order][
                "total"
            ].str(coarse.ENDPOINT_DIGITS, more=True),
        )


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def refined_tasks() -> list[dict]:
    tasks: list[dict] = []
    current = PREFIX_END
    sequence = 1
    while current < X_UPPER:
        following = min(current + REFINED_X_STEP, X_UPPER)
        tasks.append(
            {
                "sequence": sequence,
                "x_low": str(current),
                "x_high": str(following),
            }
        )
        current = following
        sequence += 1
    return tasks


def load_coarse_prefix() -> list[dict]:
    artifact = json.loads(COARSE_RESULT.read_text(encoding="utf-8"))
    records = artifact["records"]
    prefix = records[:379]
    if len(prefix) != 379:
        raise RuntimeError("coarse prefix is incomplete")
    expected = coarse.panel_tasks()[:379]
    previous_hash = "GENESIS"
    for record, task in zip(prefix, expected, strict=True):
        coarse.validate_cache_record(record, task, previous_hash)
        previous_hash = record["record_sha256"]
        if record["result"]["status"] != "certified":
            raise RuntimeError("coarse prefix contains a failed panel")
    if Fraction(prefix[-1]["task"]["x_high"]) != PREFIX_END:
        raise RuntimeError("coarse prefix endpoint mismatch")
    return prefix


def parent_references() -> list[dict]:
    references = coarse.load_reference_leaves()
    parents: list[dict] = []
    for task in refined_tasks():
        x_low = Fraction(task["x_low"])
        parent_index = int(x_low / coarse.X_STEP)
        reference = references[parent_index]
        if not (
            Fraction(reference["x_low"]) <= x_low
            and Fraction(task["x_high"])
            <= Fraction(reference["x_high"])
        ):
            raise RuntimeError("refined task escaped its reference cell")
        parents.append(reference)
    return parents


def source_hashes() -> dict:
    return {
        "coarse_builder": file_hash(Path(coarse.__file__).resolve()),
        "coarse_cache": file_hash(COARSE_CACHE),
        "coarse_result": file_hash(COARSE_RESULT),
        "q208_reference_cache": file_hash(coarse.REFERENCE_CACHE),
        "q208_reference_result": file_hash(coarse.REFERENCE_RESULT),
        "successor_lemma": file_hash(coarse.SUCCESSOR_RESULT),
    }


def contract_payload() -> dict:
    return {
        "schema": SCHEMA,
        "coarse_certified_prefix": ["0", str(PREFIX_END)],
        "refined_x_domain": [str(PREFIX_END), str(X_UPPER)],
        "refined_x_step": str(REFINED_X_STEP),
        "time_collar": [str(coarse.TIME_LOW), str(coarse.TIME_HIGH)],
        "time_step": str(coarse.TIME_STEP),
        "time_taylor_order": REFINED_TIME_TAYLOR_ORDER,
        "x_taylor_order": REFINED_X_TAYLOR_ORDER,
        "derivative_orders": list(coarse.DERIVATIVE_ORDERS),
        "jet_scale": "ell=1",
        "strict_gate": "delta*sup||(F_t,F_xt)||<dist(0,C_parent)",
        "composition": (
            "379 certified half-unit prefix panels plus 222 "
            "quarter-unit refined-tail panels cover [0,245]"
        ),
        "source_sha256": source_hashes(),
    }


def contract_hash() -> str:
    return sha256(canonical_json(contract_payload()).encode()).hexdigest()


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
        raise RuntimeError("refined cache schema mismatch")
    if record.get("contract_sha256") != contract_hash():
        raise RuntimeError("refined cache contract mismatch")
    if record.get("sequence") != task["sequence"]:
        raise RuntimeError("refined cache sequence mismatch")
    if record.get("task") != task:
        raise RuntimeError("refined cache task mismatch")
    if record.get("previous_hash") != previous_hash:
        raise RuntimeError("refined cache chain mismatch")
    claimed = record.get("record_sha256")
    replay = dict(record)
    replay.pop("record_sha256", None)
    if sha256(canonical_json(replay).encode()).hexdigest() != claimed:
        raise RuntimeError("refined cache record hash mismatch")


def load_cache(path: Path, tasks: list[dict]) -> list[dict]:
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) > len(tasks):
        raise RuntimeError("refined cache has too many rows")
    records: list[dict] = []
    previous_hash = "GENESIS"
    for index, line in enumerate(lines):
        if not line.strip():
            raise RuntimeError("refined cache contains a blank row")
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


def ratio_ball(record: dict):
    return compact.arb(
        record["result"]["transport"][
            "displacement_to_cell_distance_ratio_upper"
        ]
    )


def summarize(
    prefix: list[dict],
    refined: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    all_records = prefix + refined
    ratios = [ratio_ball(record) for record in all_records]
    worst_index = max(
        range(len(ratios)),
        key=ratios.__getitem__,
    )
    refined_failures = sum(
        record["result"]["status"] != "certified"
        for record in refined
    )
    complete = (
        len(prefix) == 379
        and len(refined) == len(tasks)
        and refined_failures == 0
    )
    worst_record = all_records[worst_index]
    return {
        "coarse_prefix_panels": len(prefix),
        "refined_tasks_total": len(tasks),
        "refined_tasks_completed": len(refined),
        "refined_certified_panels": (
            len(refined) - refined_failures
        ),
        "refined_failed_panels": refined_failures,
        "combined_cover": ["0", str(X_UPPER)],
        "combined_panel_count": len(all_records),
        "complete_collar_certificate": complete,
        "maximum_transport_ratio_upper": coarse.serialize_interval(
            ratios[worst_index]
        ),
        "maximum_transport_ratio_x": [
            worst_record["task"]["x_low"],
            worst_record["task"]["x_high"],
        ],
        "maximum_transport_ratio_source": (
            "coarse_prefix"
            if worst_index < len(prefix)
            else "refined_tail"
        ),
        "stop_reason": stop_reason,
        "resource": resource,
    }


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    return "\n".join(
        [
            "# Newman Q207-Q208 Refined Adiabatic Bottom Collar",
            "",
            f"Date: {DATE}",
            "",
            (
                "Status: rigorous complete hybrid collar certificate."
                if summary["complete_collar_certificate"]
                else "Status: resumable refined-tail computation."
            ),
            "This is one finite successor calibration, not an all-j",
            "theorem and not a proof of `Lambda<=0` or RH.",
            "",
            "## Cover",
            "",
            "```text",
            "coarse prefix: 379 half-unit panels on [0,189.5]",
            "refined tail: 222 quarter-unit panels on [189.5,245]",
            "time collar: [1/1040,1/1035]",
            "jet scale: ell=1",
            "```",
            "",
            "Each refined derivative model is centered on its quarter",
            "panel. Its transport bound is compared with the distance of",
            "the enclosing parent Q208 bottom phase cell from the origin.",
            "",
            "## Progress",
            "",
            "```text",
            (
                f"refined={summary['refined_tasks_completed']}/"
                f"{summary['refined_tasks_total']}"
            ),
            (
                f"refined certified="
                f"{summary['refined_certified_panels']}"
            ),
            f"refined failed={summary['refined_failed_panels']}",
            (
                "combined panels="
                f"{summary['combined_panel_count']}"
            ),
            (
                "maximum ratio upper="
                f"{summary['maximum_transport_ratio_upper']}"
            ),
            f"worst panel={summary['maximum_transport_ratio_x']}",
            (
                "worst source="
                f"{summary['maximum_transport_ratio_source']}"
            ),
            f"stop reason={summary['stop_reason']}",
            "```",
            "",
            "## Consequence",
            "",
            (
                "The combined cover proves the complete actual bottom "
                "collar [1/1040,1/1035]x[0,245] contact-free."
                if summary["complete_collar_certificate"]
                else "No complete collar theorem is promoted yet."
            ),
            "",
            "Together with Q207 and the independently certified",
            "derivative-positive strip `[1/1040,1/5]x[245,246]`, a",
            "complete collar result gives a second finite construction",
            "of the Q208 low rectangle.",
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
    prefix: list[dict],
    refined: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    summary = summarize(
        prefix,
        refined,
        tasks,
        stop_reason,
        resource,
    )
    complete = summary["complete_collar_certificate"]
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "rigorous complete hybrid Q207-to-Q208 adiabatic bottom-collar certificate"
            if complete
            else "resumable refined-tail Q207-to-Q208 collar computation"
        ),
        "proof_boundary": (
            "The complete status combines a certified coarse prefix and "
            "quarter-unit refined tail only for "
            "[1/1040,1/1035]x[0,245]. It validates one finite "
            "successor. No Q208-to-Q209 transport budget, no uniform "
            "all-j collar estimate, no all-j right-strip cone, no "
            "cofinal theorem, no Lambda<=0, and no RH proof is supplied."
        ),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "builder_sha256": file_hash(Path(__file__).resolve()),
        "cache": str(cache.relative_to(REPO_ROOT)).replace("\\", "/"),
        "cache_sha256": file_hash(cache) if cache.exists() else None,
        "cache_last_record_sha256": (
            refined[-1]["record_sha256"] if refined else "GENESIS"
        ),
        "summary": summary,
        "refined_records": refined,
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
        help="Maximum new refined panels; -1 means unlimited.",
    )
    parser.add_argument(
        "--runtime-limit-seconds",
        type=float,
        default=DEFAULT_RUNTIME_LIMIT_SECONDS,
    )
    args = parser.parse_args()

    args.cache.parent.mkdir(parents=True, exist_ok=True)
    args.cache.touch(exist_ok=True)
    prefix = load_coarse_prefix()
    tasks = refined_tasks()
    references = parent_references()
    records = load_cache(args.cache, tasks)
    priority_lowered = coarse.compact.request_below_normal_priority()
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
            prefix,
            records,
            tasks,
            "deferred_baseline_cpu",
            resource,
        )
        print(
            "deferred refined collar tail: "
            f"baseline mean {baseline_mean:.2f}% exceeds "
            f"{BASELINE_LIMIT:.2f}%"
        )
        return 0

    started = time.monotonic()
    certifier = RefinedCollarPanelCertifier()
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
        result = certifier.certify_panel(
            task,
            references[task["sequence"] - 1],
        )
        record = cache_record(task, result, previous_hash)
        append_cache(args.cache, record)
        records.append(record)
        previous_hash = record["record_sha256"]
        newly_completed += 1
        cpu_sample = float(psutil.cpu_percent(interval=0.25))
        resource["runtime_cpu_samples"].append(cpu_sample)
        if (
            task["sequence"] % 10 == 0
            or task["sequence"] == len(tasks)
            or result["status"] != "certified"
        ):
            print(
                f"refined collar {task['sequence']}/{len(tasks)} "
                f"x={task['x_low']}..{task['x_high']} "
                f"status={result['status']} "
                "ratio="
                f"{result['transport']['displacement_to_cell_distance_ratio_upper']} "
                f"in {result['elapsed_seconds']:.3f}s "
                f"cpu={cpu_sample:.1f}%",
                flush=True,
            )
        if result["status"] != "certified":
            stop_reason = "failed_refined_panel"
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
        prefix,
        records,
        tasks,
        stop_reason,
        resource,
    )
    summary = artifact["summary"]
    print(
        "refined collar artifact: "
        f"{summary['refined_tasks_completed']}/"
        f"{summary['refined_tasks_total']} refined panels, "
        f"failed={summary['refined_failed_panels']}, "
        f"complete={summary['complete_collar_certificate']}, "
        f"stop={summary['stop_reason']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
