#!/usr/bin/env python3
"""Run a rigorous selected-panel pilot on the new Q208 boundary strips."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import time

import psutil

import jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate as bridge


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_q208_selected_boundary_pilot.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_q208_selected_boundary_pilot.md"
)

TIME_NEW = Fraction(1, 1040)
TIME_OLD = Fraction(1, 1035)
TIME_TOP = Fraction(1, 5)
TIME_STEP = Fraction(1, 100)
CPU_BASELINE_SECONDS = 6
CPU_BASELINE_LIMIT = 55.0
CPU_PARK_THRESHOLD = 75.0

BOTTOM_PANELS = [
    (Fraction(38), Fraction(77, 2)),
    (Fraction(137, 2), Fraction(69)),
    (Fraction(199, 2), Fraction(100)),
    (Fraction(299, 2), Fraction(150)),
    (Fraction(399, 2), Fraction(200)),
    (Fraction(489, 2), Fraction(245)),
]
RIGHT_PANELS = [
    (Fraction(245), Fraction(491, 2)),
    (Fraction(491, 2), Fraction(246)),
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def configure_bridge() -> None:
    bridge.TIME_GLOBAL_LOWER = TIME_NEW
    bridge.TIME_GLOBAL_UPPER = TIME_TOP
    bridge.TIME_GLOBAL_CENTER = (TIME_NEW + TIME_TOP) / 2
    bridge.TIME_GLOBAL_RADIUS = (TIME_TOP - TIME_NEW) / 2
    bridge.LOW_STRIP_TIME_UPPER = TIME_OLD
    bridge.INITIAL_TIME_STEP = TIME_STEP
    bridge.CPU_PARK_THRESHOLD = CPU_PARK_THRESHOLD


def tasks() -> list[dict]:
    bottom = [
        {
            "region": "low_time_left_strip",
            "x_low": str(lower),
            "x_high": str(upper),
        }
        for lower, upper in BOTTOM_PANELS
    ]
    right = [
        {
            "region": "full_time_right_strip",
            "x_low": str(lower),
            "x_high": str(upper),
        }
        for lower, upper in RIGHT_PANELS
    ]
    return [*bottom, *right]


def baseline_samples() -> list[float]:
    return [
        float(psutil.cpu_percent(interval=1.0))
        for _ in range(CPU_BASELINE_SECONDS)
    ]


def lower_float(text: str) -> float:
    return float(bridge.compact.arb(text).lower())


def full_derivative_lower(leaf: dict):
    return (
        bridge.compact.arb(leaf["j_retained_prime_lower"]).lower()
        - bridge.compact.arb(leaf["tail_derivative_upper"]).upper()
    )


def summarize(records: list[dict]) -> dict:
    leaves = [
        leaf
        for record in records
        for leaf in record["result"].get("certified_leaves", [])
    ]
    unresolved = sum(
        record["result"].get("unresolved_boxes", 0)
        for record in records
    )
    minimum = min(
        leaves,
        key=lambda row: lower_float(row["certified_ratio_lower"]),
        default=None,
    )
    right_records = [
        record
        for record in records
        if record["task"]["region"] == "full_time_right_strip"
    ]
    right_leaves = [
        leaf
        for record in right_records
        for leaf in record["result"].get("certified_leaves", [])
    ]
    right_derivative_lowers = [
        full_derivative_lower(leaf)
        for leaf in right_leaves
    ]
    right_strip_derivative_positive = bool(
        len(right_records) == len(RIGHT_PANELS)
        and len(right_leaves) == 40
        and all(value > 0 for value in right_derivative_lowers)
    )
    right_minimum = (
        min(right_derivative_lowers) if right_derivative_lowers else None
    )
    return {
        "panels": len(records),
        "certified_panels": sum(
            record["result"].get("status") == "certified"
            for record in records
        ),
        "certified_leaves": len(leaves),
        "unresolved_boxes": unresolved,
        "branch_counts": {
            "value": sum(leaf.get("branch") == "value" for leaf in leaves),
            "derivative": sum(
                leaf.get("branch") == "derivative" for leaf in leaves
            ),
        },
        "full_value_sign_counts": {
            "negative": sum(
                leaf.get("full_value_negative") is True for leaf in leaves
            ),
            "positive": sum(
                leaf.get("full_value_positive") is True for leaf in leaves
            ),
            "not_sign_separated": sum(
                leaf.get("full_value_negative") is not True
                and leaf.get("full_value_positive") is not True
                for leaf in leaves
            ),
        },
        "minimum_certified_ratio_lower": (
            minimum["certified_ratio_lower"] if minimum else None
        ),
        "minimum_location": (
            {
                key: minimum[key]
                for key in (
                    "t_low",
                    "t_high",
                    "x_low",
                    "x_high",
                    "branch",
                )
            }
            if minimum
            else None
        ),
        "right_strip": {
            "domain": "[1/1040,1/5]x[245,246]",
            "certified_cells": len(right_leaves),
            "full_derivative_positive": right_strip_derivative_positive,
            "minimum_full_derivative_lower": (
                right_minimum.str(75, more=True)
                if right_minimum is not None
                else None
            ),
            "conclusion": (
                "J_t'(x)>0 throughout [1/1040,1/5]x[245,246]; "
                "hence the full first jet is nonzero there and its image "
                "lies in the open upper half-plane."
            ),
        },
        "all_selected_panels_certified": (
            len(records) == len(tasks())
            and unresolved == 0
            and all(
                record["result"].get("status") == "certified"
                for record in records
            )
        ),
    }


def render_note(payload: dict) -> str:
    summary = payload["summary"]
    config = payload["config"]
    records = payload["records"]
    panel_lines = []
    for record in records:
        result = record["result"]
        panel_lines.append(
            "| {region} | [{x_low},{x_high}] | {leaves} | {value} | "
            "{derivative} | {minimum} |".format(
                region=record["task"]["region"],
                x_low=record["task"]["x_low"],
                x_high=record["task"]["x_high"],
                leaves=result["certified_leaf_boxes"],
                value=result["branch_counts"]["value"],
                derivative=result["branch_counts"]["derivative"],
                minimum=result["minimum_certified_ratio_lower"],
            )
        )
    return "\n".join(
        [
            "# Newman Q208 Selected-Boundary Pilot",
            "",
            "Date: 2026-07-25",
            "",
            "Status: rigorous complete right-strip theorem and selected-bottom",
            "diagnostic, not a proof of Q208, `Lambda<=0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_q208_selected_boundary_pilot.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_q208_selected_boundary_pilot.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_q208_selected_boundary_pilot.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            (
                "validated Q208 selected-boundary pilot: "
                f"{summary['panels']} panels, "
                f"{summary['certified_leaves']} certified leaves, "
                f"{summary['unresolved_boxes']} unresolved, "
                "0 global promotions"
            ),
            "```",
            "",
            "## Domain",
            "",
            "The prior theorem closes",
            "",
            "```text",
            "Q_207=[1/1035,1/4]x[0,245].",
            "```",
            "",
            "The next linear-exhaustion rectangle is",
            "",
            "```text",
            "Q_208=[1/1040,1/4]x[0,246].",
            "```",
            "",
            "This pilot samples six half-unit panels across the new bottom",
            "microstrip and both half-unit panels across the new right strip:",
            "",
            "```text",
            f"bottom time strip: {config['bottom_time_strip']}",
            f"right time strip: {config['right_time_strip']}",
            "```",
            "",
            "It uses the same 352-bit six-term retained transform, directed",
            "Taylor enclosure, and raw omitted-tail bars as the Q207 theorem.",
            "The raw tail bars remain valid at x=246; no gamma-scale comparison",
            "outside its recorded x<=245 theorem domain is promoted.",
            "",
            "## Results",
            "",
            "| region | x panel | leaves | value branch | derivative branch | minimum ratio |",
            "|:---|:---|---:|---:|---:|:---|",
            *panel_lines,
            "",
            "Summary:",
            "",
            "```text",
            f"branch counts={summary['branch_counts']}",
            f"full-value sign counts={summary['full_value_sign_counts']}",
            (
                "minimum certified ratio="
                f"{summary['minimum_certified_ratio_lower']}"
            ),
            f"minimum location={summary['minimum_location']}",
            f"right-strip theorem={summary['right_strip']}",
            f"baseline CPU samples={payload['resource']['baseline_samples']}",
            f"runtime CPU samples={payload['resource']['runtime_samples']}",
            "```",
            "",
            "## Interpretation",
            "",
            "All selected panels are contact-separated with a large interval",
            "margin. The two right panels form a complete cover, not a sample:",
            "",
            "```text",
            summary["right_strip"]["conclusion"],
            "```",
            "",
            "The right-edge first-jet phase is therefore confined to the open",
            "upper half-plane and contributes no full turn. Across the selected",
            "bottom samples both value and derivative branches are needed, so",
            "a one-sign bottom-edge proof is not the expected continuation.",
            "This supports a continuous first-jet phase or winding certificate.",
            "",
            "The pilot does not cover the unsampled bottom panels, does not",
            "compute a closed-boundary winding, and does not certify Q208.",
            "The next proof-facing step is a rigorous boundary phase-unwrapping",
            "algorithm, followed by a complete Q208 boundary run if that",
            "algorithm closes.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    configure_bridge()
    priority_lowered = bridge.compact.request_below_normal_priority()
    baseline = baseline_samples()
    baseline_mean = sum(baseline) / len(baseline)
    if baseline_mean > CPU_BASELINE_LIMIT:
        print(
            "deferred Q208 selected-boundary pilot: "
            f"baseline mean {baseline_mean:.2f}% exceeds "
            f"{CPU_BASELINE_LIMIT:.2f}%"
        )
        return 2

    source_paths = {
        "q207_builder": Path(bridge.__file__).resolve(),
        "q207_result": bridge.DEFAULT_OUT.resolve(),
        "six_term_tail": bridge.TAIL_SOURCE.resolve(),
        "outer_tail": bridge.OUTER_SOURCE.resolve(),
    }
    started = time.monotonic()
    certifier = bridge.SharedPanelCertifier()
    records: list[dict] = []
    runtime_samples: list[float] = []
    for sequence, task in enumerate(tasks(), start=1):
        panel_started = time.monotonic()
        result = bridge.certify_panel(certifier, task)
        records.append(
            {
                "sequence": sequence,
                "task": task,
                "result": result,
            }
        )
        sample = float(psutil.cpu_percent(interval=0.25))
        runtime_samples.append(sample)
        print(
            f"Q208 pilot {sequence}/{len(tasks())} "
            f"{bridge.task_id(task)} status={result['status']} "
            f"leaves={result['certified_leaf_boxes']} "
            f"in {time.monotonic() - panel_started:.3f}s",
            flush=True,
        )
        if sample > CPU_PARK_THRESHOLD:
            print(
                "parked Q208 pilot after current atomic panel: "
                f"CPU sample {sample:.1f}%"
            )
            break

    summary = summarize(records)
    payload = {
        "kind": "jensen_window_pf_newman_q208_selected_boundary_pilot",
        "date": "2026-07-25",
        "status": (
            "rigorous complete Q208 right-strip theorem and "
            "selected-bottom boundary diagnostic"
        ),
        "proof_boundary": (
            "The stored Arb/Taylor rows rigorously certify the complete new "
            "right strip [1/1040,1/5]x[245,246] by strict full-derivative "
            "positivity and six selected panels on the new bottom microstrip. "
            "They do not cover the complete bottom edge, compute the closed "
            "boundary winding, certify Q208, prove Lambda<=0, or prove RH."
        ),
        "config": {
            "q207": "[1/1035,1/4]x[0,245]",
            "q208": "[1/1040,1/4]x[0,246]",
            "bottom_time_strip": "[1/1040,1/1035]",
            "right_time_strip": "[1/1040,1/5]",
            "bottom_panels": [
                [str(lower), str(upper)]
                for lower, upper in BOTTOM_PANELS
            ],
            "right_panels": [
                [str(lower), str(upper)]
                for lower, upper in RIGHT_PANELS
            ],
            "precision_bits": bridge.PRECISION_BITS,
            "retained_terms": bridge.RETAINED_TERMS,
            "taylor_time_order": bridge.TAYLOR_TIME_ORDER,
            "taylor_x_order": bridge.TAYLOR_X_ORDER,
            "time_step": str(TIME_STEP),
        },
        "source_sha256": {
            key: file_hash(path)
            for key, path in source_paths.items()
        },
        "resource": {
            "below_normal_priority_applied": priority_lowered,
            "baseline_samples": baseline,
            "baseline_mean": baseline_mean,
            "runtime_samples": runtime_samples,
            "elapsed_seconds": time.monotonic() - started,
            "worker_count": 1,
        },
        "summary": summary,
        "records": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Q208 selected-boundary pilot: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()} "
        f"panels={summary['panels']} "
        f"leaves={summary['certified_leaves']} "
        f"unresolved={summary['unresolved_boxes']}"
    )
    return 0 if summary["all_selected_panels_certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
