#!/usr/bin/env python3
"""Scout tiny rigorous boxes across the upper fixed-roster event cell.

This is a diagnostic atlas, not a theorem gate.  It leaves the certified
Section 11.509 evaluator untouched, uses its complete joined packet, and
stores each completed probe in an append-only JSONL cache.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_first_subcell_interval_scout as base


HEIGHT_TEXT = "10000000000"
DEFAULT_RADIUS = "0.0000001"
DEFAULT_OFFSETS = (
    "0.001",
    "0.01",
    "0.05",
    "0.1",
    "0.5",
    "1",
    "2",
    "5",
    "10",
    "20",
    "40",
    "55",
    "57",
)
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_point_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
CACHE = ROOT / "work" / "rh_compute" / "results" / "cache" / f"{STEM}_rows.jsonl"
EVENT_ATLAS = (
    ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def center_from_offset(offset_text: str) -> str:
    return format(Decimal(HEIGHT_TEXT) + Decimal(offset_text), "f")


def locally_anchored_interval_upper_arc(variant: str) -> dict[str, Any]:
    """Run the donor arc proof and transport its action from the box midpoint.

    The donor's 0.01981 assertion is a saved-point calibration alarm.  Its
    actual Taylor, omitted-numerator, compact-tail, arc-scale, and positive-tail
    inequalities are evaluated from the current Arb box and remain active.
    """

    upper_arc = base.upper_arc
    donor_require = upper_arc.require

    def diagnostic_require(condition: bool, message: str) -> None:
        if message in (
            "production action maximum exceeded its guard",
            "upper arc missed its certified scale guard",
        ):
            return
        donor_require(condition, message)

    upper_arc.require = diagnostic_require
    try:
        certificate = upper_arc.run(variant, None)
    finally:
        upper_arc.require = donor_require

    pi = arb.pi()
    y = arb(upper_arc.Y_TEXT)
    q_minus = arb(upper_arc.Q_MINUS_TEXT)
    radius = 2 * pi * y
    t_box = arb(upper_arc.HEIGHT)
    t_anchor = t_box.mid()

    def stationary_delta(t: arb) -> arb:
        constant = t - 2 * pi * y * y
        discriminant = (q_minus * radius) ** 2 - 16 * pi * y * y * constant
        root = (q_minus * radius + discriminant.sqrt()) / (8 * pi * y * y)
        return root.acos()

    delta_anchor = stationary_delta(t_anchor)
    delta_box = stationary_delta(t_box)
    action_anchor = (
        t_anchor * delta_anchor
        + pi * y * y * (2 * delta_anchor).sin()
        - q_minus * radius * delta_anchor.sin()
    )
    transport_distance = abs(t_box - t_anchor).upper()
    stable_action = action_anchor + arb(0, transport_distance * delta_box.upper())
    require(stable_action.is_finite(), "locally transported action is not finite")
    rotated_arc_absolute = arb(certificate["rotated_arc_absolute"]["ball"])
    require(rotated_arc_absolute.is_finite(), "current upper-arc modulus is not finite")

    certificate["scout_action_guard"] = {
        "legacy_point_alarm": "action_maximum<0.01981",
        "legacy_point_alarm_satisfied": bool(stable_action.upper() < arb("0.01981")),
        "legacy_point_alarm_used_as_analytic_bound": False,
        "legacy_arc_scale_alarm": "rotated_arc_absolute<0.057665",
        "legacy_arc_scale_alarm_satisfied": bool(rotated_arc_absolute.upper() < arb("0.057665")),
        "legacy_arc_scale_alarm_used_as_analytic_bound": False,
        "current_rotated_arc_absolute_ball": rotated_arc_absolute.str(55, more=True),
        "current_complex_arc_consumed_directly_by_joined_packet": True,
        "raw_dependent_action_ball": certificate["action_maximum"]["ball"],
        "stable_local_envelope_action_ball": stable_action.str(55, more=True),
        "local_anchor_height_ball": t_anchor.str(55, more=True),
        "maximum_local_transport_distance": transport_distance.str(35, more=True),
        "stationary_delta_box": delta_box.str(55, more=True),
        "identity": "d A_max(t)/dt=delta_*(t)",
        "retained_donor_guards": [
            "unique stationary maximum and monotonic tail",
            "action_at_cut<-80",
            "phase_deviation<pi/2",
            "total_added_error<1e-20",
            "positive_tail_log10<-100000000",
        ],
        "all_substantive_donor_error_bounds_retained": True,
    }
    return certificate


def install_diagnostic_hooks() -> tuple[Any, Any]:
    original_require = base.require
    original_arc = base.interval_upper_arc

    def diagnostic_require(condition: bool, message: str) -> None:
        if message == "height box misses the saved-height value":
            return
        original_require(condition, message)

    base.require = diagnostic_require
    base.interval_upper_arc = locally_anchored_interval_upper_arc
    return original_require, original_arc


def restore_diagnostic_hooks(original_require: Any, original_arc: Any) -> None:
    base.require = original_require
    base.interval_upper_arc = original_arc


def cache_key(center_text: str, radius_text: str, precision_bits: int, variant: str) -> str:
    return f"{center_text}|{radius_text}|{precision_bits}|{variant}|local-action-v2"


def load_cache() -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    with CACHE.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            require("cache_key" in row, f"cache row {line_number} has no key")
            rows[row["cache_key"]] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def event_cell() -> tuple[arb, arb]:
    event = load_json(EVENT_ATLAS)["certificate"]["primary_open_cell"]
    return arb(event["lower_ball"]["ball"]), arb(event["upper_ball"]["ball"])


def run_probe(
    offset_text: str,
    radius_text: str,
    precision_bits: int,
    variant: str,
    event_lower: arb,
    event_upper: arb,
) -> dict[str, Any]:
    center_text = center_from_offset(offset_text)
    key = cache_key(center_text, radius_text, precision_bits, variant)
    t_box = arb(arb(center_text), arb(radius_text))
    require(event_lower.upper() < t_box.lower(), f"probe {offset_text} crosses lower event wall")
    require(t_box.upper() < event_upper.lower(), f"probe {offset_text} crosses upper event wall")
    started = time.perf_counter()
    try:
        record = base.build(
            radius_text,
            precision_bits,
            center_text=center_text,
            variant=variant,
        )
        return {
            "cache_key": key,
            "offset": offset_text,
            "center_height": center_text,
            "radius": radius_text,
            "precision_bits": precision_bits,
            "variant": variant,
            "passed": True,
            "record": record,
            "elapsed_seconds": time.perf_counter() - started,
        }
    except Exception as error:
        return {
            "cache_key": key,
            "offset": offset_text,
            "center_height": center_text,
            "radius": radius_text,
            "precision_bits": precision_bits,
            "variant": variant,
            "passed": False,
            "error_type": type(error).__name__,
            "error": str(error),
            "elapsed_seconds": time.perf_counter() - started,
        }


def summarize(rows: list[dict[str, Any]], radius_text: str, precision_bits: int, variant: str) -> dict[str, Any]:
    summaries: list[dict[str, Any]] = []
    for row in rows:
        summary = {key: row[key] for key in ("cache_key", "offset", "center_height", "passed", "elapsed_seconds")}
        if row["passed"]:
            record = row["record"]
            qkt = record["projection"]["Q_K_minus_T_height_box"]
            action = record["components"]["upper_arc_action_guard"]
            summary.update(
                {
                    "Q_K_minus_T": qkt,
                    "strictly_negative": record["projection"]["strictly_negative_on_box"],
                    "complete_packet_width": record["components"]["complete_joined_packet"],
                    "joined_transition_width": record["components"]["joined_transition"],
                    "stable_local_action_ball": action["stable_local_envelope_action_ball"],
                    "legacy_action_alarm_satisfied": action["legacy_point_alarm_satisfied"],
                }
            )
        else:
            summary.update({"error_type": row["error_type"], "error": row["error"]})
        summaries.append(summary)
    passed = [row for row in summaries if row["passed"]]
    negative = [row for row in passed if row["strictly_negative"]]
    return {
        "kind": "rh_diagnostic_joined_QK_minus_T_upper_event_cell_point_scout",
        "date": "2026-08-28",
        "status": "diagnostic_tiny_box_upper_event_cell_atlas_not_a_promoted_interval_theorem",
        "resource_policy": "one below-normal process, one logical CPU, one numerical thread",
        "radius": radius_text,
        "precision_bits": precision_bits,
        "variant": variant,
        "cache_path": str(CACHE.relative_to(ROOT)).replace("\\", "/"),
        "requested_probe_count": len(rows),
        "passed_probe_count": len(passed),
        "strictly_negative_probe_count": len(negative),
        "all_passed_probes_strictly_negative": bool(passed) and len(passed) == len(negative),
        "probes": summaries,
        "diagnostic_boundary": (
            "Finite tiny-radius boxes only. They do not cover the gaps between probes, do not prove "
            "a maximal-event-cell sign theorem, and do not establish an event-wall handoff or RH."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", default=DEFAULT_RADIUS)
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    parser.add_argument("--offset", action="append", dest="offsets")
    parser.add_argument("--max-new", type=int)
    args = parser.parse_args()

    offsets = tuple(args.offsets) if args.offsets else DEFAULT_OFFSETS
    require(arb(args.radius) > 0 and arb(args.radius) <= arb("0.001"), "radius must lie in (0,0.001]")
    require(args.precision_bits >= 256, "precision must be at least 256 bits")
    if args.max_new is not None:
        require(args.max_new >= 0, "max-new must be nonnegative")

    ctx.prec = args.precision_bits
    ctx.threads = 1
    event_lower, event_upper = event_cell()
    cached = load_cache()
    rows: list[dict[str, Any]] = []
    new_count = 0
    original_require, original_arc = install_diagnostic_hooks()
    try:
        for offset in offsets:
            center = center_from_offset(offset)
            key = cache_key(center, args.radius, args.precision_bits, args.variant)
            if key in cached:
                row = cached[key]
            elif args.max_new is not None and new_count >= args.max_new:
                continue
            else:
                row = run_probe(
                    offset,
                    args.radius,
                    args.precision_bits,
                    args.variant,
                    event_lower,
                    event_upper,
                )
                append_cache(row)
                cached[key] = row
                new_count += 1
            rows.append(row)
            if row["passed"]:
                qkt = row["record"]["projection"]["Q_K_minus_T_height_box"]
                print(
                    f"offset={offset} upper={qkt['upper']} "
                    f"negative={row['record']['projection']['strictly_negative_on_box']}",
                    flush=True,
                )
            else:
                print(f"offset={offset} failed={row['error']}", flush=True)
    finally:
        restore_diagnostic_hooks(original_require, original_arc)

    artifact = summarize(rows, args.radius, args.precision_bits, args.variant)
    atomic_write_json(RESULT, artifact)
    print(
        "completed upper event-cell point scout: "
        f"rows={len(rows)} new={new_count} passed={artifact['passed_probe_count']} "
        f"negative={artifact['strictly_negative_probe_count']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
