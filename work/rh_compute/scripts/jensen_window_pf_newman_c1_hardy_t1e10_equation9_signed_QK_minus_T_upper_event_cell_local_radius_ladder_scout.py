#!/usr/bin/env python3
"""Measure the direct Arb sign-box radius barrier inside the upper event cell.

This is a resumable diagnostic, not a theorem gate.  It evaluates increasing
height-box radii at one or more fixed centers, keeps all analytic donor guards,
and stops each ladder at the first failed guard or loss of strict negativity.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_point_scout as point


base = point.base
DEFAULT_OFFSETS = ("20",)
DEFAULT_RADII = (
    "0.000001",
    "0.00001",
    "0.000025",
    "0.00005",
    "0.000075",
    "0.0001",
    "0.00011",
    "0.00012",
    "0.00013",
    "0.00014",
    "0.00016",
    "0.00018",
    "0.0002",
    "0.0005",
    "0.001",
)
VERSION = "local-midpoint-action-radius-ladder-v1"
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_local_radius_ladder_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
CACHE = ROOT / "work" / "rh_compute" / "results" / "cache" / f"{STEM}_rows.jsonl"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cache_key(
    center_text: str,
    radius_text: str,
    precision_bits: int,
    variant: str,
) -> str:
    return f"{center_text}|{radius_text}|{precision_bits}|{variant}|{VERSION}"


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
            require(row["cache_key"] not in rows, f"duplicate cache key at row {line_number}")
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


def center_from_offset(offset_text: str) -> str:
    return format(Decimal(point.HEIGHT_TEXT) + Decimal(offset_text), "f")


def install_ladder_hooks() -> tuple[Any, Any]:
    original_require = base.require
    original_arc = base.interval_upper_arc

    def ladder_require(condition: bool, message: str) -> None:
        # Both messages are wrapper-calibration assertions, not inequalities
        # used in any contour, roster, phase, or remainder estimate.
        if message in (
            "scout radius must lie in (0,0.001]",
            "height box misses the saved-height value",
        ):
            return
        original_require(condition, message)

    base.require = ladder_require
    base.interval_upper_arc = point.locally_anchored_interval_upper_arc
    return original_require, original_arc


def restore_ladder_hooks(original_require: Any, original_arc: Any) -> None:
    base.require = original_require
    base.interval_upper_arc = original_arc


def run_row(
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
    require(event_lower.upper() < t_box.lower(), f"box {offset_text}/{radius_text} crosses lower event wall")
    require(t_box.upper() < event_upper.lower(), f"box {offset_text}/{radius_text} crosses upper event wall")
    started = time.perf_counter()
    try:
        record = base.build(
            radius_text,
            precision_bits,
            center_text=center_text,
            variant=variant,
        )
        qkt = record["projection"]["Q_K_minus_T_height_box"]
        return {
            "cache_key": key,
            "version": VERSION,
            "offset": offset_text,
            "center_height": center_text,
            "radius": radius_text,
            "precision_bits": precision_bits,
            "variant": variant,
            "analytic_guards_passed": True,
            "strictly_negative": record["projection"]["strictly_negative_on_box"],
            "Q_K_minus_T": qkt,
            "radius_to_QKT_radius_ratio": record["projection"]["radius_to_QKT_radius_ratio"],
            "complete_packet_width": record["components"]["complete_joined_packet"],
            "joined_transition_width": record["components"]["joined_transition"],
            "upper_arc_action_guard": record["components"]["upper_arc_action_guard"],
            "height_ball": record["height_ball"],
            "elapsed_seconds": time.perf_counter() - started,
        }
    except Exception as error:
        return {
            "cache_key": key,
            "version": VERSION,
            "offset": offset_text,
            "center_height": center_text,
            "radius": radius_text,
            "precision_bits": precision_bits,
            "variant": variant,
            "analytic_guards_passed": False,
            "strictly_negative": False,
            "error_type": type(error).__name__,
            "error": str(error),
            "elapsed_seconds": time.perf_counter() - started,
        }


def summarize(
    rows: list[dict[str, Any]],
    offsets: tuple[str, ...],
    radii: tuple[str, ...],
    precision_bits: int,
    variant: str,
) -> dict[str, Any]:
    ladders: list[dict[str, Any]] = []
    for offset in offsets:
        selected = [row for row in rows if row["offset"] == offset]
        negative = [row for row in selected if row["analytic_guards_passed"] and row["strictly_negative"]]
        losses = [row for row in selected if not (row["analytic_guards_passed"] and row["strictly_negative"])]
        ladders.append(
            {
                "offset": offset,
                "center_height": center_from_offset(offset),
                "evaluated_row_count": len(selected),
                "largest_tested_strictly_negative_radius": negative[-1]["radius"] if negative else None,
                "first_tested_loss_radius": losses[0]["radius"] if losses else None,
                "first_tested_loss_kind": (
                    "analytic_guard_failure"
                    if losses and not losses[0]["analytic_guards_passed"]
                    else "sign_box_contains_nonnegative_values"
                    if losses
                    else None
                ),
                "rows": selected,
            }
        )

    event_lower, event_upper = point.event_cell()
    return {
        "kind": "rh_diagnostic_joined_QK_minus_T_upper_event_cell_local_radius_ladder_scout",
        "date": "2026-08-28",
        "status": "diagnostic_local_radius_barrier_measurement_not_an_event_cell_theorem",
        "version": VERSION,
        "resource_policy": "one below-normal process, one logical CPU, one numerical thread",
        "precision_bits": precision_bits,
        "variant": variant,
        "requested_offsets": list(offsets),
        "requested_radii_in_order": list(radii),
        "event_cell": {
            "lower_open_wall_ball": event_lower.str(60, more=True),
            "upper_open_wall_ball": event_upper.str(60, more=True),
            "every_evaluated_box_strictly_inside": True,
        },
        "cache_path": str(CACHE.relative_to(ROOT)).replace("\\", "/"),
        "cache_sha256": sha256(CACHE) if CACHE.is_file() else None,
        "source_dependencies": {
            "point_scout": {
                "path": str(Path(point.__file__).resolve().relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(Path(point.__file__).resolve()),
            },
            "direct_interval_scout": {
                "path": str(Path(base.__file__).resolve().relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(Path(base.__file__).resolve()),
            },
            "event_atlas": {
                "path": str(point.EVENT_ATLAS.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(point.EVENT_ATLAS),
            },
        },
        "bypassed_wrapper_calibrations": [
            "radius<=0.001 diagnostic policy cap",
            "overlap with the saved t=10^10 point enclosure",
            "legacy action<0.01981 point alarm",
            "legacy upper-arc modulus<0.057665 point alarm",
        ],
        "retained_substantive_guards": [
            "strict containment in the fixed-roster open event cell",
            "stationary-root and fold placement",
            "all lower and upper nonstationary gaps",
            "transition-cell parity and integral length",
            "unique upper-arc stationary maximum and monotonic tail",
            "upper-arc action_at_cut<-80 and phase_deviation<pi/2",
            "upper-arc total_added_error<1e-20",
            "upper-arc positive_tail_log10<-100000000",
            "all finite-jet, endpoint-expansion, and integration-by-parts remainder bounds",
            "finite Arb enclosures for every assembled component",
        ],
        "ladders": ladders,
        "diagnostic_boundary": (
            "A passing row is a rigorous direct enclosure only for that one closed height box. "
            "The ladder does not fill gaps between centers, prove a maximal-event-cell sign theorem, "
            "establish an event-wall handoff, or prove RH. A failed sign box is an enclosure-width "
            "obstruction and is not evidence that Q_K-T is pointwise nonnegative anywhere."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offset", action="append", dest="offsets")
    parser.add_argument("--radius", action="append", dest="radii")
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    parser.add_argument("--max-new", type=int)
    parser.add_argument("--continue-after-loss", action="store_true")
    args = parser.parse_args()

    offsets = tuple(args.offsets) if args.offsets else DEFAULT_OFFSETS
    radii = tuple(args.radii) if args.radii else DEFAULT_RADII
    require(offsets, "at least one offset is required")
    require(radii, "at least one radius is required")
    require(all(arb(radius) > 0 for radius in radii), "all radii must be positive")
    require(
        all(Decimal(radii[index]) < Decimal(radii[index + 1]) for index in range(len(radii) - 1)),
        "radii must be strictly increasing",
    )
    require(args.precision_bits >= 256, "precision must be at least 256 bits")
    if args.max_new is not None:
        require(args.max_new >= 0, "max-new must be nonnegative")

    ctx.prec = args.precision_bits
    ctx.threads = 1
    event_lower, event_upper = point.event_cell()
    cached = load_cache()
    rows: list[dict[str, Any]] = []
    new_count = 0
    budget_exhausted = False
    original_require, original_arc = install_ladder_hooks()
    try:
        for offset in offsets:
            for radius in radii:
                center = center_from_offset(offset)
                key = cache_key(center, radius, args.precision_bits, args.variant)
                if key in cached:
                    row = cached[key]
                elif args.max_new is not None and new_count >= args.max_new:
                    budget_exhausted = True
                    break
                else:
                    row = run_row(
                        offset,
                        radius,
                        args.precision_bits,
                        args.variant,
                        event_lower,
                        event_upper,
                    )
                    append_cache(row)
                    cached[key] = row
                    new_count += 1
                rows.append(row)
                if row["analytic_guards_passed"]:
                    qkt = row["Q_K_minus_T"]
                    print(
                        f"offset={offset} radius={radius} upper={qkt['upper']} "
                        f"negative={row['strictly_negative']}",
                        flush=True,
                    )
                else:
                    print(
                        f"offset={offset} radius={radius} failed={row['error']}",
                        flush=True,
                    )
                if not args.continue_after_loss and not (
                    row["analytic_guards_passed"] and row["strictly_negative"]
                ):
                    break
            if budget_exhausted:
                break
    finally:
        restore_ladder_hooks(original_require, original_arc)

    artifact = summarize(rows, offsets, radii, args.precision_bits, args.variant)
    artifact["new_row_count"] = new_count
    artifact["max_new_budget_exhausted"] = budget_exhausted
    atomic_write_json(RESULT, artifact)
    print(
        f"completed local radius ladder: rows={len(rows)} new={new_count} "
        f"budget_exhausted={budget_exhausted}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
