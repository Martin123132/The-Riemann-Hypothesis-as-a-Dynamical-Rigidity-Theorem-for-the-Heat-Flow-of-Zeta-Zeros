#!/usr/bin/env python3
"""Recenter the grouped 752-label Hardy packet panel by panel.

The direct interval integrator forgets that all 192 panels share nearly the
same height phase.  This scout stores rigorous center-height panel integrals,
factors one common logarithmic rotation, and bounds only the within-panel
logarithmic variation.  It is a diagnostic replacement for one component of
the complete K_T packet, not a theorem gate.
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

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate as upper_group
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_local_complete_K_T_scout as complete


base = upper_group.base
joined = complete.joined
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_grouped_752_panel_phase_transport_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
CACHE = ROOT / "work" / "rh_compute" / "results" / "cache" / f"{STEM}_panels.jsonl"
COMPONENT_CACHE = complete.CACHE
HEIGHT_TEXT = "10000000000"
DEFAULT_OFFSET = "20"
DEFAULT_RADIUS = "0.006"
VERSION = "grouped-752-common-phase-panel-transport-v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def center_from_offset(offset_text: str) -> str:
    return format(Decimal(HEIGHT_TEXT) + Decimal(offset_text), "f")


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return complete.complex_record(value, digits)


def real_record(value: arb, digits: int = 70) -> dict[str, str]:
    return complete.real_record(value, digits)


def complex_from(record: dict[str, Any]) -> acb:
    return complete.complex_from(record)


def real_from(record: dict[str, Any]) -> arb:
    return complete.real_from(record)


def panel_key(center_text: str, panels: int, precision_bits: int, variant: str, index: int) -> str:
    return f"{center_text}|{panels}|{precision_bits}|{variant}|{index}|{VERSION}"


def canonical_hash(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_jsonl(path: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    if not path.is_file():
        return rows
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            key = row.get("cache_key")
            require(key, f"cache row {line_number} has no key: {path}")
            require(key not in rows, f"duplicate cache key at row {line_number}: {path}")
            rows[key] = row
    return rows


def append_jsonl(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def integrate_panel(
    center_text: str,
    center: arb,
    panels: int,
    precision_bits: int,
    variant: str,
    index: int,
) -> dict[str, Any]:
    started = time.perf_counter()
    width = upper_group.SPLIT - upper_group.LOWER
    left = width * index / panels
    right = width * (index + 1) / panels
    y_left = upper_group.LOWER + left
    y_right = upper_group.LOWER + right
    log_left = y_left.log()
    log_right = y_right.log()
    log_center = (log_left + log_right) / 2
    log_half_width = ((log_right - log_left) / 2).upper()
    l1_majorant = (
        arb(upper_group.PACKET_COUNT) * (right - left) / y_left.sqrt()
    ).upper()
    base_integrand = base.grouped_integrand_factory(
        t=center,
        lower=upper_group.LOWER,
        q0=upper_group.Q_PLUS,
        count=upper_group.PACKET_COUNT,
    )

    def log_integrand(x: acb, analytic: bool) -> acb:
        y = acb(upper_group.LOWER) + x
        return y.log() * base_integrand(x, analytic)

    base_value = acb.integral(
        base_integrand,
        left,
        right,
        abs_tol=base.TOLERANCE,
        rel_tol=base.TOLERANCE,
        eval_limit=300_000,
        depth_limit=40,
    )
    log_value = acb.integral(
        log_integrand,
        left,
        right,
        abs_tol=base.TOLERANCE,
        rel_tol=base.TOLERANCE,
        eval_limit=300_000,
        depth_limit=40,
    )
    require(base_value.is_finite() and log_value.is_finite(), f"nonfinite panel {index}")
    return {
        "cache_key": panel_key(center_text, panels, precision_bits, variant, index),
        "version": VERSION,
        "center_height": center_text,
        "panels": panels,
        "precision_bits": precision_bits,
        "variant": variant,
        "index": index,
        "left_offset": left.str(45, more=True),
        "right_offset": right.str(45, more=True),
        "log_center": real_record(log_center),
        "log_half_width": real_record(log_half_width),
        "l1_majorant": real_record(l1_majorant),
        "base_panel_integral": complex_record(base_value),
        "log_panel_integral": complex_record(log_value),
        "elapsed_seconds": time.perf_counter() - started,
        "passed": True,
    }


def selected_panels(
    rows: dict[str, dict[str, Any]],
    center_text: str,
    panels: int,
    precision_bits: int,
    variant: str,
) -> list[dict[str, Any]]:
    selected = []
    for index in range(panels):
        key = panel_key(center_text, panels, precision_bits, variant, index)
        require(key in rows, f"missing panel row: {key}")
        require(rows[key].get("passed") is True, f"failed panel row: {key}")
        selected.append(rows[key])
    return selected


def evaluate_model(
    panel_rows: list[dict[str, Any]],
    component_rows: dict[str, dict[str, Any]],
    center_text: str,
    radius_text: str,
    precision_bits: int,
    variant: str,
) -> dict[str, Any]:
    lower_key = complete.cache_key(center_text, radius_text, precision_bits, variant, "lower_ordinary")
    require(lower_key in component_rows, f"missing lower-ordinary component row: {lower_key}")
    lower_row = component_rows[lower_key]
    require(lower_row.get("passed") is True, "lower-ordinary component row failed")
    record = lower_row["record"]

    t_box = arb(arb(center_text), arb(radius_text))
    h = t_box - arb(center_text)
    h_abs = abs(h).upper()
    global_log_center = (
        upper_group.LOWER.log() + upper_group.SPLIT.log()
    ) / 2
    imaginary = acb(0, 1)
    phase_sum = acb(0)
    log_phase_sum = acb(0)
    base_error = arb(0)
    log_error = arb(0)
    row_hashes: list[str] = []
    for row in panel_rows:
        log_center = real_from(row["log_center"])
        log_half_width = real_from(row["log_half_width"]).upper()
        l1_majorant = real_from(row["l1_majorant"]).upper()
        phase = (-imaginary * h * (log_center - global_log_center)).exp()
        phase_sum += phase * complex_from(row["base_panel_integral"])
        log_phase_sum += phase * complex_from(row["log_panel_integral"])
        base_error += h_abs * log_half_width * l1_majorant
        log_error += (
            upper_group.SPLIT.log().upper()
            * h_abs
            * log_half_width
            * l1_majorant
        )
        row_hashes.append(canonical_hash(row))

    hardy_shift = complex_from(record["phase_and_prefactor"]["Hardy_shift"])
    transported = hardy_shift * phase_sum - imaginary * log_phase_sum
    transport_error = (
        abs(hardy_shift).upper() * base_error.upper() + log_error.upper()
    ).upper()
    transported = complete.lower_ordinary.ibp.add_complex_error(transported, transport_error)
    common_phase = (-imaginary * h * global_log_center).exp()
    grouped_K = common_phase * transported

    theta = real_from(record["phase_and_prefactor"]["theta"])
    theta_prime = real_from(record["phase_and_prefactor"]["theta_prime"])
    theta_center, _ = joined.theta_and_derivative(arb(center_text))
    H = real_from(record["phase_and_prefactor"]["H"])
    require(H.lower() > 0, "H lost positivity")
    stable_theta = theta_center + h * theta_prime
    stable_effective_theta = theta_center + h * (theta_prime - global_log_center)
    direct_effective_theta = theta - h * global_log_center
    require(stable_theta.overlaps(theta), "stable theta transport misses direct theta")
    require(
        stable_effective_theta.overlaps(direct_effective_theta),
        "stable effective phase misses direct subtraction",
    )
    recentered_projection = (
        joined.hardy_projection(stable_effective_theta, transported) / H
    )
    direct_projection = joined.hardy_projection(theta, grouped_K) / H
    require(recentered_projection.overlaps(direct_projection), "phase-factor projections miss")

    old_grouped_K = complex_from(record["grouped_752_K"]["weighted_grouped_752_ball"])
    require(grouped_K.overlaps(old_grouped_K), "recentered grouped K misses direct enclosure")

    lower_finite = complex_from(record["lower_finite_cell_K"])
    lower_complement = complex_from(record["lower_complement_K"]["weighted_tail_ball"])
    finite_upper = complex_from(record["finite_upper_K"]["weighted_tail_ball"])
    remote_upper = complex_from(record["remote_upper_K"]["weighted_tail_ball"])
    gamma_radius = real_from(record["Gamma_defect_K"]["Hardy_operator_component_radius"]).upper()
    gamma_rectangle = acb(arb(0, gamma_radius), arb(0, gamma_radius))
    non_grouped_lower = (
        lower_finite
        - lower_complement
        - finite_upper
        - remote_upper
        + gamma_rectangle
    )
    replacement_lower = non_grouped_lower - grouped_K
    replacement_lower_projection_direct = (
        joined.hardy_projection(theta, replacement_lower) / H
    )
    non_grouped_projection = joined.hardy_projection(stable_theta, non_grouped_lower) / H
    replacement_lower_projection = non_grouped_projection - recentered_projection
    require(
        replacement_lower_projection.overlaps(replacement_lower_projection_direct),
        "stable split lower projection misses direct projection",
    )
    old_lower_projection = real_from(record["joined_lower_ordinary_Q_derivative_contribution"])
    require(replacement_lower_projection.overlaps(old_lower_projection), "replacement lower projection misses old enclosure")

    selected = {
        component: component_rows[
            complete.cache_key(center_text, radius_text, precision_bits, variant, component)
        ]
        for component in complete.COMPONENTS
    }
    require(all(row.get("passed") is True for row in selected.values()), "a complete component row failed")
    transition_record = selected["transition"]["record"]
    arc_record = selected["upper_arc"]["record"]
    tail_record = selected["positive_tail"]["record"]
    correction_record = selected["tiny_correction"]["record"]
    transition_K = complex_from(transition_record["transition_K_direct"])
    upper_arc_K = complex_from(arc_record["actual_K_direct"])
    correction_K = complex_from(correction_record["compact_weighted_mode_sum_K_corr"])
    tail_modulus = real_from(tail_record["modulus_bound"]).upper()
    tail_K = acb(arb(0, tail_modulus), arb(0, tail_modulus))
    replacement_complete_K = (
        replacement_lower + transition_K + upper_arc_K + correction_K + tail_K
    )
    replacement_direct = joined.hardy_projection(theta, replacement_complete_K) / H
    transition_projection = real_from(transition_record["transition_Q_derivative_contribution"])
    upper_arc_projection = real_from(arc_record["upper_arc_Q_derivative_contribution"])
    correction_projection = real_from(correction_record["Hardy_projection_over_H"])
    tail_projection = real_from(tail_record["Hardy_projection_over_H_absolute_bound"]).upper()
    replacement_components = (
        replacement_lower_projection
        + transition_projection
        + upper_arc_projection
        + correction_projection
        + arb(0, tail_projection)
    )
    require(replacement_direct.overlaps(replacement_components), "replacement complete projections miss")
    old_assembly = complete.assemble(
        component_rows, center_text, radius_text, precision_bits, variant
    )
    require(
        replacement_direct.overlaps(real_from(old_assembly["Q_K_minus_T_derivative_direct"])),
        "replacement direct projection misses old complete enclosure",
    )
    require(
        replacement_components.overlaps(
            real_from(old_assembly["Q_K_minus_T_derivative_component_sum"])
        ),
        "replacement component projection misses old complete enclosure",
    )

    return {
        "center_height": center_text,
        "radius": radius_text,
        "height_ball": t_box.str(60, more=True),
        "panels": len(panel_rows),
        "precision_bits": precision_bits,
        "variant": variant,
        "global_log_center": real_record(global_log_center),
        "theta_center": real_record(theta_center),
        "stable_theta": real_record(stable_theta),
        "stable_effective_theta": real_record(stable_effective_theta),
        "direct_effective_theta": real_record(direct_effective_theta),
        "base_phase_replacement_error": real_record(base_error),
        "log_phase_replacement_error": real_record(log_error),
        "grouped_K_transport_error": real_record(transport_error),
        "grouped_752_K_recentered": complex_record(grouped_K),
        "grouped_752_K_direct_overlap": True,
        "grouped_752_Q_derivative_recentered": real_record(recentered_projection),
        "replacement_lower_ordinary_K": complex_record(replacement_lower),
        "replacement_lower_ordinary_Q_derivative_direct": real_record(
            replacement_lower_projection_direct
        ),
        "replacement_lower_ordinary_Q_derivative": real_record(replacement_lower_projection),
        "replacement_complete_K_T": complex_record(replacement_complete_K),
        "replacement_Q_K_minus_T_derivative_direct": real_record(replacement_direct),
        "replacement_Q_K_minus_T_derivative_component_sum": real_record(replacement_components),
        "replacement_direct_component_overlap": True,
        "old_complete_derivative_direct": old_assembly["Q_K_minus_T_derivative_direct"],
        "old_complete_derivative_component_sum": old_assembly[
            "Q_K_minus_T_derivative_component_sum"
        ],
        "replacement_old_overlaps": True,
        "strictly_positive_replacement_derivative": bool(
            replacement_components.lower() > 0
        ),
        "strictly_positive_replacement_direct_projection": bool(
            replacement_direct.lower() > 0
        ),
        "panel_row_hashes": row_hashes,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offset", default=DEFAULT_OFFSET)
    parser.add_argument("--radius", default=DEFAULT_RADIUS)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    parser.add_argument("--panels", type=int)
    parser.add_argument("--precision-bits", type=int)
    parser.add_argument("--max-new-panels", type=int)
    args = parser.parse_args()

    panels = args.panels or (192 if args.variant == "production" else 256)
    precision_bits = args.precision_bits or (384 if args.variant == "production" else 448)
    require(panels >= 32, "at least 32 panels required")
    require(precision_bits >= 256, "precision must be at least 256 bits")
    if args.max_new_panels is not None:
        require(args.max_new_panels >= 0, "max-new-panels must be nonnegative")

    ctx.prec = precision_bits
    ctx.threads = 1
    resource_mode = base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    center_text = center_from_offset(args.offset)
    center = arb(center_text)
    cached = load_jsonl(CACHE)
    order = list(range(panels))
    if args.variant == "independent":
        order.reverse()
    new_count = 0
    for index in order:
        key = panel_key(center_text, panels, precision_bits, args.variant, index)
        if key in cached:
            continue
        if args.max_new_panels is not None and new_count >= args.max_new_panels:
            break
        row = integrate_panel(
            center_text, center, panels, precision_bits, args.variant, index
        )
        append_jsonl(CACHE, row)
        cached[key] = row
        new_count += 1
        print(
            f"panel={index + 1}/{panels} elapsed={row['elapsed_seconds']:.3f}s",
            flush=True,
        )

    complete_panel_set = all(
        panel_key(center_text, panels, precision_bits, args.variant, index) in cached
        for index in range(panels)
    )
    model = None
    if complete_panel_set:
        panels_selected = selected_panels(
            cached, center_text, panels, precision_bits, args.variant
        )
        component_rows = load_jsonl(COMPONENT_CACHE)
        model = evaluate_model(
            panels_selected,
            component_rows,
            center_text,
            args.radius,
            384,
            args.variant,
        )

    source_paths = (
        Path(__file__).resolve(),
        Path(upper_group.__file__).resolve(),
        Path(base.__file__).resolve(),
        Path(complete.__file__).resolve(),
    )
    artifact = {
        "kind": "rh_diagnostic_grouped_752_common_phase_panel_transport_scout",
        "date": "2026-08-28",
        "status": (
            "diagnostic_grouped_752_panel_transport_complete"
            if complete_panel_set
            else "diagnostic_grouped_752_panel_cache_incomplete"
        ),
        "passed": complete_panel_set,
        "version": VERSION,
        "resource_policy": "one below-normal process, one logical CPU, one numerical thread",
        "offset": args.offset,
        "center_height": center_text,
        "radius": args.radius,
        "panels": panels,
        "panel_precision_bits": precision_bits,
        "component_precision_bits": 384,
        "variant": args.variant,
        "new_panel_count": new_count,
        "panel_cache_path": CACHE.relative_to(ROOT).as_posix(),
        "component_cache_path": COMPONENT_CACHE.relative_to(ROOT).as_posix(),
        "completed_panel_count": sum(
            panel_key(center_text, panels, precision_bits, args.variant, index) in cached
            for index in range(panels)
        ),
        "model": model,
        "sources": {path.relative_to(ROOT).as_posix(): sha256(path) for path in source_paths},
        "diagnostic_boundary": (
            "This scout replaces only the grouped 752-label collar by a common-phase panel "
            "transport. It is not independently promoted, does not cover a height gap or event "
            "wall, and does not imply RH."
        ),
    }
    atomic_write_json(RESULT, artifact)
    if model is not None:
        print(
            "assembled panel-phase replacement: "
            f"lower={model['replacement_Q_K_minus_T_derivative_component_sum']['lower']} "
            f"positive={model['strictly_positive_replacement_derivative']}",
            flush=True,
        )
    else:
        print(f"checkpointed {artifact['completed_panel_count']}/{panels} panels", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
