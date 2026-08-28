#!/usr/bin/env python3
"""Transport one grouped-752 panel atlas to arbitrary event-cell boxes.

The original grouped transport scout anchors its panel integrals at the center
of the target height box.  That is ideal for a single theorem box but repeats
the expensive panel integrations when an adaptive cover needs many boxes.
Here the panel anchor and target boxes are independent.  The exact identity

    F_t(y) = F_c(y) exp(-i (t-c) log y)

lets one anchor atlas serve every target box in the same fixed-roster event
cell.  This remains a diagnostic scout until a finite cover and an independent
configuration have been promoted by a separate theorem gate.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_grouped_752_panel_phase_transport_scout as transport


complete = transport.complete
joined = transport.joined
upper_group = transport.upper_group
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_anchored_grouped_752_phase_transport_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
HEIGHT_TEXT = "10000000000"
DEFAULT_ANCHOR_OFFSET = "12.75"
COMPONENT_PRECISION_BITS = 384
VERSION = "anchored-grouped-752-common-phase-transport-v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def center_from_offset(offset_text: str) -> str:
    return format(Decimal(HEIGHT_TEXT) + Decimal(offset_text), "f")


def parse_target(specification: str) -> tuple[str, str, str]:
    try:
        offset_text, radius_text = specification.split(":", maxsplit=1)
        offset = Decimal(offset_text)
        radius = Decimal(radius_text)
    except Exception as error:
        raise RuntimeError(
            f"target must have OFFSET:RADIUS form, got {specification!r}"
        ) from error
    require(radius > 0, f"target radius must be positive: {specification}")
    return format(offset, "f"), center_from_offset(format(offset, "f")), format(radius, "f")


def evaluate_anchored_model(
    panel_rows: list[dict[str, Any]],
    component_rows: dict[str, dict[str, Any]],
    *,
    anchor_text: str,
    target_text: str,
    target_offset: str,
    radius_text: str,
    precision_bits: int,
    variant: str,
) -> dict[str, Any]:
    lower_key = complete.cache_key(
        target_text, radius_text, precision_bits, variant, "lower_ordinary"
    )
    require(lower_key in component_rows, f"missing lower-ordinary component row: {lower_key}")
    lower_row = component_rows[lower_key]
    require(lower_row.get("passed") is True, "lower-ordinary component row failed")
    record = lower_row["record"]

    t_box = arb(arb(target_text), arb(radius_text))
    anchor = arb(anchor_text)
    h = t_box - anchor
    h_abs = abs(h).upper()
    global_log_center = (upper_group.LOWER.log() + upper_group.SPLIT.log()) / 2
    imaginary = acb(0, 1)
    phase_sum = acb(0)
    log_phase_sum = acb(0)
    base_error = arb(0)
    log_error = arb(0)
    row_hashes: list[str] = []
    for row in panel_rows:
        log_center = transport.real_from(row["log_center"])
        log_half_width = transport.real_from(row["log_half_width"]).upper()
        l1_majorant = transport.real_from(row["l1_majorant"]).upper()
        phase = (-imaginary * h * (log_center - global_log_center)).exp()
        phase_sum += phase * transport.complex_from(row["base_panel_integral"])
        log_phase_sum += phase * transport.complex_from(row["log_panel_integral"])
        base_error += h_abs * log_half_width * l1_majorant
        log_error += (
            upper_group.SPLIT.log().upper()
            * h_abs
            * log_half_width
            * l1_majorant
        )
        row_hashes.append(transport.canonical_hash(row))

    hardy_shift = transport.complex_from(record["phase_and_prefactor"]["Hardy_shift"])
    transported = hardy_shift * phase_sum - imaginary * log_phase_sum
    transport_error = (
        abs(hardy_shift).upper() * base_error.upper() + log_error.upper()
    ).upper()
    transported = complete.lower_ordinary.ibp.add_complex_error(
        transported, transport_error
    )
    common_phase = (-imaginary * h * global_log_center).exp()
    grouped_K = common_phase * transported

    theta = transport.real_from(record["phase_and_prefactor"]["theta"])
    theta_prime = transport.real_from(record["phase_and_prefactor"]["theta_prime"])
    theta_anchor, _ = joined.theta_and_derivative(anchor)
    H = transport.real_from(record["phase_and_prefactor"]["H"])
    require(H.lower() > 0, "H lost positivity")
    stable_theta = theta_anchor + h * theta_prime
    stable_effective_theta = theta_anchor + h * (theta_prime - global_log_center)
    direct_effective_theta = theta - h * global_log_center
    require(stable_theta.overlaps(theta), "anchor theta transport misses direct theta")
    require(
        stable_effective_theta.overlaps(direct_effective_theta),
        "anchor effective phase misses direct subtraction",
    )
    recentered_projection = joined.hardy_projection(
        stable_effective_theta, transported
    ) / H
    direct_projection = joined.hardy_projection(theta, grouped_K) / H
    require(
        recentered_projection.overlaps(direct_projection),
        "anchor phase-factor projections miss",
    )

    old_grouped_K = transport.complex_from(
        record["grouped_752_K"]["weighted_grouped_752_ball"]
    )
    require(grouped_K.overlaps(old_grouped_K), "anchored grouped K misses direct enclosure")

    lower_finite = transport.complex_from(record["lower_finite_cell_K"])
    lower_complement = transport.complex_from(
        record["lower_complement_K"]["weighted_tail_ball"]
    )
    finite_upper = transport.complex_from(
        record["finite_upper_K"]["weighted_tail_ball"]
    )
    remote_upper = transport.complex_from(
        record["remote_upper_K"]["weighted_tail_ball"]
    )
    gamma_radius = transport.real_from(
        record["Gamma_defect_K"]["Hardy_operator_component_radius"]
    ).upper()
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
    non_grouped_projection = joined.hardy_projection(
        stable_theta, non_grouped_lower
    ) / H
    replacement_lower_projection = non_grouped_projection - recentered_projection
    require(
        replacement_lower_projection.overlaps(replacement_lower_projection_direct),
        "stable anchored lower split misses direct projection",
    )
    old_lower_projection = transport.real_from(
        record["joined_lower_ordinary_Q_derivative_contribution"]
    )
    require(
        replacement_lower_projection.overlaps(old_lower_projection),
        "anchored lower projection misses old enclosure",
    )

    selected = {
        component: component_rows[
            complete.cache_key(
                target_text, radius_text, precision_bits, variant, component
            )
        ]
        for component in complete.COMPONENTS
    }
    require(all(row.get("passed") is True for row in selected.values()), "a component row failed")
    transition_record = selected["transition"]["record"]
    arc_record = selected["upper_arc"]["record"]
    tail_record = selected["positive_tail"]["record"]
    correction_record = selected["tiny_correction"]["record"]
    transition_K = transport.complex_from(transition_record["transition_K_direct"])
    upper_arc_K = transport.complex_from(arc_record["actual_K_direct"])
    correction_K = transport.complex_from(
        correction_record["compact_weighted_mode_sum_K_corr"]
    )
    tail_modulus = transport.real_from(tail_record["modulus_bound"]).upper()
    tail_K = acb(arb(0, tail_modulus), arb(0, tail_modulus))
    replacement_complete_K = (
        replacement_lower + transition_K + upper_arc_K + correction_K + tail_K
    )
    replacement_direct = joined.hardy_projection(theta, replacement_complete_K) / H
    transition_projection = transport.real_from(
        transition_record["transition_Q_derivative_contribution"]
    )
    upper_arc_projection = transport.real_from(
        arc_record["upper_arc_Q_derivative_contribution"]
    )
    correction_projection = transport.real_from(
        correction_record["Hardy_projection_over_H"]
    )
    tail_projection = transport.real_from(
        tail_record["Hardy_projection_over_H_absolute_bound"]
    ).upper()
    replacement_components = (
        replacement_lower_projection
        + transition_projection
        + upper_arc_projection
        + correction_projection
        + arb(0, tail_projection)
    )
    require(
        replacement_direct.overlaps(replacement_components),
        "anchored complete projections miss",
    )
    old_assembly = complete.assemble(
        component_rows, target_text, radius_text, precision_bits, variant
    )
    require(
        replacement_direct.overlaps(
            transport.real_from(old_assembly["Q_K_minus_T_derivative_direct"])
        ),
        "anchored direct projection misses old complete enclosure",
    )
    require(
        replacement_components.overlaps(
            transport.real_from(
                old_assembly["Q_K_minus_T_derivative_component_sum"]
            )
        ),
        "anchored component projection misses old complete enclosure",
    )

    return {
        "anchor_height": anchor_text,
        "target_offset": target_offset,
        "target_center_height": target_text,
        "radius": radius_text,
        "closed_target_interval": [
            t_box.lower().str(55, more=True),
            t_box.upper().str(55, more=True),
        ],
        "variant": variant,
        "global_log_center": transport.real_record(global_log_center),
        "theta_anchor": transport.real_record(theta_anchor),
        "stable_theta": transport.real_record(stable_theta),
        "stable_effective_theta": transport.real_record(stable_effective_theta),
        "direct_effective_theta": transport.real_record(direct_effective_theta),
        "base_phase_replacement_error": transport.real_record(base_error),
        "log_phase_replacement_error": transport.real_record(log_error),
        "grouped_K_transport_error": transport.real_record(transport_error),
        "grouped_752_K_recentered": transport.complex_record(grouped_K),
        "grouped_752_K_direct_overlap": True,
        "grouped_752_Q_derivative_recentered": transport.real_record(
            recentered_projection
        ),
        "replacement_lower_ordinary_K": transport.complex_record(replacement_lower),
        "replacement_lower_ordinary_Q_derivative": transport.real_record(
            replacement_lower_projection
        ),
        "replacement_complete_K_T": transport.complex_record(
            replacement_complete_K
        ),
        "replacement_Q_K_minus_T_derivative_direct": transport.real_record(
            replacement_direct
        ),
        "replacement_Q_K_minus_T_derivative_component_sum": transport.real_record(
            replacement_components
        ),
        "strictly_positive_replacement_derivative": bool(
            replacement_components.lower() > 0
        ),
        "strictly_negative_replacement_derivative": bool(
            replacement_components.upper() < 0
        ),
        "absolute_derivative_upper": abs(replacement_components).upper().str(
            55, more=True
        ),
        "replacement_direct_component_overlap": True,
        "replacement_old_overlaps": True,
        "panel_row_hashes": row_hashes,
        "component_row_hashes": {
            component: complete.canonical_hash(row)
            for component, row in selected.items()
        },
    }


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchor-offset", default=DEFAULT_ANCHOR_OFFSET)
    parser.add_argument("--target", action="append", required=True)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    parser.add_argument("--panels", type=int)
    parser.add_argument("--panel-precision-bits", type=int)
    parser.add_argument("--max-new-panels", type=int)
    args = parser.parse_args()

    panels = args.panels or (192 if args.variant == "production" else 256)
    panel_precision_bits = args.panel_precision_bits or (
        384 if args.variant == "production" else 448
    )
    require(panels >= 32, "at least 32 panels required")
    require(panel_precision_bits >= 256, "panel precision must be at least 256 bits")
    if args.max_new_panels is not None:
        require(args.max_new_panels >= 0, "max-new-panels must be nonnegative")

    targets = [parse_target(specification) for specification in args.target]
    anchor_text = center_from_offset(args.anchor_offset)
    ctx.prec = panel_precision_bits
    ctx.threads = 1
    resource_mode = transport.base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")

    panel_cache = transport.load_jsonl(transport.CACHE)
    new_count = 0
    order = list(range(panels))
    if args.variant == "independent":
        order.reverse()
    for index in order:
        key = transport.panel_key(
            anchor_text, panels, panel_precision_bits, args.variant, index
        )
        if key in panel_cache:
            continue
        if args.max_new_panels is not None and new_count >= args.max_new_panels:
            break
        row = transport.integrate_panel(
            anchor_text,
            arb(anchor_text),
            panels,
            panel_precision_bits,
            args.variant,
            index,
        )
        transport.append_jsonl(transport.CACHE, row)
        panel_cache[key] = row
        new_count += 1
        print(
            f"anchor_panel={index + 1}/{panels} elapsed={row['elapsed_seconds']:.3f}s",
            flush=True,
        )

    complete_panel_set = all(
        transport.panel_key(
            anchor_text, panels, panel_precision_bits, args.variant, index
        )
        in panel_cache
        for index in range(panels)
    )
    models: list[dict[str, Any]] = []
    if complete_panel_set:
        selected_panels = transport.selected_panels(
            panel_cache,
            anchor_text,
            panels,
            panel_precision_bits,
            args.variant,
        )
        component_rows = transport.load_jsonl(complete.CACHE)
        for target_offset, target_text, radius_text in targets:
            model = evaluate_anchored_model(
                selected_panels,
                component_rows,
                anchor_text=anchor_text,
                target_text=target_text,
                target_offset=target_offset,
                radius_text=radius_text,
                precision_bits=COMPONENT_PRECISION_BITS,
                variant=args.variant,
            )
            models.append(model)
            derivative = model["replacement_Q_K_minus_T_derivative_component_sum"]
            print(
                f"target={target_offset}:{radius_text} "
                f"lower={derivative['lower']} upper={derivative['upper']} "
                f"positive={model['strictly_positive_replacement_derivative']} "
                f"negative={model['strictly_negative_replacement_derivative']}",
                flush=True,
            )

    source_paths = (
        Path(__file__).resolve(),
        Path(transport.__file__).resolve(),
        Path(complete.__file__).resolve(),
        Path(upper_group.__file__).resolve(),
        Path(transport.base.__file__).resolve(),
    )
    artifact = {
        "kind": "rh_diagnostic_anchored_grouped_752_common_phase_transport_scout",
        "date": "2026-08-28",
        "status": (
            "diagnostic_anchored_grouped_752_transport_complete"
            if complete_panel_set
            else "diagnostic_anchored_grouped_752_panel_cache_incomplete"
        ),
        "passed": complete_panel_set and len(models) == len(targets),
        "version": VERSION,
        "resource_policy": "one below-normal process, one logical CPU, one numerical thread",
        "anchor_offset": args.anchor_offset,
        "anchor_height": anchor_text,
        "targets": [
            {"offset": offset, "center_height": center, "radius": radius}
            for offset, center, radius in targets
        ],
        "variant": args.variant,
        "panels": panels,
        "panel_precision_bits": panel_precision_bits,
        "component_precision_bits": COMPONENT_PRECISION_BITS,
        "new_panel_count": new_count,
        "completed_panel_count": sum(
            transport.panel_key(
                anchor_text, panels, panel_precision_bits, args.variant, index
            )
            in panel_cache
            for index in range(panels)
        ),
        "panel_cache_path": transport.CACHE.relative_to(ROOT).as_posix(),
        "component_cache_path": complete.CACHE.relative_to(ROOT).as_posix(),
        "models": models,
        "sources": {
            path.relative_to(ROOT).as_posix(): sha256(path) for path in source_paths
        },
        "diagnostic_boundary": (
            "This scout reuses one grouped-752 anchor atlas on specified target boxes. "
            "It does not itself prove that the targets form a cover, certify a value bound, "
            "cross an event wall, or imply RH."
        ),
    }
    atomic_write_json(RESULT, artifact)
    if not complete_panel_set:
        print(
            f"checkpointed {artifact['completed_panel_count']}/{panels} anchor panels",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
