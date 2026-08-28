#!/usr/bin/env python3
"""Scout the complete joined height derivative at an arbitrary event-cell center.

Each expensive ownership component is stored as one fsynced JSONL row before
the next component starts.  This is a diagnostic derivative scout, not a
connected-cover or maximal-event-cell theorem.
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
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_ordinary_K_first_subcell_scout as lower_ordinary
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_height_derivative_scout as transition
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_upper_arc_K_first_subcell_scout as upper_arc
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_positive_real_tail_K_first_subcell_gate as positive_tail
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_tiny_target_correction_K_first_subcell_gate as correction
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_recentered_upper_arc_K as recentered_upper_arc


joined = transition.joined
HEIGHT_TEXT = "10000000000"
DEFAULT_OFFSET = "20"
DEFAULT_RADIUS = "0.0000001"
DEFAULT_PRECISION = 384
VERSION = "recentered-complete-KT-components-v1"
UPPER_ARC_VERSION = "recentered-stable-action-upper-arc-v2"
COMPONENTS = (
    "lower_ordinary",
    "transition",
    "upper_arc",
    "positive_tail",
    "tiny_correction",
)
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_local_complete_K_T_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
CACHE = ROOT / "work" / "rh_compute" / "results" / "cache" / f"{STEM}_components.jsonl"
EVENT_ATLAS = transition.EVENT_ATLAS


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def center_from_offset(offset_text: str) -> str:
    return format(Decimal(HEIGHT_TEXT) + Decimal(offset_text), "f")


def cache_key(
    center_text: str,
    radius_text: str,
    precision_bits: int,
    variant: str,
    component: str,
) -> str:
    algorithm_version = UPPER_ARC_VERSION if component == "upper_arc" else VERSION
    return (
        f"{center_text}|{radius_text}|{precision_bits}|{variant}|"
        f"{component}|{algorithm_version}"
    )


def canonical_hash(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_cache() -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    if not CACHE.is_file():
        return rows
    with CACHE.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            key = row.get("cache_key")
            require(key, f"cache row {line_number} has no key")
            require(key not in rows, f"duplicate cache key at row {line_number}")
            rows[key] = row
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


def complex_from(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def real_from(record: dict[str, Any]) -> arb:
    if "lower" not in record or "upper" not in record:
        return arb(record["ball"])
    lower = arb(record["lower"]).lower()
    upper = arb(record["upper"]).upper()
    midpoint = (lower + upper) / 2
    return arb(midpoint, ((upper - lower) / 2).upper())


def real_record(value: arb, digits: int = 70) -> dict[str, str]:
    return {
        "ball": value.str(digits, more=True),
        "lower": value.lower().str(55, more=True),
        "upper": value.upper().str(55, more=True),
        "radius": value.rad().str(40, more=True),
    }


def complex_record(value: acb, digits: int = 70) -> dict[str, str]:
    return {
        "real_ball": value.real.str(digits, more=True),
        "imag_ball": value.imag.str(digits, more=True),
        "absolute_ball": abs(value).str(55, more=True),
    }


def install_center(center_text: str) -> list[tuple[Any, Any]]:
    modules = (
        transition,
        upper_arc,
        lower_ordinary,
        lower_ordinary.lower_k,
        positive_tail.lower_k,
        correction.lower_k,
    )
    originals: list[tuple[Any, Any]] = []
    seen: set[int] = set()
    for module in modules:
        if id(module) in seen:
            continue
        seen.add(id(module))
        originals.append((module, module.HEIGHT))
        module.HEIGHT = center_text
    return originals


def restore_centers(originals: list[tuple[Any, Any]]) -> None:
    for module, value in originals:
        module.HEIGHT = value


def install_radius_hooks() -> list[tuple[Any, Any]]:
    wrappers = (
        (lower_ordinary, "radius must lie in [0,0.001]"),
        (transition, "radius must lie in [0,0.001]"),
        (upper_arc, "radius must lie in (0,0.001]"),
    )
    originals: list[tuple[Any, Any]] = []
    for module, policy_message in wrappers:
        original = module.require

        def extended_require(
            condition: bool,
            message: str,
            *,
            original_require: Any = original,
            skipped_message: str = policy_message,
        ) -> None:
            if message == skipped_message:
                return
            original_require(condition, message)

        originals.append((module, original))
        module.require = extended_require
    return originals


def restore_radius_hooks(originals: list[tuple[Any, Any]]) -> None:
    for module, original in originals:
        module.require = original


def component_builder(
    component: str,
    center_text: str,
    t_box: arb,
    radius_text: str,
    precision_bits: int,
    variant: str,
) -> Callable[[], dict[str, Any]]:
    independent = variant == "independent"
    if component == "lower_ordinary":
        return lambda: lower_ordinary.build(radius_text, precision_bits, variant)
    if component == "transition":
        return lambda: transition.build(radius_text, precision_bits, variant)
    if component == "upper_arc":
        return lambda: recentered_upper_arc.build(
            center_text,
            radius_text,
            precision_bits,
            variant,
        )
    if component == "positive_tail":
        return lambda: positive_tail.tail_certificate(
            t_box,
            use_duplication_h=independent,
        )
    if component == "tiny_correction":
        return lambda: correction.correction_certificate(
            t_box,
            reverse=independent,
            use_duplication_h=independent,
        )
    raise RuntimeError(f"unknown component: {component}")


def run_component(
    component: str,
    center_text: str,
    radius_text: str,
    precision_bits: int,
    variant: str,
    t_box: arb,
) -> dict[str, Any]:
    key = cache_key(center_text, radius_text, precision_bits, variant, component)
    started = time.perf_counter()
    try:
        record = component_builder(
            component,
            center_text,
            t_box,
            radius_text,
            precision_bits,
            variant,
        )()
        return {
            "cache_key": key,
            "version": UPPER_ARC_VERSION if component == "upper_arc" else VERSION,
            "component": component,
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
            "version": UPPER_ARC_VERSION if component == "upper_arc" else VERSION,
            "component": component,
            "center_height": center_text,
            "radius": radius_text,
            "precision_bits": precision_bits,
            "variant": variant,
            "passed": False,
            "error_type": type(error).__name__,
            "error": str(error),
            "elapsed_seconds": time.perf_counter() - started,
        }


def assemble(rows: dict[str, dict[str, Any]], center_text: str, radius_text: str, precision_bits: int, variant: str) -> dict[str, Any]:
    selected = {
        component: rows[cache_key(center_text, radius_text, precision_bits, variant, component)]
        for component in COMPONENTS
    }
    require(all(row["passed"] for row in selected.values()), "a K_T component did not pass")

    lower = selected["lower_ordinary"]["record"]
    trans = selected["transition"]["record"]
    arc = selected["upper_arc"]["record"]
    tail = selected["positive_tail"]["record"]
    corr = selected["tiny_correction"]["record"]

    lower_K = complex_from(lower["joined_lower_ordinary_K"])
    transition_K = complex_from(trans["transition_K_direct"])
    upper_arc_K = complex_from(arc["actual_K_direct"])
    correction_K = complex_from(corr["compact_weighted_mode_sum_K_corr"])
    tail_modulus = real_from(tail["modulus_bound"]).upper()
    tail_K = acb(arb(0, tail_modulus), arb(0, tail_modulus))
    complete_K = lower_K + (transition_K + upper_arc_K) + correction_K + tail_K

    theta = real_from(lower["phase_and_prefactor"]["theta"])
    H = real_from(lower["phase_and_prefactor"]["H"])
    require(H.lower() > 0, "assembled H lost positivity")
    direct_projection = joined.hardy_projection(theta, complete_K) / H

    lower_projection = real_from(lower["joined_lower_ordinary_Q_derivative_contribution"])
    transition_projection = real_from(trans["transition_Q_derivative_contribution"])
    upper_arc_projection = real_from(arc["upper_arc_Q_derivative_contribution"])
    correction_projection = real_from(corr["Hardy_projection_over_H"])
    tail_projection_bound = real_from(tail["Hardy_projection_over_H_absolute_bound"]).upper()
    component_projection = (
        lower_projection
        + transition_projection
        + upper_arc_projection
        + correction_projection
        + arb(0, tail_projection_bound)
    )
    require(
        direct_projection.overlaps(component_projection),
        "direct complete K_T projection misses component projection sum",
    )

    return {
        "center_height": center_text,
        "radius": radius_text,
        "height_ball": arb(arb(center_text), arb(radius_text)).str(60, more=True),
        "variant": variant,
        "complete_K_T": complex_record(complete_K),
        "component_K_T": {
            "lower_ordinary": complex_record(lower_K),
            "transition": complex_record(transition_K),
            "upper_arc": complex_record(upper_arc_K),
            "positive_tail_rectangle": complex_record(tail_K),
            "tiny_correction": complex_record(correction_K),
        },
        "Q_K_minus_T_derivative_direct": real_record(direct_projection),
        "Q_K_minus_T_derivative_component_sum": real_record(component_projection),
        "direct_component_projection_overlap": True,
        "strictly_positive_derivative": bool(
            direct_projection.lower() > 0 and component_projection.lower() > 0
        ),
        "strictly_negative_derivative": bool(
            direct_projection.upper() < 0 and component_projection.upper() < 0
        ),
        "component_row_hashes": {
            component: canonical_hash(row) for component, row in selected.items()
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offset", default=DEFAULT_OFFSET)
    parser.add_argument("--radius", default=DEFAULT_RADIUS)
    parser.add_argument("--precision-bits", type=int, default=DEFAULT_PRECISION)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    parser.add_argument("--max-new-components", type=int)
    args = parser.parse_args()

    require(arb(args.radius) > 0, "radius must be positive")
    require(args.precision_bits >= 256, "precision must be at least 256 bits")
    if args.max_new_components is not None:
        require(args.max_new_components >= 0, "max-new-components must be nonnegative")

    ctx.prec = args.precision_bits
    ctx.threads = 1
    resource_mode = joined.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    center_text = center_from_offset(args.offset)
    t_box = arb(arb(center_text), arb(args.radius))
    event = json.loads(EVENT_ATLAS.read_text(encoding="utf-8"))["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    require(event_lower.upper() < t_box.lower(), "local K_T box crosses lower event wall")
    require(t_box.upper() < event_upper.lower(), "local K_T box crosses upper event wall")

    cached = load_cache()
    new_count = 0
    originals = install_center(center_text)
    radius_hooks = install_radius_hooks()
    try:
        for component in COMPONENTS:
            key = cache_key(center_text, args.radius, args.precision_bits, args.variant, component)
            if key in cached:
                row = cached[key]
            elif args.max_new_components is not None and new_count >= args.max_new_components:
                break
            else:
                row = run_component(
                    component,
                    center_text,
                    args.radius,
                    args.precision_bits,
                    args.variant,
                    t_box,
                )
                append_cache(row)
                cached[key] = row
                new_count += 1
            if row["passed"]:
                print(
                    f"component={component} passed elapsed={row['elapsed_seconds']:.3f}s",
                    flush=True,
                )
            else:
                print(f"component={component} failed={row['error']}", flush=True)
                break
    finally:
        restore_radius_hooks(radius_hooks)
        restore_centers(originals)

    keys = [
        cache_key(center_text, args.radius, args.precision_bits, args.variant, component)
        for component in COMPONENTS
    ]
    complete = all(key in cached and cached[key]["passed"] for key in keys)
    assembly = (
        assemble(cached, center_text, args.radius, args.precision_bits, args.variant)
        if complete
        else None
    )
    source_paths = (
        Path(__file__).resolve(),
        Path(lower_ordinary.__file__).resolve(),
        Path(transition.__file__).resolve(),
        Path(upper_arc.__file__).resolve(),
        Path(positive_tail.__file__).resolve(),
        Path(correction.__file__).resolve(),
        Path(recentered_upper_arc.__file__).resolve(),
    )
    artifact = {
        "kind": "rh_diagnostic_upper_event_cell_local_complete_K_T_scout",
        "date": "2026-08-28",
        "status": (
            "diagnostic_local_complete_K_T_derivative_assembled"
            if complete
            else "diagnostic_local_complete_K_T_component_checkpoint_incomplete"
        ),
        "passed": complete,
        "version": VERSION,
        "resource_policy": "one below-normal process, one logical CPU, one numerical thread",
        "offset": args.offset,
        "center_height": center_text,
        "radius": args.radius,
        "precision_bits": args.precision_bits,
        "variant": args.variant,
        "wrapper_radius_caps_bypassed": bool(arb(args.radius) > arb("0.001")),
        "retained_extended_radius_guards": [
            "strict fixed-roster event-cell containment",
            "stationary-root and nonstationary-gap inequalities",
            "phase-sector and compact-tail inequalities",
            "upper-arc analytic and Hardy-operator error budgets",
            "finite-jet, endpoint-recurrence, and integration-by-parts remainder bounds",
            "direct/assembled component and complete-projection overlaps",
        ],
        "new_component_count": new_count,
        "component_cache_path": CACHE.relative_to(ROOT).as_posix(),
        "completed_components": [
            component
            for component, key in zip(COMPONENTS, keys)
            if key in cached and cached[key]["passed"]
        ],
        "assembly": assembly,
        "sources": {
            path.relative_to(ROOT).as_posix(): sha256(path) for path in source_paths
        },
        "diagnostic_boundary": (
            "A complete row encloses (Q_K-T)' only on one local box. It does not connect the "
            "known sign islands, bound a second derivative, cover the event cell, cross an event "
            "wall, or imply RH."
        ),
    }
    atomic_write_json(RESULT, artifact)
    if complete:
        derivative = assembly["Q_K_minus_T_derivative_component_sum"]
        print(
            "assembled local complete K_T: "
            f"derivative_lower={derivative['lower']} derivative_upper={derivative['upper']}",
            flush=True,
        )
    else:
        print(
            f"checkpointed local K_T components: completed={artifact['completed_components']}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
