#!/usr/bin/env python3
"""Certify beta^-4 transport throughout all 399 exact event cells."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, wait
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any, Iterable, Iterator


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate"
ATLAS_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_atlas_gate"
ATLAS_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{ATLAS_STEM}.py"
ATLAS_RESULT = REPO_ROOT / f"work/rh_compute/results/{ATLAS_STEM}.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
CACHE = REPO_ROOT / f"work/rh_compute/results/cache/{STEM}_rows.jsonl"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 70
EVENT_COUNT = 399
MOMENT_ORDER = 6
CONFIG_VERSION = "shifted-moments-v1"
TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def segment_ball(center: acb, half_width: arb, direction: acb) -> acb:
    if direction.imag.contains(0):
        return acb(arb(center.real, half_width), center.imag)
    return acb(center.real, arb(center.imag, half_width))


def contour_panels(atlas, module, face) -> Iterator[tuple[str, acb, arb, acb, int]]:
    real_panels, _ = atlas.shifted_real_panels(face)
    for center, half_width in real_panels:
        yield "real_core", acb(atlas.fraction_arb(center)), atlas.fraction_arb(half_width), acb(1), 1
    for side, orientation in ((1, 1), (-1, -1)):
        for index in range(module.VERTICAL_PANELS):
            center_y = module.VERTICAL_HALF_WIDTH + 2 * module.VERTICAL_HALF_WIDTH * index
            label = "right_connector" if side > 0 else "left_connector"
            yield label, acb(side * module.CONTOUR_RADIUS, center_y), module.VERTICAL_HALF_WIDTH, face.i, orientation
    for side in (-1, 1):
        start = -module.HORIZONTAL_END if side < 0 else module.CONTOUR_RADIUS
        for index in range(module.HORIZONTAL_PANELS_PER_SIDE):
            center_x = start + module.REAL_HALF_WIDTH + 2 * module.REAL_HALF_WIDTH * index
            yield "horizontal_near", acb(center_x, module.CONTOUR_HEIGHT), module.REAL_HALF_WIDTH, acb(1), 1


def panel_error_bases(module, face, center: acb, half_width: arb, direction: acb) -> tuple[arb, arb, arb]:
    is_lifted_horizontal = direction.imag.contains(0) and center.imag.lower() > arb("0.9")
    central = False
    if is_lifted_horizontal:
        disk_multiplier = arb("1.2")
    else:
        q0 = center / face.beta
        tanh0 = -face.i * (face.i * q0).tan()
        exact_b0 = -(face.detuning + face.beta * tanh0)
        canonical_b0 = -(center + face.detuning)
        central = min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.002")
        disk_multiplier = arb("1.8" if central else "3.0")

    disk_radius = disk_multiplier * half_width
    z_disk = acb(arb(center.real, disk_radius), arb(center.imag, disk_radius))
    q_disk = z_disk / face.beta
    tanh_disk, _ = face.tanh_ball(q_disk)
    exact_b_disk = -(face.detuning + face.beta * tanh_disk)
    canonical_b_disk = -(z_disk + face.detuning)
    use_removable = not (exact_b_disk.abs_lower() > 0 and canonical_b_disk.abs_lower() > 0)
    if use_removable and not central and not is_lifted_horizontal:
        disk_multiplier = arb("1.8")
        disk_radius = disk_multiplier * half_width
        z_disk = acb(arb(center.real, disk_radius), arb(center.imag, disk_radius))
        q_disk = z_disk / face.beta
        tanh_disk, _ = face.tanh_ball(q_disk)
        exact_b_disk = -(face.detuning + face.beta * tanh_disk)
        canonical_b_disk = -(z_disk + face.detuning)
        use_removable = not (exact_b_disk.abs_lower() > 0 and canonical_b_disk.abs_lower() > 0)
    if use_removable:
        q0 = center / face.beta
        tanh0 = -face.i * (face.i * q0).tan()
        exact_b0 = -(face.detuning + face.beta * tanh0)
        canonical_b0 = -(center + face.detuning)
        require(
            min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.002"),
            "a noncentral moment panel entered the removable chart",
        )

    inner_error = face.inner_truncation_bound(z_disk, use_removable)
    maximum_base = face.true_integrand_bound(z_disk) + inner_error
    ratio = half_width / disk_radius
    outer_base = 2 * maximum_base * half_width * ratio**module.SERIES_ORDER / (
        (module.SERIES_ORDER + 1) * (1 - ratio)
    )
    inner_base = 2 * half_width * inner_error
    return outer_base, abs(z_disk).upper(), inner_base * arb(1)


def certified_moments(atlas, module, face, order: int) -> tuple[list[acb], arb, dict[str, arb]]:
    moments = [acb(0) for _ in range(order + 1)]
    errors = [arb(0) for _ in range(order + 1)]
    absolute_mass = arb(0)
    mass_by_piece: dict[str, arb] = {}
    for piece, center, half_width, direction, orientation in contour_panels(atlas, module, face):
        base = face.difference_series(center, direction)
        variable = acb_series([center, direction], prec=module.SERIES_ORDER)
        power = acb_series([1], prec=module.SERIES_ORDER)
        outer_base, disk_power_base, inner_base = panel_error_bases(module, face, center, half_width, direction)
        path_power_base = abs(segment_ball(center, half_width, direction)).upper()
        disk_power = arb(1)
        path_power = arb(1)
        for degree in range(order + 1):
            moments[degree] += orientation * module.integrate_symmetric_series(base * power, half_width)
            errors[degree] += outer_base * disk_power + inner_base * path_power
            power *= variable
            disk_power *= disk_power_base
            path_power *= path_power_base
        segment = segment_ball(center, half_width, direction)
        panel_mass = 2 * half_width * face.true_integrand_bound(segment)
        absolute_mass += panel_mass
        mass_by_piece[piece] = mass_by_piece.get(piece, arb(0)) + panel_mass
    return [module.add_error(value, error) for value, error in zip(moments, errors)], absolute_mass, mass_by_piece


def polynomial_value(coefficients: list[acb], theta: int) -> acb:
    total = acb(0)
    power = arb(1)
    for coefficient in coefficients:
        total += power * coefficient
        power *= theta
    return total


def integrate_height_event(atlas, second, first, module, face_class, event_index: int, theta: int) -> acb:
    mode = atlas.mode_for_event(event_index)
    face = first.event_face(face_class, module, event_index, mode)
    h = arb(theta) * face.pi / 16
    face.t = face.tau + h
    face.eta = face.tstar / face.t
    face.lam = (face.eta - 1) * face.t / face.beta
    core, _, _ = atlas.shifted_real_core(module, face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far = second.second_order_far_bound(first, module, face)
    return module.add_error(core + right - left + horizontal, far)


def integrate_event(atlas, second, first, module, face_class, center_rows: list[dict[str, Any]], event_index: int) -> dict[str, Any]:
    mode = atlas.mode_for_event(event_index)
    face = first.event_face(face_class, module, event_index, mode)
    started = time.perf_counter()
    moments, absolute_mass, mass_by_piece = certified_moments(atlas, module, face, MOMENT_ORDER)
    kappa = face.pi / (16 * face.beta)
    coefficients = [
        (face.i * kappa) ** degree / arb(math.factorial(degree)) * moment
        for degree, moment in enumerate(moments)
    ]
    z_max = (module.HORIZONTAL_END**2 + module.CONTOUR_HEIGHT**2).sqrt()
    expansion_size = kappa * z_max
    finite_remainder = (
        absolute_mass
        * expansion_size.exp()
        * expansion_size ** (MOMENT_ORDER + 1)
        / arb(math.factorial(MOMENT_ORDER + 1))
    )
    far_at_event = second.second_order_far_bound(first, module, face)
    uniform_far = kappa.exp() * far_at_event
    uniform_bound = sum((abs(value).upper() for value in coefficients), arb(0)) + finite_remainder + uniform_far
    physical_bound = arb(2).sqrt() / face.pi * uniform_bound
    require(uniform_bound < TARGET_NORMALIZED, f"event {event_index} uniform normalized target failed with {uniform_bound}")
    require(physical_bound < TARGET_PHYSICAL, f"event {event_index} uniform physical target failed with {physical_bound}")

    event_reconstruction = module.add_error(moments[0], far_at_event)
    stored_event = parse_complex(center_rows[event_index]["normalized_second_order_difference_ball"])
    require(event_reconstruction.real.overlaps(stored_event.real), f"event {event_index} real-center overlap failed")
    require(event_reconstruction.imag.overlaps(stored_event.imag), f"event {event_index} imaginary-center overlap failed")
    endpoint_balls = []
    for theta in (-1, 1):
        endpoint_balls.append({
            "theta": theta,
            "transport_ball": complex_record(
                module.add_error(polynomial_value(coefficients, theta), finite_remainder + uniform_far)
            ),
        })
    return {
        "config_version": CONFIG_VERSION,
        "event_index": event_index,
        "mode": mode,
        "signed_odd": atlas.signed_odd_for_event(event_index),
        "detuning_ball": face.detuning.str(PRECISION, more=True),
        "moments": [
            {
                "degree": degree,
                "moment_ball": complex_record(moment),
                "transport_coefficient_ball": complex_record(coefficients[degree]),
            }
            for degree, moment in enumerate(moments)
        ],
        "finite_contour_absolute_mass_bound": absolute_mass.str(PRECISION, more=True),
        "absolute_mass_by_piece": {key: value.str(PRECISION, more=True) for key, value in sorted(mass_by_piece.items())},
        "maximum_transport_argument": expansion_size.str(PRECISION, more=True),
        "finite_exponential_remainder_bound": finite_remainder.str(PRECISION, more=True),
        "event_far_tail_absolute_bound": far_at_event.str(PRECISION, more=True),
        "uniform_transported_far_tail_bound": uniform_far.str(PRECISION, more=True),
        "uniform_normalized_absolute_bound": uniform_bound.str(PRECISION, more=True),
        "uniform_normalized_absolute_upper": uniform_bound.upper().str(PRECISION, more=True),
        "uniform_physical_absolute_bound": physical_bound.str(PRECISION, more=True),
        "uniform_physical_absolute_upper": physical_bound.upper().str(PRECISION, more=True),
        "event_center_overlap": True,
        "endpoint_transport_balls": endpoint_balls,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def load_cache() -> dict[int, dict[str, Any]]:
    if not CACHE.is_file():
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line_number, line in enumerate(CACHE.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        row = json.loads(line)
        require(row.get("config_version") == CONFIG_VERSION, f"cache config mismatch on line {line_number}")
        index = int(row["event_index"])
        require(index not in rows, f"duplicate cached event {index}")
        rows[index] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def parse_events(value: str) -> list[int]:
    events = [int(part.strip()) for part in value.split(",") if part.strip()]
    require(events and all(0 <= index < EVENT_COUNT for index in events), "invalid event selection")
    return events


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# All-event continuous-height second-order atlas

Date: 2026-08-12

Status: rigorous beta-minus-four transport through all 399 exact event cells;
not a proof of the ordinary-mode join, complete T_upper, or RH

For each exact event center `tau_n`, both the exact and beta-minus-four
corrected common-chart integrands satisfy

```text
D_(n,tau_n+h)(z)=exp(i h z/beta)D_(n,tau_n)(z).
```

Writing `h=(pi/16)theta`, `|theta|<=1`, and
`kappa=pi/(16 beta)`, the finite-contour integral is enclosed by moments
through degree `{artifact['method']['moment_order']}` and an explicit
exponential remainder.  The lifted far tails acquire at most `exp(kappa)`.

All `{c['event_count']}` exact event cells are certified.  The largest
normalized upper bound is `{c['maximum_normalized_absolute_upper']}` at event
`{c['maximum_normalized_event']}`; the corresponding all-cell physical upper
bound is below `{c['maximum_physical_absolute_upper']}`.  These are below the
targets `0.000019` and `0.0000086`.

This closes the corrected lower-fold comparison only on the 399 listed event
cells.  The ordinary interior-mode join and full `T_upper` comparison remain
open.  `Lambda<=0`, PF-infinity, RH, and a prize-level conclusion remain
unproved.
"""


def finalize(rows: list[dict[str, Any]], priority: str, workers: int, elapsed: float) -> None:
    require(len(rows) == EVENT_COUNT, "cannot finalize an incomplete height atlas")
    require([row["event_index"] for row in rows] == list(range(EVENT_COUNT)), "height cache is not ordered")
    norm_max = max(rows, key=lambda row: float(arb(row["uniform_normalized_absolute_upper"])))
    phys_max = max(rows, key=lambda row: float(arb(row["uniform_physical_absolute_upper"])))
    artifact: dict[str, Any] = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "beta_minus_4_continuous_height_certified_on_all_399_exact_event_cells",
        "resource_policy": {"workers": workers, "process_priority": priority, "resumable_append_only_cache": True},
        "method": {
            "exact_transport_identity": "D_(n,tau_n+h)=exp(i*h*z/beta)*D_(n,tau_n)",
            "theta_interval": "[-1,1]",
            "moment_order": MOMENT_ORDER,
            "moment_definition": "M_(n,k)=int_Gamma z^k D_(n,tau_n)(z) dz",
            "far_tail_transport_factor": "exp(pi/(16*beta)) on Im(z)=1",
            "event_center_partition": "w=z+d_n denominator-2048 rational panel",
        },
        "certificate": {
            "event_count": EVENT_COUNT,
            "maximum_normalized_event": norm_max["event_index"],
            "maximum_normalized_absolute_upper": norm_max["uniform_normalized_absolute_upper"],
            "maximum_physical_event": phys_max["event_index"],
            "maximum_physical_absolute_upper": phys_max["uniform_physical_absolute_upper"],
            "normalized_target": "0.000019",
            "physical_target": "0.0000086",
            "all_rows": rows,
        },
        "claims": {
            "exact_common_height_transport_identity_proved": True,
            "all_399_exact_event_cells_uniformly_certified": True,
            "normalized_target_proved_on_all_cells": True,
            "physical_target_proved_on_all_cells": True,
            "ordinary_mode_join_proved": False,
            "complete_T_upper_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": "The 399 exact lower-fold event cells only. No ordinary-mode join, complete T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
        "next_action": "Join the certified lower-fold event cells to ordinary interior modes 622..39852 without overlap or omission, then assemble the complete T_upper comparison.",
        "dependencies": {
            "all_event_center_atlas": {"path": relative(ATLAS_RESULT), "sha256": file_hash(ATLAS_RESULT)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "row_cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "runtime": {"elapsed_seconds_this_invocation": round(elapsed, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def load_sources(worker_name: str):
    atlas = load_module(ATLAS_BUILDER, f"height_atlas_{worker_name}")
    second = atlas.load_module(atlas.SECOND_BUILDER, f"height_second_{worker_name}")
    first = second.load_first()
    first.set_low_priority()
    flint.ctx.dps = PRECISION
    module = first.load_full_module()
    atlas.configure(module)
    face_class = atlas.make_atlas_face(second, first, module)
    center_artifact = json.loads(ATLAS_RESULT.read_text(encoding="utf-8"))
    center_rows = center_artifact["certificate"]["all_rows"]
    return atlas, second, first, module, face_class, center_rows


def run_events(indices: Iterable[int]) -> list[dict[str, Any]]:
    state = load_sources(f"scout_{os.getpid()}")
    rows = []
    for index in indices:
        row = integrate_event(*state, index)
        rows.append(row)
        print(
            f"event {index:03d} uniform<={row['uniform_normalized_absolute_upper']} "
            f"mass={row['finite_contour_absolute_mass_bound']} elapsed={row['elapsed_seconds']}s",
            flush=True,
        )
    return rows


_WORKER_STATE: tuple[Any, Any, Any, Any, Any, list[dict[str, Any]]] | None = None


def initialize_worker() -> None:
    global _WORKER_STATE
    _WORKER_STATE = load_sources(f"worker_{os.getpid()}")


def compute_worker(event_index: int) -> dict[str, Any]:
    require(_WORKER_STATE is not None, "height worker was not initialized")
    return integrate_event(*_WORKER_STATE, event_index)


def run_parallel(
    pending: list[int],
    cached: dict[int, dict[str, Any]],
    workers: int,
    started: float,
    runtime_limit: float,
) -> None:
    try:
        import psutil
    except Exception:  # pragma: no cover
        psutil = None
    high_cpu_samples = 0
    park = False
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize_worker) as executor:
        for offset in range(0, len(pending), workers):
            if time.perf_counter() - started >= runtime_limit:
                break
            indices = pending[offset : offset + workers]
            futures = [executor.submit(compute_worker, index) for index in indices]
            incomplete = set(futures)
            while incomplete:
                _, incomplete = wait(incomplete, timeout=5)
                if incomplete and psutil is not None:
                    total_cpu = psutil.cpu_percent(interval=0.2)
                    high_cpu_samples = high_cpu_samples + 1 if total_cpu > 75 else 0
                    if high_cpu_samples >= 2:
                        park = True
            rows = sorted((future.result() for future in futures), key=lambda row: row["event_index"])
            for row in rows:
                append_cache(row)
                cached[row["event_index"]] = row
                print(
                    f"event {row['event_index']:03d} uniform<={row['uniform_normalized_absolute_upper']} "
                    f"elapsed={row['elapsed_seconds']}s",
                    flush=True,
                )
            if park:
                break


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", help="comma-separated scout events; do not modify the canonical cache")
    parser.add_argument("--reset-cache", action="store_true")
    parser.add_argument("--runtime-limit-seconds", type=float, default=10_800)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    args = parser.parse_args()
    require(ATLAS_BUILDER.is_file() and ATLAS_RESULT.is_file() and CHECKER.is_file(), "missing dependency or checker")

    if args.events:
        run_events(parse_events(args.events))
        return

    atlas = load_module(ATLAS_BUILDER, "height_atlas_controller")
    second = atlas.load_module(atlas.SECOND_BUILDER, "height_second_controller")
    first = second.load_first()
    priority = first.set_low_priority()
    if args.reset_cache and CACHE.exists():
        CACHE.unlink()
    cached = load_cache()
    started = time.perf_counter()
    pending = [index for index in range(EVENT_COUNT) if index not in cached]
    if args.workers == 1:
        state = load_sources(f"serial_{os.getpid()}")
        for index in pending:
            if time.perf_counter() - started >= args.runtime_limit_seconds:
                break
            row = integrate_event(*state, index)
            append_cache(row)
            cached[index] = row
            print(
                f"event {index:03d} uniform<={row['uniform_normalized_absolute_upper']} "
                f"elapsed={row['elapsed_seconds']}s",
                flush=True,
            )
    else:
        run_parallel(pending, cached, args.workers, started, args.runtime_limit_seconds)

    rows = [cached[index] for index in sorted(cached)]
    if len(rows) == EVENT_COUNT:
        finalize(rows, priority, args.workers, time.perf_counter() - started)
        print(f"certified continuous height on all {EVENT_COUNT} exact event cells", flush=True)
    else:
        print(
            f"checkpointed {len(rows)}/{EVENT_COUNT} height rows in {relative(CACHE)}; rerun to resume",
            flush=True,
        )


if __name__ == "__main__":
    main()
