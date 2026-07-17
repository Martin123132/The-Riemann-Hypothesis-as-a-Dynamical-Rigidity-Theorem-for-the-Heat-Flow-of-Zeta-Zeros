#!/usr/bin/env python3
"""Certify order-seven first-summand curvature on the ray 2<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
import jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate as corridor_core  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate as order6_ray  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7_compact  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
    stable_log_series,
)
from jensen_window_pf_compound_order7_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as HIGH_CUMULANT_SOURCE,
    FINITE_EXACT_CAP as HIGH_CUMULANT_CAP,
)
from jensen_window_pf_compound_order7_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
)


DEFAULT_EXTENSION_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_u2_h11_h12_extension_tiles.jsonl"
)
DEFAULT_RAY_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_u2_u20_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md"
)
MODE_TWO = order6_ray.MODE_TWO
COLLAR_END = order6_ray.COLLAR_END
BASE_TAIL_START = order6_ray.BASE_TAIL_START
RAY_WIDTH = order6_ray.RAY_WIDTH
COLLAR_T = 5
MID_CUMULANT_CAP = order6_ray.HIGH_CUMULANT_CAP
PRECISION_BITS = order6_ray.PRECISION_BITS
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))


@dataclass(frozen=True)
class RayRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def extension_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order6_ray.extension_tasks()


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order6_ray.ray_tasks()


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def extension_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4_compact.PANELS,
        window_y=order4_compact.WINDOW_Y,
        eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
        max_moment=12,
    )
    if not result.get("passed"):
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": result.get("failure"),
        }
    return {
        "kind": "order7_u2_extension_h11_h12_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "h_derivatives": {
            str(order): order4_compact.interval_text(
                result["h_derivatives"][order]
            )
            for order in (11, 12)
        },
    }


def choose_pad(left: Fraction, right: Fraction) -> tuple[Fraction, flint.arb]:
    pad = RAY_WIDTH
    while True:
        candidate = pad / 10
        mode = order6_ray.arb_interval(left - candidate, right + candidate)
        curvature = potential_jet_arb(mode, 2)[2]
        if bool(arb_rational(candidate) * curvature > 2 * COLLAR_T):
            pad = candidate
            continue
        mode = order6_ray.arb_interval(left - pad, right + pad)
        curvature = potential_jet_arb(mode, 2)[2]
        if not bool(arb_rational(pad) * curvature > COLLAR_T):
            raise RuntimeError("adaptive mode pad does not cover the t+-5 collar")
        return pad, curvature


def corridor_h_derivatives(
    left: Fraction,
    right: Fraction,
) -> tuple[dict[int, flint.arb], dict]:
    pad, _ = choose_pad(left, right)
    outer_left = left - pad
    outer_right = right + pad
    mode = order6_ray.arb_interval(outer_left, outer_right)
    jet = potential_jet_arb(mode, 2)
    t, curvature = jet[1], jet[2]
    q_value = flint.arb.pi() * (4 * mode).exp()
    if not bool(
        curvature > 0 and arb_rational(pad) * curvature > COLLAR_T
    ):
        raise RuntimeError("finite-ray mode geometry failed")

    kappa2 = corridor_core.exact_interval(
        1 + arb_rational(order6_ray.CORRIDOR_K2[0]) / q_value,
        1 + arb_rational(order6_ray.CORRIDOR_K2[1]) / q_value,
    )
    derivatives = {
        2: corridor_core.signed_hurwitz_gamma_derivative(
            2,
            t + flint.arb(1) / 2,
        )
        - kappa2 / curvature
    }
    for order, cap in order6_ray.CORRIDOR_CAPS.items():
        baseline = math.factorial(order - 2) * q_value ** (
            1 - flint.arb(order) / 2
        )
        magnitude = corridor_core.exact_interval(
            baseline,
            arb_rational(cap) * baseline,
        )
        cumulant = magnitude if order % 2 == 0 else -magnitude
        derivatives[order] = corridor_core.signed_hurwitz_gamma_derivative(
            order,
            t + flint.arb(1) / 2,
        ) - cumulant / corridor_core.derivative_power(curvature, order)
    for order, cap in (
        (9, MID_CUMULANT_CAP),
        (10, MID_CUMULANT_CAP),
        (11, HIGH_CUMULANT_CAP),
        (12, HIGH_CUMULANT_CAP),
    ):
        magnitude = (
            cap
            * math.factorial(order - 2)
            * q_value ** (1 - flint.arb(order) / 2)
        )
        cumulant = flint.arb(0, magnitude.upper())
        derivatives[order] = corridor_core.signed_hurwitz_gamma_derivative(
            order,
            t + flint.arb(1) / 2,
        ) - cumulant / corridor_core.derivative_power(curvature, order)
    return derivatives, {
        "mode_pad": str(pad),
        "outer_mode": [str(outer_left), str(outer_right)],
        "t_collar_product_lower": arb_lower_text(
            arb_rational(pad) * curvature
        ),
        "cumulant_caps": {
            "9": MID_CUMULANT_CAP,
            "10": MID_CUMULANT_CAP,
            "11": HIGH_CUMULANT_CAP,
            "12": HIGH_CUMULANT_CAP,
        },
    }


def evaluate_dimensionless(
    left: Fraction,
    right: Fraction,
    derivatives: dict[int, flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    central_mode = order6_ray.arb_interval(left, right)
    t = potential_jet_arb(central_mode, 1)[1]
    z = 1 / t
    b = [
        t ** (order + 1)
        * derivatives[order + 2]
        / math.factorial(order)
        for order in range(11)
    ]
    try:
        ell = stable_log_series(series_scale(b, z), 10)
        J = [
            2 * b[order]
            - z * (order + 1) * (order + 2) * ell[order + 2]
            for order in range(9)
        ]
        h = series_add(
            series_scale(ell[:9], 2),
            stable_log_series(series_scale(J, z), 8),
        )
        R = [
            3 * b[order]
            - z * (order + 1) * (order + 2) * h[order + 2]
            for order in range(7)
        ]
        q_layer = series_add(
            series_sub(series_scale(h[:7], 2), ell[:7]),
            stable_log_series(series_scale(R, z), 6),
        )
        S = [
            4 * b[order]
            - z * (order + 1) * (order + 2) * q_layer[order + 2]
            for order in range(5)
        ]
        p = series_add(
            series_sub(series_scale(q_layer[:5], 2), h[:5]),
            stable_log_series(series_scale(S, z), 4),
        )
        T = [
            5 * b[order]
            - z * (order + 1) * (order + 2) * p[order + 2]
            for order in range(3)
        ]
        r_layer = series_add(
            series_sub(series_scale(p[:3], 2), q_layer[:3]),
            stable_log_series(series_scale(T, z), 2),
        )
    except Exception as exc:
        return {"passed": False, "failure": "nested-series", "detail": repr(exc)}

    scaled = 2 * r_layer[2]
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    coordinates = {"J": J[0], "R": R[0], "S": S[0], "T": T[0]}
    if not bool(
        all(value > 0 for value in coordinates.values()) and margin > 0
    ):
        return {
            "passed": False,
            "failure": "dimensionless-margin",
            "mode": [str(left), str(right)],
            **{
                f"{name}_ball": value.str(35).replace("e", "E")
                for name, value in coordinates.items()
            },
            "scaled_ball": scaled.str(35).replace("e", "E"),
        }
    return {
        "passed": True,
        "failure": None,
        "mode": [str(left), str(right)],
        "t_ball": t.str(35).replace("e", "E"),
        **{
            f"{name}_lower": arb_lower_text(value)
            for name, value in coordinates.items()
        },
        "scaled_curvature_upper": arb_upper_text(scaled),
        "scaled_margin_lower": arb_lower_text(margin),
        "diagnostics": diagnostics,
    }


def ray_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    try:
        derivatives, diagnostics = corridor_h_derivatives(left, right)
        result = evaluate_dimensionless(
            left,
            right,
            derivatives,
            diagnostics=diagnostics,
        )
    except Exception as exc:
        result = {"passed": False, "failure": "exception", "detail": repr(exc)}
    return {
        "kind": "order7_nested_curvature_dimensionless_finite_ray_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def load_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError(f"cache {path.name} has too many rows")
    for record, (index, left, right) in zip(records, tasks):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"cache prefix mismatch at {path.name}:{index}")
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    worker,
    *,
    workers: int,
    overwrite: bool,
    max_rows: int | None = None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, tasks)
    stop = len(tasks) if max_rows is None else min(len(tasks), max_rows)
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    if workers == 1:
        iterator = map(worker, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
        )
        iterator = executor.map(worker, remaining, chunksize=16)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"cache row {record.get('index')} failed: "
                        f"{record.get('failure')} {record.get('detail', '')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 500 == 0:
                    handle.flush()
                    print(
                        f"{path.stem}: {len(records)}/{stop} "
                        f"({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def compact_tail_records() -> list[dict]:
    tasks = order7_compact.deterministic_tasks()
    low = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
    middle = order6_compact.load_extension_cache(
        order6_compact.DEFAULT_EXTENSION_CACHE,
        tasks,
    )
    high = order7_compact.load_extension_cache(
        order7_compact.DEFAULT_EXTENSION_CACHE,
        tasks,
    )
    rows = []
    for low_row, middle_row, high_row in zip(low, middle, high):
        if Fraction(low_row["mode_left"]) < BASE_TAIL_START:
            continue
        merged = dict(low_row)
        merged["h_derivatives"] = {
            **low_row["h_derivatives"],
            **middle_row["h_derivatives"],
            **high_row["h_derivatives"],
        }
        rows.append(merged)
    if not rows or rows[0]["mode_left"] != str(BASE_TAIL_START):
        raise RuntimeError("order-seven compact tail starts at the wrong mode")
    return rows


def merged_extension(top: list[dict]) -> list[dict]:
    tasks = extension_tasks()
    low = order6_ray.load_cache(order6_ray.DEFAULT_EXTENSION_CACHE, tasks)
    if len(low) != len(tasks) or len(top) != len(tasks):
        raise RuntimeError("u=2 extension caches are incomplete")
    rows = []
    for low_row, top_row in zip(low, top):
        if (
            low_row["mode_left"] != top_row["mode_left"]
            or low_row["mode_right"] != top_row["mode_right"]
        ):
            raise RuntimeError("u=2 extension caches are not aligned")
        merged = dict(low_row)
        merged["h_derivatives"] = {
            **low_row["h_derivatives"],
            **top_row["h_derivatives"],
        }
        rows.append(merged)
    return rows


def initial_collar(extension: list[dict]) -> dict:
    flint.ctx.prec = PRECISION_BITS
    records = compact_tail_records() + merged_extension(extension)
    derivatives = {
        order: order4_compact.hull(
            [
                order4_compact.interval_from_text(
                    row["h_derivatives"][str(order)]
                )
                for row in records
            ]
        )
        for order in range(2, 13)
    }
    left_t = order4_compact.interval_from_text(records[0]["t_left"])
    right_t = order4_compact.interval_from_text(records[-1]["t_right"])
    central_left_t = potential_jet_arb(arb_rational(MODE_TWO), 1)[1]
    central_right_t = potential_jet_arb(arb_rational(COLLAR_END), 1)[1]
    if not bool(
        left_t < central_left_t - COLLAR_T
        and right_t > central_right_t + COLLAR_T
    ):
        raise RuntimeError("u=2 extension cache does not contain a t+-5 collar")
    result = evaluate_dimensionless(
        MODE_TWO,
        COLLAR_END,
        derivatives,
        diagnostics={
            "record_count": len(records),
            "left_t_collar": (central_left_t - left_t)
            .str(35)
            .replace("e", "E"),
            "right_t_collar": (right_t - central_right_t)
            .str(35)
            .replace("e", "E"),
        },
    )
    if not result.get("passed"):
        raise RuntimeError(f"initial order-seven u=2 collar failed: {result}")
    return result


def source_contract() -> dict:
    compact = load_json(order7_compact.DEFAULT_OUT)
    high = load_json(HIGH_CUMULANT_SOURCE)
    if compact.get("compact", {}).get("all_blocks_passed") is not True:
        raise RuntimeError("order-seven compact handoff changed")
    if high.get("exact", {}).get("finite_exact_corridor_cap") != HIGH_CUMULANT_CAP:
        raise RuntimeError("order-seven high-cumulant source changed")
    order6_sources = order6_ray.source_contract()
    return {
        "exact_corridors": order6_sources["exact_corridors"],
        "mid_cumulants": order6_sources["high_cumulants"],
        "high_cumulants": high["exact"]["finite_exact_corridor"],
        "compact_handoff": compact["compact"]["theorem"],
        "order6_extension_sha256": sha256(order6_ray.DEFAULT_EXTENSION_CACHE),
        "order7_high_cumulant_sha256": sha256(HIGH_CUMULANT_SOURCE),
        "order7_compact_sha256": sha256(order7_compact.DEFAULT_OUT),
    }


def finite_certificate(extension: list[dict], ray: list[dict]) -> dict:
    collar = initial_collar(extension)
    all_rows = [collar, *ray]
    largest = max(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["scaled_curvature_upper"]),
    )
    weakest_t = min(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["T_lower"]),
    )
    selected_indices = sorted(
        {
            0,
            len(all_rows) // 4,
            len(all_rows) // 2,
            3 * len(all_rows) // 4,
            len(all_rows) - 1,
            largest,
            weakest_t,
        }
    )
    return {
        "mode_range": ["2", "20"],
        "block_width": str(RAY_WIDTH),
        "initial_collar": collar,
        "ray_blocks": len(ray),
        "all_blocks_passed": True,
        "largest_scaled_index": largest,
        "largest_scaled_curvature_upper": all_rows[largest][
            "scaled_curvature_upper"
        ],
        "weakest_T_index": weakest_t,
        "weakest_T_lower": all_rows[weakest_t]["T_lower"],
        "selected": [all_rows[index] for index in selected_indices],
        "theorem": "r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20",
    }


def build_artifact(extension: list[dict], ray: list[dict]) -> dict:
    sources = source_contract()
    finite = finite_certificate(extension, ray)
    rows = [
        RayRow(
            "co7ncfr_01_exact_derivatives",
            "exact_theorem_composition",
            "ready_to_apply",
            "Exact and coarse cumulant corridors through order twelve enclose every required H derivative on the finite ray.",
            "H^(r) enclosed for r=2,...,12 on a t+-5 collar",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co7ncfr_02_initial_collar",
            "interval_certificate",
            "ready_to_apply",
            "The extended compact quadrature cache bridges mode u=2 without using a corridor below its domain.",
            "r_1''(t)<=600/t^2 for 2<=u<=2001/1000",
            "Finite aligned H2-H12 cache only.",
            finite["initial_collar"],
        ),
        RayRow(
            "co7ncfr_03_dimensionless_cover",
            "interval_theorem",
            "ready_to_apply",
            "A rational 1/1000-mode grid and dimensionless stable arithmetic prove the finite-ray curvature theorem.",
            finite["theorem"],
            "Outward-rounded Arb cover; no point sampling.",
            {
                "blocks": finite["ray_blocks"],
                "largest_scaled_upper": finite[
                    "largest_scaled_curvature_upper"
                ],
                "weakest_T_lower": finite["weakest_T_lower"],
            },
        ),
        RayRow(
            "co7ncfr_04_asymptotic_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the dimensionless theorem on the asymptotic mode ray.",
            "r_1''(t)<=600/t^2 for u>=20",
            "Open asymptotic ray only.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate",
        "date": "2026-07-14",
        "status": "rigorous order-seven nested first-summand curvature theorem on 2<=u<=20",
        "proof_boundary": (
            "This artifact proves continuous first-summand curvature only on "
            "2<=u<=20. It does not prove the u>=20 ray, order-seven entry, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "finite_ray": finite,
        "cache": {
            "extension_path": DEFAULT_EXTENSION_CACHE.relative_to(
                REPO_ROOT
            ).as_posix(),
            "extension_rows": len(extension),
            "extension_sha256": sha256(DEFAULT_EXTENSION_CACHE),
            "ray_path": DEFAULT_RAY_CACHE.relative_to(REPO_ROOT).as_posix(),
            "ray_rows": len(ray),
            "ray_sha256": sha256(DEFAULT_RAY_CACHE),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 3,
            "open_rows": 1,
            "initial_collar_blocks": 1,
            "finite_ray_blocks": len(ray),
            "finite_ray_theorems": 1,
            "open_asymptotic_rays": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_ray"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous first-summand order-seven curvature theorem on",
        "`2<=u<=20`. This is not a proof of the asymptotic ray, order-seven",
        "entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Dimensionless Cover",
        "",
        "Exact signed cumulant corridors through order eight, coarse global",
        "corridors for orders nine and ten, and the finite 14001 corridors for",
        "orders eleven and twelve enclose `H^(2),...,H^(12)` on a strict",
        "`t+-5` collar. A short aligned quadrature extension handles `u=2`;",
        "the remaining interval uses rational mode width `1/1000`.",
        "",
        "```text",
        finite["theorem"],
        f"ray blocks={finite['ray_blocks']},",
        f"largest t^2*r_1'' upper={finite['largest_scaled_curvature_upper']},",
        f"weakest T_1 lower={finite['weakest_T_lower']}.",
        "```",
        "",
        "## Remaining Ray",
        "",
        "```text",
        "Prove r_1''(t)<=600/t^2 for u>=20.",
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension-cache", type=Path, default=DEFAULT_EXTENSION_CACHE)
    parser.add_argument("--ray-cache", type=Path, default=DEFAULT_RAY_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-caches", action="store_true")
    parser.add_argument("--max-ray-blocks", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()
    workers = max(1, args.workers)
    extension = build_cache(
        args.extension_cache,
        extension_tasks(),
        extension_task,
        workers=workers,
        overwrite=args.overwrite_caches,
    )
    tasks = ray_tasks()
    ray = build_cache(
        args.ray_cache,
        tasks,
        ray_task,
        workers=workers,
        overwrite=args.overwrite_caches,
        max_rows=args.max_ray_blocks,
    )
    print(
        f"order-seven finite-ray caches: extension={len(extension)}, "
        f"ray={len(ray)}/{len(tasks)}"
    )
    if args.cache_only or len(ray) != len(tasks):
        return 0
    artifact = build_artifact(extension, ray)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven nested finite-ray certificate: "
        f"{summary['finite_ray_blocks']} ray blocks, "
        f"{summary['finite_ray_theorems']} theorem, "
        f"{summary['open_asymptotic_rays']} open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
