#!/usr/bin/env python3
"""Certify the order-six nested first-summand curvature on 2<=u<=20."""

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
from jensen_window_pf_compound_order4_gaussian_cumulant_ray_target import (  # noqa: E402
    CORRIDOR_CAPS,
    CORRIDOR_K2,
)
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
    stable_log_series,
)
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as compact_source  # noqa: E402
from jensen_window_pf_compound_order6_nested_curvature_interval_core import (  # noqa: E402
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
    "jensen_window_pf_compound_order6_nested_curvature_u2_extension_tiles.jsonl"
)
DEFAULT_RAY_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_u2_u20_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md"
)
SOURCE_EXACT_CORRIDORS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json"
)
SOURCE_HIGH_CUMULANTS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.json"
)
SOURCE_COMPACT = compact_source.DEFAULT_OUT
MODE_TWO = Fraction(2)
COLLAR_END = Fraction(2001, 1000)
EXTENSION_START = order4_compact.OUTER_END
EXTENSION_END = Fraction(20011, 10000)
BASE_TAIL_START = Fraction(19999, 10000)
RAY_END = Fraction(20)
RAY_WIDTH = Fraction(1, 1000)
COLLAR_T = 4
HIGH_CUMULANT_CAP = 50000
PRECISION_BITS = 384
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
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def arb_interval(left: Fraction, right: Fraction) -> flint.arb:
    lower = arb_rational(left)
    upper = arb_rational(right)
    return (lower + upper) / 2 + flint.arb(0, (upper - lower) / 2)


def source_contract() -> dict:
    exact = load_json(SOURCE_EXACT_CORRIDORS)
    high = load_json(SOURCE_HIGH_CUMULANTS)
    compact = load_json(SOURCE_COMPACT)
    if exact.get("summary", {}).get("global_exact_corridors_closed") is not True:
        raise RuntimeError("exact cumulant corridors through order eight changed")
    if high.get("exact", {}).get("exact_corridor_cap") != HIGH_CUMULANT_CAP:
        raise RuntimeError("high-cumulant corridor source changed")
    if compact.get("compact", {}).get("certified_t_range") != "321<=t<=V'(2)":
        raise RuntimeError("order-six compact handoff changed")
    return {
        "exact_corridors": "signed exact cumulant corridors r=2,...,8 on u>=2",
        "high_cumulants": high["exact"]["exact_corridor"],
        "compact_handoff": compact["compact"]["certified_t_range"],
    }


def extension_tasks() -> list[tuple[int, Fraction, Fraction]]:
    tasks = []
    left = EXTENSION_START
    index = 0
    while left < EXTENSION_END:
        right = min(left + order4_compact.TILE_WIDTH, EXTENSION_END)
        tasks.append((index, left, right))
        left = right
        index += 1
    return tasks


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    tasks = []
    left = COLLAR_END
    index = 0
    while left < RAY_END:
        right = min(left + RAY_WIDTH, RAY_END)
        tasks.append((index, left, right))
        left = right
        index += 1
    return tasks


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
        max_moment=10,
    )
    if not result.get("passed"):
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": result.get("failure"),
        }
    t_left = potential_jet_arb(arb_rational(left), 1)[1]
    t_right = potential_jet_arb(arb_rational(right), 1)[1]
    return {
        "kind": "order6_u2_extension_h2_h10_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "t_left": order4_compact.interval_text(t_left),
        "t_right": order4_compact.interval_text(t_right),
        "h_derivatives": {
            str(order): order4_compact.interval_text(result["h_derivatives"][order])
            for order in range(2, 11)
        },
        "normalizer_lower": result["normalizer_lower"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
    }


def choose_pad(left: Fraction, right: Fraction) -> tuple[Fraction, flint.arb]:
    pad = RAY_WIDTH
    while True:
        candidate = pad / 10
        mode = arb_interval(left - candidate, right + candidate)
        curvature = potential_jet_arb(mode, 2)[2]
        if bool(arb_rational(candidate) * curvature > 8):
            pad = candidate
            continue
        mode = arb_interval(left - pad, right + pad)
        curvature = potential_jet_arb(mode, 2)[2]
        if not bool(arb_rational(pad) * curvature > COLLAR_T):
            raise RuntimeError("adaptive mode pad does not cover the t+-4 collar")
        return pad, curvature


def corridor_h_derivatives(
    left: Fraction, right: Fraction
) -> tuple[dict[int, flint.arb], dict]:
    pad, _ = choose_pad(left, right)
    outer_left = left - pad
    outer_right = right + pad
    mode = arb_interval(outer_left, outer_right)
    jet = potential_jet_arb(mode, 10)
    t, curvature = jet[1], jet[2]
    q = flint.arb.pi() * (4 * mode).exp()
    if not bool(curvature > 0 and arb_rational(pad) * curvature > COLLAR_T):
        raise RuntimeError("finite-ray mode geometry failed")

    kappa2 = corridor_core.exact_interval(
        1 + arb_rational(CORRIDOR_K2[0]) / q,
        1 + arb_rational(CORRIDOR_K2[1]) / q,
    )
    derivatives = {
        2: corridor_core.signed_hurwitz_gamma_derivative(
            2, t + flint.arb(1) / 2
        )
        - kappa2 / curvature
    }
    for order, cap in CORRIDOR_CAPS.items():
        baseline = math.factorial(order - 2) * q ** (
            1 - flint.arb(order) / 2
        )
        magnitude = corridor_core.exact_interval(
            baseline, arb_rational(cap) * baseline
        )
        cumulant = magnitude if order % 2 == 0 else -magnitude
        derivatives[order] = corridor_core.signed_hurwitz_gamma_derivative(
            order, t + flint.arb(1) / 2
        ) - cumulant / corridor_core.derivative_power(curvature, order)
    for order in (9, 10):
        magnitude = (
            HIGH_CUMULANT_CAP
            * math.factorial(order - 2)
            * q ** (1 - flint.arb(order) / 2)
        )
        cumulant = flint.arb(0, magnitude.upper())
        derivatives[order] = corridor_core.signed_hurwitz_gamma_derivative(
            order, t + flint.arb(1) / 2
        ) - cumulant / corridor_core.derivative_power(curvature, order)
    return derivatives, {
        "mode_pad": str(pad),
        "outer_mode": [str(outer_left), str(outer_right)],
        "t_collar_product_lower": arb_lower_text(arb_rational(pad) * curvature),
    }


def evaluate_dimensionless(
    left: Fraction,
    right: Fraction,
    derivatives: dict[int, flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    central_mode = arb_interval(left, right)
    t = potential_jet_arb(central_mode, 1)[1]
    z = 1 / t
    b = [
        t ** (order + 1) * derivatives[order + 2] / math.factorial(order)
        for order in range(9)
    ]
    try:
        ell = stable_log_series(series_scale(b, z), 8)
        J = [
            2 * b[order]
            - z * (order + 1) * (order + 2) * ell[order + 2]
            for order in range(7)
        ]
        h = series_add(
            series_scale(ell[:7], 2),
            stable_log_series(series_scale(J, z), 6),
        )
        R = [
            3 * b[order]
            - z * (order + 1) * (order + 2) * h[order + 2]
            for order in range(5)
        ]
        q = series_add(
            series_sub(series_scale(h[:5], 2), ell[:5]),
            stable_log_series(series_scale(R, z), 4),
        )
        S = [
            4 * b[order]
            - z * (order + 1) * (order + 2) * q[order + 2]
            for order in range(3)
        ]
        p = series_add(
            series_sub(series_scale(q[:3], 2), h[:3]),
            stable_log_series(series_scale(S, z), 2),
        )
    except Exception as exc:
        return {"passed": False, "failure": "nested-series", "detail": repr(exc)}
    scaled = 2 * p[2]
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    if not bool(J[0] > 0 and R[0] > 0 and S[0] > 0 and margin > 0):
        return {
            "passed": False,
            "failure": "dimensionless-margin",
            "mode": [str(left), str(right)],
            "J_ball": J[0].str(35).replace("e", "E"),
            "R_ball": R[0].str(35).replace("e", "E"),
            "S_ball": S[0].str(35).replace("e", "E"),
            "scaled_ball": scaled.str(35).replace("e", "E"),
        }
    return {
        "passed": True,
        "failure": None,
        "mode": [str(left), str(right)],
        "t_ball": t.str(35).replace("e", "E"),
        "J_lower": arb_lower_text(J[0]),
        "R_lower": arb_lower_text(R[0]),
        "S_lower": arb_lower_text(S[0]),
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
            left, right, derivatives, diagnostics=diagnostics
        )
    except Exception as exc:
        result = {"passed": False, "failure": "exception", "detail": repr(exc)}
    return {
        "kind": "order6_nested_curvature_dimensionless_finite_ray_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def load_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
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
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_cache(path, tasks)
    remaining = tasks[len(records) :]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    if workers == 1:
        iterator = map(worker, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers, initializer=initialize_worker
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
                        f"{path.stem}: {len(records)}/{len(tasks)} "
                        f"({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def compact_tail_records() -> list[dict]:
    base_tasks = compact_source.deterministic_tasks()
    base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, base_tasks)
    high = compact_source.load_extension_cache(
        compact_source.DEFAULT_EXTENSION_CACHE, base_tasks
    )
    rows = []
    for low_row, high_row in zip(base, high):
        if Fraction(low_row["mode_left"]) < BASE_TAIL_START:
            continue
        merged = dict(low_row)
        merged["h_derivatives"] = {
            **low_row["h_derivatives"],
            **high_row["h_derivatives"],
        }
        rows.append(merged)
    if not rows or rows[0]["mode_left"] != str(BASE_TAIL_START):
        raise RuntimeError("compact tail records start at the wrong mode")
    return rows


def initial_collar(extension: list[dict]) -> dict:
    flint.ctx.prec = PRECISION_BITS
    records = compact_tail_records() + extension
    derivatives = {
        order: order4_compact.hull(
            [
                order4_compact.interval_from_text(row["h_derivatives"][str(order)])
                for row in records
            ]
        )
        for order in range(2, 11)
    }
    left_t = order4_compact.interval_from_text(records[0]["t_left"])
    right_t = order4_compact.interval_from_text(records[-1]["t_right"])
    central_left_t = potential_jet_arb(arb_rational(MODE_TWO), 1)[1]
    central_right_t = potential_jet_arb(arb_rational(COLLAR_END), 1)[1]
    if not bool(
        left_t < central_left_t - COLLAR_T
        and right_t > central_right_t + COLLAR_T
    ):
        raise RuntimeError("u=2 extension cache does not contain a t+-4 collar")
    result = evaluate_dimensionless(
        MODE_TWO,
        COLLAR_END,
        derivatives,
        diagnostics={
            "record_count": len(records),
            "left_t_collar": (central_left_t - left_t).str(35).replace("e", "E"),
            "right_t_collar": (right_t - central_right_t).str(35).replace("e", "E"),
        },
    )
    if not result.get("passed"):
        raise RuntimeError(f"initial u=2 collar failed: {result}")
    return result


def finite_certificate(extension: list[dict], ray: list[dict]) -> dict:
    collar = initial_collar(extension)
    all_rows = [collar, *ray]
    largest = max(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["scaled_curvature_upper"]),
    )
    weakest_s = min(
        range(len(all_rows)), key=lambda index: Decimal(all_rows[index]["S_lower"])
    )
    selected_indices = sorted(
        {0, len(all_rows) // 4, len(all_rows) // 2, 3 * len(all_rows) // 4, len(all_rows) - 1, largest, weakest_s}
    )
    return {
        "mode_range": ["2", "20"],
        "block_width": str(RAY_WIDTH),
        "initial_collar": collar,
        "ray_blocks": len(ray),
        "all_blocks_passed": True,
        "largest_scaled_index": largest,
        "largest_scaled_curvature_upper": all_rows[largest]["scaled_curvature_upper"],
        "weakest_S_index": weakest_s,
        "weakest_S_lower": all_rows[weakest_s]["S_lower"],
        "selected": [all_rows[index] for index in selected_indices],
        "theorem": "p_1''(t)<=200/t^2 for every saddle mode 2<=u<=20",
    }


def build_artifact(extension: list[dict], ray: list[dict]) -> dict:
    sources = source_contract()
    finite = finite_certificate(extension, ray)
    rows = [
        RayRow(
            "co6ncfr_01_exact_derivatives",
            "exact_theorem_composition",
            "ready_to_apply",
            "Exact cumulant corridors through order ten enclose every required H derivative on the finite ray.",
            "H^(r) enclosed for r=2,...,10 on a t+-4 collar",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co6ncfr_02_initial_collar",
            "interval_certificate",
            "ready_to_apply",
            "The extended compact quadrature cache bridges the mode-u=2 endpoint without using a corridor below its domain.",
            "p_1''(t)<=200/t^2 for 2<=u<=2001/1000",
            "Finite aligned H2-H10 cache only.",
            finite["initial_collar"],
        ),
        RayRow(
            "co6ncfr_03_dimensionless_cover",
            "interval_theorem",
            "ready_to_apply",
            "A rational 1/1000-mode grid and dimensionless stable arithmetic prove the full finite-ray curvature theorem.",
            finite["theorem"],
            "Outward-rounded Arb cover; no point sampling.",
            {
                "blocks": finite["ray_blocks"],
                "largest_scaled_upper": finite["largest_scaled_curvature_upper"],
                "weakest_S_lower": finite["weakest_S_lower"],
            },
        ),
        RayRow(
            "co6ncfr_04_asymptotic_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the dimensionless theorem on the asymptotic mode ray.",
            "p_1''(t)<=200/t^2 for u>=20",
            "Open asymptotic ray only.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate",
        "date": "2026-07-13",
        "status": "rigorous order-six nested first-summand curvature theorem on 2<=u<=20",
        "proof_boundary": (
            "This artifact proves the continuous first-summand curvature only on "
            "the finite saddle ray 2<=u<=20. It does not prove the u>=20 ray, "
            "order-six entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "finite_ray": finite,
        "cache": {
            "extension_path": DEFAULT_EXTENSION_CACHE.relative_to(REPO_ROOT).as_posix(),
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
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_ray"]
    lines = [
        "# Jensen-Window PF Compound Order-Six Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous first-summand order-six curvature theorem on",
        "`2<=u<=20`. This is not a proof of the asymptotic ray, order-six entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py",
        "```",
        "",
        "## Dimensionless Cover",
        "",
        "The exact signed cumulant corridors through order eight and the coarse",
        "absolute corridors for orders nine and ten produce `H^(2),...,H^(10)`",
        "on a strict `t+-4` collar. The stable hierarchy is evaluated after",
        "rescaling by `z=1/t`, preserving the cancellations explicitly.",
        "",
        "A short aligned quadrature extension handles the endpoint `u=2`; the",
        "remaining interval uses a rational mode grid of width `1/1000`.",
        "",
        "```text",
        finite["theorem"],
        f"ray blocks={finite['ray_blocks']},",
        f"largest t^2*p_1'' upper={finite['largest_scaled_curvature_upper']},",
        f"weakest S_1 lower={finite['weakest_S_lower']}.",
        "```",
        "",
        "## Remaining Ray",
        "",
        "```text",
        "Prove p_1''(t)<=200/t^2 for u>=20.",
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md",
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
    ray = build_cache(
        args.ray_cache,
        ray_tasks(),
        ray_task,
        workers=workers,
        overwrite=args.overwrite_caches,
    )
    print(
        f"order-six finite-ray caches: extension={len(extension)}, ray={len(ray)}"
    )
    if args.cache_only:
        return 0
    artifact = build_artifact(extension, ray)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six nested finite-ray certificate: "
        f"{summary['finite_ray_blocks']} ray blocks, "
        f"{summary['finite_ray_theorems']} theorem, "
        f"{summary['open_asymptotic_rays']} open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
