#!/usr/bin/env python3
"""Certify order-nine first-summand curvature on the ray 2<=u<=20."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate as order7_ray  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate as order8_ray  # noqa: E402
import jensen_window_pf_compound_order9_nested_curvature_compact_certificate as order9_compact  # noqa: E402
import jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_cache as h15_h16  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
    stable_log_series,
)
from jensen_window_pf_compound_order9_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as HIGH_CUMULANT_SOURCE,
    EXACT_CORRIDOR_CAP as NEW_CUMULANT_CAP,
)
from jensen_window_pf_compound_order9_nested_curvature_interval_core import (  # noqa: E402
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
    "jensen_window_pf_compound_order9_nested_curvature_u2_h15_h16_extension_tiles.jsonl"
)
DEFAULT_RAY_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_u2_u20_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.md"
)
MODE_TWO = order8_ray.MODE_TWO
COLLAR_END = order8_ray.COLLAR_END
BASE_TAIL_START = order8_ray.BASE_TAIL_START
RAY_WIDTH = order8_ray.RAY_WIDTH
COLLAR_T = 7
PRECISION_BITS = order8_ray.PRECISION_BITS
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def extension_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order8_ray.extension_tasks()


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order8_ray.ray_tasks()


def load_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
    return order7_ray.load_cache(path, tasks)


def initialize_worker() -> None:
    flint.ctx.prec = PRECISION_BITS


def extension_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4.PANELS,
        window_y=order4.WINDOW_Y,
        eighth_envelope_bound=order4.EIGHTH_ENVELOPE,
        max_moment=16,
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
        "kind": "order9_u2_extension_h15_h16_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "h_derivatives": {
            str(degree): order4.interval_text(result["h_derivatives"][degree])
            for degree in (15, 16)
        },
    }


def choose_pad(left: Fraction, right: Fraction) -> tuple[Fraction, flint.arb]:
    pad = RAY_WIDTH
    while True:
        candidate = pad / 10
        mode = order8_ray.order6_ray.arb_interval(
            left - candidate,
            right + candidate,
        )
        curvature = potential_jet_arb(mode, 2)[2]
        if bool(arb_rational(candidate) * curvature > 2 * COLLAR_T):
            pad = candidate
            continue
        mode = order8_ray.order6_ray.arb_interval(left - pad, right + pad)
        curvature = potential_jet_arb(mode, 2)[2]
        if not bool(arb_rational(pad) * curvature > COLLAR_T):
            raise RuntimeError("adaptive mode pad does not cover t+-7")
        return pad, curvature


def corridor_h_derivatives(
    left: Fraction,
    right: Fraction,
) -> tuple[dict[int, flint.arb], dict]:
    pad, _ = choose_pad(left, right)
    mode = order8_ray.order6_ray.arb_interval(left - pad, right + pad)
    t, curvature = potential_jet_arb(mode, 2)[1:3]
    q_value = flint.arb.pi() * (4 * mode).exp()
    if not bool(
        curvature > 0 and arb_rational(pad) * curvature > COLLAR_T
    ):
        raise RuntimeError("finite-ray mode geometry failed")
    kappa2 = order8_ray.corridor_core.exact_interval(
        1 + arb_rational(order8_ray.order6_ray.CORRIDOR_K2[0]) / q_value,
        1 + arb_rational(order8_ray.order6_ray.CORRIDOR_K2[1]) / q_value,
    )
    derivatives = {
        2: order8_ray.corridor_core.signed_hurwitz_gamma_derivative(
            2,
            t + flint.arb(1) / 2,
        )
        - kappa2 / curvature
    }
    for degree, cap in order8_ray.order6_ray.CORRIDOR_CAPS.items():
        baseline = math.factorial(degree - 2) * q_value ** (
            1 - flint.arb(degree) / 2
        )
        magnitude = order8_ray.corridor_core.exact_interval(
            baseline,
            arb_rational(cap) * baseline,
        )
        cumulant = magnitude if degree % 2 == 0 else -magnitude
        derivatives[degree] = (
            order8_ray.corridor_core.signed_hurwitz_gamma_derivative(
                degree,
                t + flint.arb(1) / 2,
            )
            - cumulant
            / order8_ray.corridor_core.derivative_power(curvature, degree)
        )
    high_caps = (
        (9, order8_ray.MID_CUMULANT_CAP),
        (10, order8_ray.MID_CUMULANT_CAP),
        (11, order8_ray.TOP_CUMULANT_CAP),
        (12, order8_ray.TOP_CUMULANT_CAP),
        (13, order8_ray.ULTRA_CUMULANT_CAP),
        (14, order8_ray.ULTRA_CUMULANT_CAP),
        (15, NEW_CUMULANT_CAP),
        (16, NEW_CUMULANT_CAP),
    )
    for degree, cap in high_caps:
        magnitude = (
            arb_rational(cap)
            * math.factorial(degree - 2)
            * q_value ** (1 - flint.arb(degree) / 2)
        )
        cumulant = flint.arb(0, magnitude.upper())
        derivatives[degree] = (
            order8_ray.corridor_core.signed_hurwitz_gamma_derivative(
                degree,
                t + flint.arb(1) / 2,
            )
            - cumulant
            / order8_ray.corridor_core.derivative_power(curvature, degree)
        )
    return derivatives, {
        "mode_pad": str(pad),
        "outer_mode": [str(left - pad), str(right + pad)],
        "t_collar_product_lower": arb_lower_text(
            arb_rational(pad) * curvature
        ),
        "cumulant_caps": {str(degree): cap for degree, cap in high_caps},
        "collar_t": COLLAR_T,
    }


def evaluate_dimensionless(
    left: Fraction,
    right: Fraction,
    derivatives: dict[int, flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    central_mode = order8_ray.order6_ray.arb_interval(left, right)
    t = potential_jet_arb(central_mode, 1)[1]
    z = 1 / t
    b = [
        t ** (degree + 1)
        * derivatives[degree + 2]
        / math.factorial(degree)
        for degree in range(15)
    ]
    try:
        ell = stable_log_series(series_scale(b, z), 14)
        J = [
            2 * b[degree]
            - z * (degree + 1) * (degree + 2) * ell[degree + 2]
            for degree in range(13)
        ]
        h = series_add(
            series_scale(ell[:13], 2),
            stable_log_series(series_scale(J, z), 12),
        )
        R = [
            3 * b[degree]
            - z * (degree + 1) * (degree + 2) * h[degree + 2]
            for degree in range(11)
        ]
        q = series_add(
            series_sub(series_scale(h[:11], 2), ell[:11]),
            stable_log_series(series_scale(R, z), 10),
        )
        S = [
            4 * b[degree]
            - z * (degree + 1) * (degree + 2) * q[degree + 2]
            for degree in range(9)
        ]
        p = series_add(
            series_sub(series_scale(q[:9], 2), h[:9]),
            stable_log_series(series_scale(S, z), 8),
        )
        T = [
            5 * b[degree]
            - z * (degree + 1) * (degree + 2) * p[degree + 2]
            for degree in range(7)
        ]
        r = series_add(
            series_sub(series_scale(p[:7], 2), q[:7]),
            stable_log_series(series_scale(T, z), 6),
        )
        U = [
            6 * b[degree]
            - z * (degree + 1) * (degree + 2) * r[degree + 2]
            for degree in range(5)
        ]
        s = series_add(
            series_sub(series_scale(r[:5], 2), p[:5]),
            stable_log_series(series_scale(U, z), 4),
        )
        V = [
            7 * b[degree]
            - z * (degree + 1) * (degree + 2) * s[degree + 2]
            for degree in range(3)
        ]
        w = series_add(
            series_sub(series_scale(s[:3], 2), r[:3]),
            stable_log_series(series_scale(V, z), 2),
        )
    except Exception as exc:
        return {"passed": False, "failure": "nested-series", "detail": repr(exc)}
    scaled = 2 * w[2]
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    coordinates = {
        "J": J[0],
        "R": R[0],
        "S": S[0],
        "T": T[0],
        "U": U[0],
        "V": V[0],
    }
    if not bool(all(value > 0 for value in coordinates.values()) and margin > 0):
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
        "kind": "order9_nested_curvature_dimensionless_finite_ray_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def compact_tail_records() -> list[dict]:
    tasks = order9_compact.deterministic_tasks()
    caches = [
        order9_compact.order4.load_cache(order9_compact.order4.DEFAULT_CACHE, tasks),
        order9_compact.order6.load_extension_cache(
            order9_compact.order6.DEFAULT_EXTENSION_CACHE, tasks
        ),
        order9_compact.order7.load_extension_cache(
            order9_compact.order7.DEFAULT_EXTENSION_CACHE, tasks
        ),
        order9_compact.order8.load_cache(order9_compact.order8.DEFAULT_CACHE, tasks),
        h15_h16.load_cache(h15_h16.DEFAULT_CACHE, tasks),
    ]
    collar = order9_compact.build_right_collar(order9_compact.DEFAULT_RIGHT_COLLAR)
    caches = order9_compact.append_right_collar(caches, collar)
    rows = []
    for aligned in zip(*caches):
        if Fraction(aligned[0]["mode_left"]) < BASE_TAIL_START:
            continue
        merged = dict(aligned[0])
        merged["h_derivatives"] = {}
        for record in aligned:
            merged["h_derivatives"].update(record["h_derivatives"])
        rows.append(merged)
    if not rows or rows[0]["mode_left"] != str(BASE_TAIL_START):
        raise RuntimeError("order-nine compact tail starts at the wrong mode")
    return rows


def merged_extension(new: list[dict]) -> list[dict]:
    tasks = extension_tasks()
    inherited_ultra = order8_ray.load_cache(order8_ray.DEFAULT_EXTENSION_CACHE, tasks)
    inherited = order8_ray.merged_extension(inherited_ultra)
    if len(inherited) != len(tasks) or len(new) != len(tasks):
        raise RuntimeError("order-nine u=2 extension caches are incomplete")
    rows = []
    for base, extension in zip(inherited, new):
        if (base["mode_left"], base["mode_right"]) != (
            extension["mode_left"],
            extension["mode_right"],
        ):
            raise RuntimeError("order-nine u=2 extension is not aligned")
        merged = dict(base)
        merged["h_derivatives"] = {
            **base["h_derivatives"],
            **extension["h_derivatives"],
        }
        rows.append(merged)
    return rows


def initial_collar(extension: list[dict]) -> dict:
    flint.ctx.prec = PRECISION_BITS
    records = compact_tail_records() + merged_extension(extension)
    derivatives = {
        degree: order4.hull(
            [
                order4.interval_from_text(record["h_derivatives"][str(degree)])
                for record in records
            ]
        )
        for degree in range(2, 17)
    }
    left_t = order4.interval_from_text(records[0]["t_left"])
    right_t = order4.interval_from_text(records[-1]["t_right"])
    central_left_t = potential_jet_arb(arb_rational(MODE_TWO), 1)[1]
    central_right_t = potential_jet_arb(arb_rational(COLLAR_END), 1)[1]
    if not bool(
        left_t < central_left_t - COLLAR_T
        and right_t > central_right_t + COLLAR_T
    ):
        raise RuntimeError("u=2 extension cache misses the t+-7 collar")
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
        raise RuntimeError(f"initial order-nine u=2 collar failed: {result}")
    return result


def source_contract() -> dict:
    compact = load_json(order9_compact.DEFAULT_OUT)
    high = load_json(HIGH_CUMULANT_SOURCE)
    if compact.get("compact", {}).get("all_blocks_passed") is not True:
        raise RuntimeError("order-nine compact handoff changed")
    if high.get("exact", {}).get("exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=15,16, u>=2"
    ):
        raise RuntimeError("order-nine high-cumulant source changed")
    return {
        "inherited_order8": order8_ray.source_contract(),
        "new_cumulants": high["exact"]["finite_exact_corridor"],
        "compact_handoff": compact["compact"]["theorem"],
        "high_cumulant_sha256": sha256(HIGH_CUMULANT_SOURCE),
        "compact_sha256": sha256(order9_compact.DEFAULT_OUT),
    }


def finite_certificate(extension: list[dict], ray: list[dict]) -> dict:
    collar = initial_collar(extension)
    all_rows = [collar, *ray]
    largest = max(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["scaled_curvature_upper"]),
    )
    weakest = min(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["V_lower"]),
    )
    selected_indices = sorted(
        {0, len(all_rows) // 4, len(all_rows) // 2, 3 * len(all_rows) // 4,
         len(all_rows) - 1, largest, weakest}
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
        "weakest_V_index": weakest,
        "weakest_V_lower": all_rows[weakest]["V_lower"],
        "selected": [all_rows[index] for index in selected_indices],
        "theorem": "w_1''(t)<=4200/t^2 for every saddle mode 2<=u<=20",
    }


def build_artifact(extension: list[dict], ray: list[dict]) -> dict:
    finite = finite_certificate(extension, ray)
    rows = [
        RayRow(
            "co9ncfr_01_exact_derivatives",
            "exact_theorem_composition",
            "ready_to_apply",
            "Exact and coarse cumulant corridors through order sixteen enclose every required H derivative.",
            "H^(r) enclosed for r=2,...,16 on t+-7 collars",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co9ncfr_02_initial_collar",
            "interval_certificate",
            "ready_to_apply",
            "Aligned compact and extension caches bridge mode u=2.",
            "w_1''(t)<=4200/t^2 for 2<=u<=2001/1000",
            "Finite aligned H2-H16 cache only.",
            finite["initial_collar"],
        ),
        RayRow(
            "co9ncfr_03_dimensionless_cover",
            "interval_theorem",
            "ready_to_apply",
            "A rational 1/1000-mode grid proves the finite-ray curvature theorem.",
            finite["theorem"],
            "Outward-rounded Arb cover; no point sampling.",
            {
                "blocks": finite["ray_blocks"],
                "largest_scaled_upper": finite["largest_scaled_curvature_upper"],
                "weakest_V_lower": finite["weakest_V_lower"],
            },
        ),
        RayRow(
            "co9ncfr_04_asymptotic_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the dimensionless theorem on the asymptotic ray.",
            "w_1''(t)<=4200/t^2 for u>=20",
            "Open asymptotic ray only.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate",
        "date": "2026-07-14",
        "status": "rigorous order-nine first-summand curvature theorem on 2<=u<=20",
        "proof_boundary": (
            "This artifact proves only the displayed finite first-summand ray. "
            "It does not prove the asymptotic ray, full-kernel transfer, entry, "
            "PF-infinity, or RH."
        ),
        "source_contract": source_contract(),
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
            "rows": 4,
            "ready_rows": 3,
            "open_rows": 1,
            "initial_collar_blocks": 1,
            "finite_ray_blocks": len(ray),
            "finite_ray_theorems": 1,
            "open_asymptotic_rays": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_ray"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous first-summand order-nine theorem on `2<=u<=20`.",
        "This is not a proof of the asymptotic ray, entry, PF-infinity, or RH.",
        "",
        "H2-H16 corridors cover strict `t+-7` collars on a `1/1000` mode grid.",
        "",
        "```text",
        finite["theorem"],
        f"ray blocks={finite['ray_blocks']}",
        f"largest t^2*w_1'' upper={finite['largest_scaled_curvature_upper']}",
        f"weakest V lower={finite['weakest_V_lower']}",
        "```",
        "",
        "Remaining ray: prove `w_1''(t)<=4200/t^2` for `u>=20`.",
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
    extension = order7_ray.build_cache(
        args.extension_cache,
        extension_tasks(),
        extension_task,
        workers=workers,
        overwrite=args.overwrite_caches,
    )
    tasks = ray_tasks()
    ray = order7_ray.build_cache(
        args.ray_cache,
        tasks,
        ray_task,
        workers=workers,
        overwrite=args.overwrite_caches,
        max_rows=args.max_ray_blocks,
    )
    print(
        f"order-nine finite-ray caches: extension={len(extension)}, "
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
    print(
        "wrote order-nine nested finite-ray certificate: "
        f"{artifact['summary']['finite_ray_blocks']} ray blocks"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
