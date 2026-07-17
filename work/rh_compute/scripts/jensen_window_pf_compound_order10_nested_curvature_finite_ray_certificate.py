#!/usr/bin/env python3
"""Certify order-ten first-summand curvature on saddle modes 2<=u<=20."""

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

import jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate as order9_global  # noqa: E402
import jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate as order9  # noqa: E402
from jensen_window_pf_compound_order10_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as HIGH_CUMULANT_SOURCE,
)
from jensen_window_pf_compound_order10_nested_curvature_interval_core import (  # noqa: E402
    evaluate_dimensionless_seventh_curvature,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_COLLAR_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_u2_h17_h18_collar_tiles.jsonl"
)
DEFAULT_RAY_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_u2_u20_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.md"
)
MODE_TWO = order9.MODE_TWO
COLLAR_END = order9.COLLAR_END
BASE_TAIL_START = order9.BASE_TAIL_START
RAY_WIDTH = order9.RAY_WIDTH
COLLAR_T = order9.COLLAR_T
PRECISION_BITS = order9.PRECISION_BITS
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))


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


def collar_tasks() -> list[tuple[int, Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for _, left, right in order9.order9_compact.deterministic_tasks():
        if left >= BASE_TAIL_START:
            intervals.append((left, right))
    for offset in range(order9.order9_compact.RIGHT_COLLAR_TILE_COUNT):
        left = order9.order4.OUTER_END + offset * order9.order4.TILE_WIDTH
        intervals.append((left, left + order9.order4.TILE_WIDTH))
    intervals.extend((left, right) for _, left, right in order9.extension_tasks())
    unique = list(dict.fromkeys(intervals))
    return [(index, left, right) for index, (left, right) in enumerate(unique)]


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order9.ray_tasks()


def collar_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    result = order9.integrate_h_derivatives(
        left,
        right,
        order9.order4.PANELS,
        window_y=order9.order4.WINDOW_Y,
        eighth_envelope_bound=order9.order4.EIGHTH_ENVELOPE,
        max_moment=18,
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
        "kind": "order10_u2_collar_h17_h18_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "h_derivatives": {
            str(degree): order9.order4.interval_text(result["h_derivatives"][degree])
            for degree in (17, 18)
        },
    }


def corridor_h_derivatives(
    left: Fraction,
    right: Fraction,
) -> tuple[dict[int, flint.arb], dict]:
    derivatives, inherited = order9.corridor_h_derivatives(left, right)
    pad = Fraction(inherited["mode_pad"])
    mode = order9.order8_ray.order6_ray.arb_interval(left - pad, right + pad)
    t, curvature = potential_jet_arb(mode, 2)[1:3]
    q_value = flint.arb.pi() * (4 * mode).exp()
    for degree in (17, 18):
        magnitude = (
            math.factorial(degree - 2)
            * q_value ** (1 - flint.arb(degree) / 2)
        )
        cumulant = flint.arb(0, magnitude.upper())
        derivatives[degree] = (
            order9.order8_ray.corridor_core.signed_hurwitz_gamma_derivative(
                degree,
                t + flint.arb(1) / 2,
            )
            - cumulant
            / order9.order8_ray.corridor_core.derivative_power(curvature, degree)
        )
    return derivatives, {
        "mode_pad": inherited["mode_pad"],
        "outer_mode": inherited["outer_mode"],
        "t_collar_product_lower": inherited["t_collar_product_lower"],
        "new_cumulant_caps": {"17": 1, "18": 1},
        "collar_t": COLLAR_T,
    }


def ray_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    try:
        derivatives, diagnostics = corridor_h_derivatives(left, right)
        result = evaluate_dimensionless_seventh_curvature(
            left,
            right,
            derivatives,
            diagnostics=diagnostics,
        )
    except Exception as exc:
        result = {"passed": False, "failure": "exception", "detail": repr(exc)}
    return {
        "kind": "order10_nested_curvature_dimensionless_finite_ray_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def initial_collar(supplemental: list[dict]) -> dict:
    flint.ctx.prec = PRECISION_BITS
    order9_extension = order9.load_cache(
        order9.DEFAULT_EXTENSION_CACHE,
        order9.extension_tasks(),
    )
    records = order9.compact_tail_records() + order9.merged_extension(order9_extension)
    lookup = {
        (record["mode_left"], record["mode_right"]): record
        for record in supplemental
    }
    merged = []
    for record in records:
        key = (record["mode_left"], record["mode_right"])
        if key not in lookup:
            raise RuntimeError(f"missing order-ten collar supplement {key}")
        combined = dict(record)
        combined["h_derivatives"] = {
            **record["h_derivatives"],
            **lookup[key]["h_derivatives"],
        }
        merged.append(combined)
    derivatives = {
        degree: order9.order4.hull(
            [
                order9.order4.interval_from_text(record["h_derivatives"][str(degree)])
                for record in merged
            ]
        )
        for degree in range(2, 19)
    }
    left_t = order9.order4.interval_from_text(merged[0]["t_left"])
    right_t = order9.order4.interval_from_text(merged[-1]["t_right"])
    central_left_t = potential_jet_arb(arb_rational(MODE_TWO), 1)[1]
    central_right_t = potential_jet_arb(arb_rational(COLLAR_END), 1)[1]
    if not bool(
        left_t < central_left_t - COLLAR_T
        and right_t > central_right_t + COLLAR_T
    ):
        raise RuntimeError("order-ten u=2 cache misses the t+-7 collar")
    result = evaluate_dimensionless_seventh_curvature(
        MODE_TWO,
        COLLAR_END,
        derivatives,
        diagnostics={
            "record_count": len(merged),
            "left_t_collar": (central_left_t - left_t).str(35).replace("e", "E"),
            "right_t_collar": (right_t - central_right_t).str(35).replace("e", "E"),
            "new_derivative_orders": [17, 18],
        },
    )
    if not result.get("passed"):
        raise RuntimeError(f"initial order-ten u=2 collar failed: {result}")
    return result


def source_contract() -> dict:
    high = load_json(HIGH_CUMULANT_SOURCE)
    global_order9 = load_json(order9_global.DEFAULT_OUT)
    if high.get("exact", {}).get("exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=17,18, u>=2"
    ):
        raise RuntimeError("order-ten high-cumulant source changed")
    if global_order9.get("summary", {}).get("global_above_5700_compositions") != 1:
        raise RuntimeError("global order-nine curvature source changed")
    return {
        "inherited_order9_finite": order9.source_contract(),
        "order9_global_curvature": (
            "w_1''(t)<=4200/t^2 for every real t>=5700"
        ),
        "new_cumulants": high["exact"]["finite_exact_corridor"],
        "high_cumulant_sha256": sha256(HIGH_CUMULANT_SOURCE),
        "order9_global_sha256": sha256(order9_global.DEFAULT_OUT),
    }


def finite_certificate(ray: list[dict]) -> dict:
    all_rows = ray
    largest = max(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["scaled_curvature_upper"]),
    )
    weakest = min(
        range(len(all_rows)),
        key=lambda index: Decimal(all_rows[index]["W_lower"]),
    )
    selected_indices = sorted(
        {
            0,
            len(all_rows) // 4,
            len(all_rows) // 2,
            3 * len(all_rows) // 4,
            len(all_rows) - 1,
            largest,
            weakest,
        }
    )
    return {
        "mode_range": [str(COLLAR_END), "20"],
        "block_width": str(RAY_WIDTH),
        "ray_blocks": len(ray),
        "all_blocks_passed": True,
        "largest_scaled_index": largest,
        "largest_scaled_curvature_upper": all_rows[largest][
            "scaled_curvature_upper"
        ],
        "weakest_W_index": weakest,
        "weakest_W_lower": all_rows[weakest]["W_lower"],
        "selected": [all_rows[index] for index in selected_indices],
        "theorem": (
            "z_1''(t)<=5500/t^2 for every saddle mode 2001/1000<=u<=20"
        ),
    }


def build_artifact(ray: list[dict]) -> dict:
    finite = finite_certificate(ray)
    rows = [
        RayRow(
            "co10ncfr_01_exact_derivatives",
            "exact_theorem_composition",
            "ready_to_apply",
            "Exact and coarse cumulant corridors through order eighteen enclose every required H derivative.",
            "H^(r) enclosed for r=2,...,18 on t+-7 collars",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co10ncfr_02_dimensionless_cover",
            "interval_theorem",
            "ready_to_apply",
            "A rational 1/1000-mode grid proves the finite-ray curvature theorem.",
            finite["theorem"],
            "Outward-rounded Arb cover; no point sampling.",
            {
                "blocks": finite["ray_blocks"],
                "largest_scaled_upper": finite["largest_scaled_curvature_upper"],
                "weakest_W_lower": finite["weakest_W_lower"],
            },
        ),
        RayRow(
            "co10ncfr_03_compact_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Extend the compact theorem through the exact mode-two collar.",
            "z_1''(t)<=5500/t^2 for 5700<=t<=V'(2001/1000)",
            "Open compact range and collar only.",
        ),
        RayRow(
            "co10ncfr_04_asymptotic_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the dimensionless theorem on the asymptotic ray.",
            "z_1''(t)<=5500/t^2 for u>=20",
            "Open asymptotic ray only.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate",
        "date": "2026-07-16",
        "status": (
            "rigorous order-ten first-summand curvature theorem on "
            "2001/1000<=u<=20"
        ),
        "proof_boundary": (
            "This artifact proves only the displayed finite first-summand ray. "
            "It does not prove the compact range below mode two, the asymptotic "
            "ray, full-kernel transfer, entry, PF-infinity, or RH."
        ),
        "source_contract": source_contract(),
        "finite_ray": finite,
        "cache": {
            "ray_path": DEFAULT_RAY_CACHE.relative_to(REPO_ROOT).as_posix(),
            "ray_rows": len(ray),
            "ray_sha256": sha256(DEFAULT_RAY_CACHE),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": 4,
            "ready_rows": 2,
            "open_rows": 2,
            "initial_collar_blocks": 0,
            "finite_ray_blocks": len(ray),
            "finite_ray_theorems": 1,
            "open_asymptotic_rays": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_ray"]
    lines = [
        "# Jensen-Window PF Compound Order-Ten Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous first-summand order-ten theorem on",
        "`2001/1000<=u<=20`.",
        "This is not a proof of the compact range, asymptotic ray, full-kernel",
        "transfer, entry, PF-infinity, or RH.",
        "",
        "H2-H18 corridors cover strict `t+-7` collars on a `1/1000` mode grid.",
        "",
        "```text",
        finite["theorem"],
        f"ray blocks={finite['ray_blocks']}",
        f"largest t^2*z_1'' upper={finite['largest_scaled_curvature_upper']}",
        f"weakest W lower={finite['weakest_W_lower']}",
        "```",
        "",
        "Open ranges: the compact handoff through `u=2001/1000` and `u>=20`.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ray-cache", type=Path, default=DEFAULT_RAY_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-caches", action="store_true")
    parser.add_argument("--max-ray-blocks", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()
    workers = max(1, args.workers)
    tasks = ray_tasks()
    ray = order9.order7_ray.build_cache(
        args.ray_cache,
        tasks,
        ray_task,
        workers=workers,
        overwrite=args.overwrite_caches,
        max_rows=args.max_ray_blocks,
    )
    print(
        f"order-ten finite-ray cache: ray={len(ray)}/{len(tasks)}"
    )
    if args.cache_only or len(ray) != len(tasks):
        return 0
    artifact = build_artifact(ray)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-ten nested finite-ray certificate: "
        f"{artifact['summary']['finite_ray_blocks']} ray blocks"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
