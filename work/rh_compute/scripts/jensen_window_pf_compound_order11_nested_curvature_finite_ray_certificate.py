#!/usr/bin/env python3
"""Certify order-eleven first-summand curvature on 2.001<=u<=20."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate as order9  # noqa: E402
from jensen_window_pf_compound_order11_nested_curvature_interval_core import (  # noqa: E402
    evaluate_dimensionless_eighth_curvature,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
)


RESULTS = REPO_ROOT / "work/rh_compute/results"
DEFAULT_CACHE = RESULTS / "jensen_window_pf_compound_order11_nested_curvature_u2_u20_blocks.jsonl"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate.md"
SOURCE_HIGH = RESULTS / "jensen_window_pf_compound_order11_high_cumulant_coarse_corridor.json"
SOURCE_ORDER10_FINITE = RESULTS / "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.json"
SOURCE_ORDER10_GLOBAL = RESULTS / "jensen_window_pf_compound_order10_first_summand_curvature_certificate.json"

MODE_START = Fraction(2001, 1000)
MODE_END = Fraction(20)
RAY_WIDTH = Fraction(1, 1000)
COLLAR_T = 8
PRECISION_BITS = order9.PRECISION_BITS
DEFAULT_WORKERS = max(1, min(4, (os.cpu_count() or 4) - 1))
THEOREM = "y_1''(t)<=6000/t^2 for every saddle mode 2001/1000<=u<=20"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_contract() -> dict:
    high = load_json(SOURCE_HIGH)
    inherited = load_json(SOURCE_ORDER10_FINITE)
    global_order10 = load_json(SOURCE_ORDER10_GLOBAL)
    if high.get("exact", {}).get("exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=19,20, u>=2"
    ):
        raise RuntimeError("order-eleven high-cumulant source changed")
    if (
        inherited.get("status")
        != "rigorous order-ten first-summand curvature theorem on 2001/1000<=u<=20"
        or inherited.get("finite_ray", {}).get("all_blocks_passed") is not True
    ):
        raise RuntimeError("inherited finite derivative-corridor source changed")
    if global_order10.get("theorem") != "z_1''(t)<=4200/t^2 for every real t>=1251":
        raise RuntimeError("global order-ten curvature source changed")
    return {
        "inherited_H2_H18": inherited["finite_ray"]["theorem"],
        "new_H19_H20": high["exact"]["exact_corridor"],
        "order10_global_curvature": global_order10["theorem"],
        "high_cumulant_sha256": sha256(SOURCE_HIGH),
        "order10_finite_sha256": sha256(SOURCE_ORDER10_FINITE),
        "order10_global_sha256": sha256(SOURCE_ORDER10_GLOBAL),
    }


def ray_tasks() -> list[tuple[int, Fraction, Fraction]]:
    count = int((MODE_END - MODE_START) / RAY_WIDTH)
    return [
        (index, MODE_START + index * RAY_WIDTH, MODE_START + (index + 1) * RAY_WIDTH)
        for index in range(count)
    ]


def choose_pad(left: Fraction, right: Fraction) -> tuple[Fraction, flint.arb]:
    pad = RAY_WIDTH
    while True:
        candidate = pad / 10
        mode = order9.order8_ray.order6_ray.arb_interval(left - candidate, right + candidate)
        curvature = potential_jet_arb(mode, 2)[2]
        if bool(arb_rational(candidate) * curvature > 2 * COLLAR_T):
            pad = candidate
            continue
        mode = order9.order8_ray.order6_ray.arb_interval(left - pad, right + pad)
        curvature = potential_jet_arb(mode, 2)[2]
        if not bool(arb_rational(pad) * curvature > COLLAR_T):
            raise RuntimeError("adaptive mode pad does not cover t+-8")
        return pad, curvature


def corridor_h_derivatives(
    left: Fraction,
    right: Fraction,
    *,
    extra_absolute_caps: tuple[tuple[int, int], ...] = (),
) -> tuple[dict[int, flint.arb], dict]:
    pad, _ = choose_pad(left, right)
    mode = order9.order8_ray.order6_ray.arb_interval(left - pad, right + pad)
    t, curvature = potential_jet_arb(mode, 2)[1:3]
    q_value = flint.arb.pi() * (4 * mode).exp()
    if not bool(curvature > 0 and arb_rational(pad) * curvature > COLLAR_T):
        raise RuntimeError("finite-ray mode geometry failed")
    corridor_core = order9.order8_ray.corridor_core
    order6_ray = order9.order8_ray.order6_ray
    kappa2 = corridor_core.exact_interval(
        1 + arb_rational(order6_ray.CORRIDOR_K2[0]) / q_value,
        1 + arb_rational(order6_ray.CORRIDOR_K2[1]) / q_value,
    )
    derivatives = {
        2: corridor_core.signed_hurwitz_gamma_derivative(2, t + flint.arb(1) / 2)
        - kappa2 / curvature
    }
    signed_caps = list(order6_ray.CORRIDOR_CAPS.items())
    absolute_caps = [
        (9, order9.order8_ray.MID_CUMULANT_CAP),
        (10, order9.order8_ray.MID_CUMULANT_CAP),
        (11, order9.order8_ray.TOP_CUMULANT_CAP),
        (12, order9.order8_ray.TOP_CUMULANT_CAP),
        (13, order9.order8_ray.ULTRA_CUMULANT_CAP),
        (14, order9.order8_ray.ULTRA_CUMULANT_CAP),
        (15, order9.NEW_CUMULANT_CAP),
        (16, order9.NEW_CUMULANT_CAP),
        (17, 1),
        (18, 1),
        (19, 1),
        (20, 1),
    ]
    if any(degree <= 20 or cap <= 0 for degree, cap in extra_absolute_caps):
        raise ValueError("extra absolute cumulant caps must have degree>20 and cap>0")
    absolute_caps.extend(extra_absolute_caps)
    for degree, cap in signed_caps:
        baseline = math.factorial(degree - 2) * q_value ** (1 - flint.arb(degree) / 2)
        magnitude = corridor_core.exact_interval(baseline, arb_rational(cap) * baseline)
        cumulant = magnitude if degree % 2 == 0 else -magnitude
        derivatives[degree] = (
            corridor_core.signed_hurwitz_gamma_derivative(degree, t + flint.arb(1) / 2)
            - cumulant / corridor_core.derivative_power(curvature, degree)
        )
    for degree, cap in absolute_caps:
        magnitude = (
            arb_rational(cap)
            * math.factorial(degree - 2)
            * q_value ** (1 - flint.arb(degree) / 2)
        )
        cumulant = flint.arb(0, magnitude.upper())
        derivatives[degree] = (
            corridor_core.signed_hurwitz_gamma_derivative(degree, t + flint.arb(1) / 2)
            - cumulant / corridor_core.derivative_power(curvature, degree)
        )
    return derivatives, {
        "mode_pad": str(pad),
        "outer_mode": [str(left - pad), str(right + pad)],
        "t_collar_product_lower": arb_lower_text(arb_rational(pad) * curvature),
        "new_cumulant_caps": {
            str(degree): cap for degree, cap in absolute_caps if degree >= 19
        },
        "collar_t": COLLAR_T,
    }


def ray_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = PRECISION_BITS
    try:
        derivatives, diagnostics = corridor_h_derivatives(left, right)
        result = evaluate_dimensionless_eighth_curvature(
            left,
            right,
            derivatives,
            diagnostics=diagnostics,
        )
    except Exception as exc:
        result = {"passed": False, "failure": "exception", "detail": repr(exc)}
    return {
        "kind": "order11_nested_curvature_dimensionless_finite_ray_block",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        **result,
    }


def load_cache(path: Path, tasks: list[tuple[int, Fraction, Fraction]]) -> list[dict]:
    if not path.exists():
        return []
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            index = len(records)
            if index >= len(tasks):
                raise RuntimeError("finite-ray cache has extra rows")
            expected = tasks[index]
            if (
                record.get("kind") != "order11_nested_curvature_dimensionless_finite_ray_block"
                or record.get("index") != expected[0]
                or record.get("mode_left") != str(expected[1])
                or record.get("mode_right") != str(expected[2])
                or record.get("passed") is not True
            ):
                raise RuntimeError(f"invalid finite-ray cache row {line_number}")
            records.append(record)
    return records


def build_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    runtime_seconds: float,
) -> list[dict]:
    records = load_cache(path, tasks)
    path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with path.open("a", encoding="utf-8") as handle, ProcessPoolExecutor(max_workers=workers) as pool:
        while len(records) < len(tasks) and time.monotonic() - started < runtime_seconds:
            batch = tasks[len(records) : len(records) + workers]
            results = list(pool.map(ray_task, batch))
            for record in results:
                if record.get("passed") is not True:
                    raise RuntimeError(f"finite-ray block failed: {record}")
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                handle.flush()
                records.append(record)
            if len(records) % 1000 == 0 or len(records) == len(tasks):
                print(f"order-eleven finite ray: {len(records)}/{len(tasks)}")
    return records


def build_artifact(records: list[dict]) -> dict:
    tasks = ray_tasks()
    if len(records) != len(tasks):
        raise RuntimeError(f"incomplete finite-ray cache: {len(records)}/{len(tasks)}")
    largest = max(records, key=lambda row: Fraction(row["scaled_curvature_upper"]))
    weakest = min(records, key=lambda row: Fraction(row["scaled_margin_lower"]))
    weakest_x = min(records, key=lambda row: Fraction(row["X_lower"]))
    return {
        "kind": "jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate",
        "date": "2026-07-18",
        "status": "rigorous order-eleven first-summand curvature theorem on 2001/1000<=u<=20",
        "theorem": THEOREM,
        "proof_boundary": (
            "This artifact proves only the finite first-summand saddle ray. The compact "
            "handoff, global composition, full-kernel transfer, PF-infinity, and RH remain separate."
        ),
        "source_contract": source_contract(),
        "finite_ray": {
            "mode_range": [str(MODE_START), str(MODE_END)],
            "block_width": str(RAY_WIDTH),
            "blocks": len(records),
            "all_blocks_passed": True,
            "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
            "largest_scaled_curvature_block": largest["index"],
            "smallest_margin_lower": weakest["scaled_margin_lower"],
            "weakest_X_lower": weakest_x["X_lower"],
            "cache_sha256": sha256(DEFAULT_CACHE),
        },
        "summary": {
            "finite_ray_theorems": 1,
            "blocks": len(records),
            "failed_blocks": 0,
            "open_compact_ranges": 1,
            "global_first_summand_theorems": 0,
            "rh_claims": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite_ray"]
    lines = [
        "# Order-Eleven Nested Curvature Finite-Ray Certificate",
        "",
        "Date: 2026-07-18",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH.",
        "",
        "```text",
        artifact["theorem"],
        f"blocks={finite['blocks']}",
        f"largest scaled upper={finite['largest_scaled_curvature_upper']}",
        f"smallest margin={finite['smallest_margin_lower']}",
        f"weakest X lower={finite['weakest_X_lower']}",
        "```",
        "",
        artifact["proof_boundary"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--runtime-seconds", type=float, default=7200)
    args = parser.parse_args()
    source_contract()
    tasks = ray_tasks()
    records = build_cache(
        args.cache,
        tasks,
        workers=max(1, args.workers),
        runtime_seconds=max(1.0, args.runtime_seconds),
    )
    if len(records) != len(tasks):
        print(f"order-eleven finite-ray checkpoint: {len(records)}/{len(tasks)}")
        return 0
    artifact = build_artifact(records)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "wrote order-eleven finite-ray theorem: "
        f"{len(records)} blocks, largest scaled upper "
        f"{artifact['finite_ray']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
