#!/usr/bin/env python3
"""Scout sparse H0-H23 point jets with rigorous order-24 propagation."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
import math
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_cache  # noqa: E402
import jensen_window_pf_compound_order10_compact_point_h0_h7_cache as point_cache  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    hull,
    integrate_h_jet_taylor_quadrature,
)
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    _symmetric,
    localized_seventh_formula_continuation_row,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
    upper_absolute,
)


SPARSE_ORIGIN = Fraction(5692)
DEFAULT_SPARSE_STEP = 8
EXACT_MAX_MOMENT = 23
REMAINDER_MOMENT = 24


def sparse_base(target: Fraction, sparse_step: int = DEFAULT_SPARSE_STEP) -> Fraction:
    quotient = (target - SPARSE_ORIGIN) // sparse_step
    return SPARSE_ORIGIN + quotient * sparse_step


def selected_tasks(
    anchors: list[Fraction],
    profile_name: str,
    sparse_step: int,
) -> list[tuple[int, Fraction, str]]:
    targets = {
        anchor + shift
        for anchor in anchors
        for shift in range(-8, 9)
    }
    bases = sorted({sparse_base(target, sparse_step) for target in targets})
    return [
        (index, target, profile_name)
        for index, target in enumerate(bases)
    ]


def initialize_worker(profile_name: str) -> None:
    flint.ctx.prec = point_cache.PROFILE_SPECS[profile_name]["precision_bits"]


def exact_task(task: tuple[int, Fraction, str]) -> dict:
    index, target, profile_name = task
    profile = point_cache.PROFILE_SPECS[profile_name]
    flint.ctx.prec = profile["precision_bits"]
    mode_left, mode_right = point_cache.mode_bracket(target, profile_name)
    result = integrate_h_jet_taylor_quadrature(
        mode_left,
        mode_right,
        profile["panels"],
        window_y=point_cache.WINDOW_Y,
        taylor_order=point_cache.TAYLOR_ORDER,
        max_moment=EXACT_MAX_MOMENT,
    )
    if not result.get("passed"):
        return {
            "index": index,
            "target_t": str(target),
            "profile": profile_name,
            "passed": False,
            "failure": result.get("failure"),
        }
    mode_ball = (
        (arb_rational(mode_left) + arb_rational(mode_right)) / 2
        + flint.arb(
            0,
            (arb_rational(mode_right) - arb_rational(mode_left)) / 2,
        )
    )
    target_ball = potential_jet_arb(mode_ball, 1)[1]
    if not bool(target_ball.contains(arb_rational(target))):
        return {
            "index": index,
            "target_t": str(target),
            "profile": profile_name,
            "passed": False,
            "failure": "target-not-in-mode-bracket",
        }
    return {
        "kind": "order10_sparse_exact_point_h0_h23_jet",
        "index": index,
        "target_t": str(target),
        "profile": profile_name,
        "mode_left": str(mode_left),
        "mode_right": str(mode_right),
        "target_t_ball": target_ball.str(50).replace("e", "E"),
        "passed": True,
        "h_derivatives": {
            str(order): compact.interval_text(result["h_derivatives"][order])
            for order in range(EXACT_MAX_MOMENT + 1)
        },
        "maximum_panel_error_upper": result["maximum_panel_error_upper"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
    }


def build_exact_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, str]],
    *,
    workers: int,
    overwrite: bool,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    if path.exists():
        records = [
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    else:
        records = []
    for record, task in zip(records, tasks):
        index, target, profile_name = task
        if (
            record.get("kind") != "order10_sparse_exact_point_h0_h23_jet"
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("profile") != profile_name
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"sparse exact cache mismatch at row {index}")
    if len(records) > len(tasks):
        raise RuntimeError("sparse exact cache has too many rows")
    remaining = tasks[len(records) :]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    profile_name = tasks[0][2]
    started = perf_counter()
    if workers == 1:
        initialize_worker(profile_name)
        iterator = map(exact_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=(profile_name,),
        )
        iterator = executor.map(exact_task, remaining, chunksize=1)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for record in iterator:
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"sparse exact row failed: {json.dumps(record)}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                handle.flush()
                records.append(record)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    print(
        f"sparse exact H0-H23 rows: {len(records)}/{len(tasks)} "
        f"({perf_counter() - started:.1f}s)"
    )
    return records


def load_h_map(path: Path) -> dict[Fraction, dict]:
    result = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            derivatives = record["h_derivatives"]
            result[left] = {
                "target_t_left": left,
                "target_t_right": Fraction(record["target_t_right"]),
                "H": {
                    order: compact.interval_from_text(derivatives[str(order)])
                    for order in range(2, 25)
                },
            }
    if len(result) != 32336:
        raise RuntimeError("compact H cache is incomplete")
    return result


def load_exact_coefficients(records: list[dict]) -> tuple[dict, dict]:
    coefficients = {}
    diagnostics = {}
    for record in records:
        target = Fraction(record["target_t"])
        derivatives = record["h_derivatives"]
        coefficients[target] = [
            compact.interval_from_text(derivatives[str(order)])
            / math.factorial(order)
            for order in range(EXACT_MAX_MOMENT + 1)
        ]
        diagnostics[target] = {
            "target_t": str(target),
            "mode_bracket": [record["mode_left"], record["mode_right"]],
            "maximum_panel_error_upper": record["maximum_panel_error_upper"],
            "maximum_tail_moment_upper": record["maximum_tail_moment_upper"],
            "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
        }
    return coefficients, diagnostics


def propagated_point_source(
    anchors: list[Fraction],
    exact_records: list[dict],
    h_map: dict[Fraction, dict],
    sparse_step: int,
) -> tuple[dict, dict]:
    exact, exact_diagnostics = load_exact_coefficients(exact_records)
    targets = sorted(
        {
            anchor + shift
            for anchor in anchors
            for shift in range(-8, 9)
        }
    )
    result = {}
    remainder_diagnostics = {}
    for target in targets:
        base = sparse_base(target, sparse_step)
        displacement = int(target - base)
        base_coefficients = exact[base]
        propagated = []
        per_order = {}
        if displacement:
            h24 = hull(
                [
                    h_map[base + offset]["H"][REMAINDER_MOMENT]
                    for offset in range(displacement)
                ]
            )
            h24_abs = upper_absolute(h24)
        else:
            h24_abs = flint.arb(0)
        for derivative in range(8):
            value = flint.arb(0)
            for moment in range(EXACT_MAX_MOMENT - derivative, -1, -1):
                value = value * displacement + (
                    math.comb(derivative + moment, derivative)
                    * base_coefficients[derivative + moment]
                )
            remainder = (
                h24_abs
                * displacement ** (REMAINDER_MOMENT - derivative)
                / math.factorial(REMAINDER_MOMENT - derivative)
                / math.factorial(derivative)
            )
            propagated.append(value + _symmetric(remainder))
            per_order[str(derivative)] = compact.interval_text(remainder)
        result[target] = (
            propagated,
            {
                **exact_diagnostics[base],
                "target_t": str(target),
                "sparse_base_t": str(base),
                "sparse_displacement": displacement,
            },
        )
        remainder_diagnostics[str(target)] = per_order
    return result, remainder_diagnostics


def block_h_rows(
    h_map: dict[Fraction, dict],
    anchor: Fraction,
    width: Fraction,
) -> list[dict]:
    return [
        h_map[Fraction(target)]
        for target in range(int(anchor - 8), int(anchor + width + 8))
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor", action="append", type=Fraction, required=True)
    parser.add_argument("--width", type=Fraction, default=Fraction(2))
    parser.add_argument(
        "--profile",
        choices=tuple(point_cache.PROFILE_SPECS),
        required=True,
    )
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--sparse-step", type=int, default=DEFAULT_SPARSE_STEP)
    parser.add_argument("--exact-cache", type=Path, required=True)
    parser.add_argument("--h-cache", type=Path, default=h_cache.DEFAULT_CACHE)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    anchors = sorted(set(args.anchor))
    if args.width not in (Fraction(1), Fraction(2)):
        raise ValueError("scout width must be one or two")
    if not 1 <= args.sparse_step <= 32:
        raise ValueError("sparse step must lie in 1..32")
    tasks = selected_tasks(anchors, args.profile, args.sparse_step)
    records = build_exact_cache(
        args.exact_cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
    )
    flint.ctx.prec = PRECISION_BITS
    h_map = load_h_map(args.h_cache)
    points, remainders = propagated_point_source(
        anchors,
        records,
        h_map,
        args.sparse_step,
    )
    blocks = []
    for anchor in anchors:
        row = localized_seventh_formula_continuation_row(
            anchor,
            anchor + args.width,
            block_h_rows(h_map, anchor, args.width),
            point_order=7,
            remainder_order=8,
            point_h_source=points,
            require_pass=False,
        )
        blocks.append(
            {
                "anchor": str(anchor),
                "right": str(anchor + args.width),
                "profile": args.profile,
                "scaled_curvature_upper": row["scaled_curvature_upper"],
                "curvature_margin_lower": row["curvature_margin_lower"],
                "W_lower": row["W_lower"],
                "passed": row["passed"],
            }
        )
    print(
        json.dumps(
            {
                "profile": args.profile,
                "sparse_step": args.sparse_step,
                "exact_max_moment": EXACT_MAX_MOMENT,
                "exact_rows": len(records),
                "blocks": blocks,
                "all_blocks_passed": all(row["passed"] for row in blocks),
                "maximum_propagation_remainder_upper": {
                    str(order): max(
                        (
                            float(values[str(order)][1])
                            for values in remainders.values()
                        ),
                        default=0.0,
                    )
                    for order in range(8)
                },
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if all(row["passed"] for row in blocks) else 2


if __name__ == "__main__":
    raise SystemExit(main())
