#!/usr/bin/env python3
"""Test named exact-point profiles on selected order-ten compact blocks."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


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
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    localized_seventh_formula_continuation_row,
)


DEFAULT_ANCHORS = (
    Fraction(5700),
    Fraction(10000),
    Fraction(20000),
    Fraction(30000),
    Fraction(38018),
)


def selected_point_tasks(
    anchors: list[Fraction],
    profile_name: str,
) -> list[tuple[int, Fraction, str]]:
    targets = sorted(
        {
            anchor + shift
            for anchor in anchors
            for shift in range(-8, 9)
        }
    )
    if any(target.denominator != 1 for target in targets):
        raise ValueError("profile scout anchors must be integers")
    return [
        (index, target, profile_name)
        for index, target in enumerate(targets)
    ]


def load_selected_h_rows(
    path: Path,
    anchors: list[Fraction],
    block_width: Fraction,
) -> dict[Fraction, dict]:
    support = {
        (anchor - 8, anchor + 8 + block_width)
        for anchor in anchors
    }
    selected = {}
    total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            total += 1
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            right = Fraction(record["target_t_right"])
            if not any(left >= start and right <= end for start, end in support):
                continue
            derivatives = record.get("h_derivatives", {})
            if (
                record.get("kind") != "order10_compact_h2_h24_unit_tile"
                or record.get("contract_id") != h_cache.ROW_CONTRACT
                or record.get("passed") is not True
                or set(derivatives)
                != {str(order) for order in range(2, h_cache.MAX_MOMENT + 1)}
            ):
                raise RuntimeError(f"invalid selected H row at {left}")
            selected[left] = {
                "target_t_left": left,
                "target_t_right": right,
                "H": {
                    order: compact.interval_from_text(derivatives[str(order)])
                    for order in range(2, h_cache.MAX_MOMENT + 1)
                },
            }
    expected_total = len(
        h_cache.deterministic_tasks(
            h_cache.DEFAULT_START_T,
            h_cache.DEFAULT_END_T,
            h_cache.DEFAULT_TILE_WIDTH_T,
        )
    )
    if total != expected_total:
        raise RuntimeError(
            f"compact H cache is incomplete: {total}/{expected_total} rows"
        )
    return selected


def point_source(records: list[dict]) -> dict[Fraction, tuple[list, dict]]:
    result = {}
    for record in records:
        target = Fraction(record["target_t"])
        derivatives = record["h_derivatives"]
        result[target] = (
            [
                compact.interval_from_text(derivatives[str(order)])
                / math.factorial(order)
                for order in range(point_cache.MAX_MOMENT + 1)
            ],
            {
                "target_t": str(target),
                "mode_bracket": [record["mode_left"], record["mode_right"]],
                "maximum_panel_error_upper": record[
                    "maximum_panel_error_upper"
                ],
                "maximum_tail_moment_upper": record[
                    "maximum_tail_moment_upper"
                ],
                "minimum_tail_slope_lower": record[
                    "minimum_tail_slope_lower"
                ],
            },
        )
    return result


def block_h_rows(
    selected: dict[Fraction, dict],
    anchor: Fraction,
    block_width: Fraction,
) -> list[dict]:
    rows = [
        selected[target]
        for target in range(int(anchor - 8), int(anchor + 8 + block_width))
    ]
    if (
        len(rows) != int(16 + block_width)
        or rows[0]["target_t_left"] != anchor - 8
        or rows[-1]["target_t_right"] != anchor + 8 + block_width
    ):
        raise RuntimeError(f"incomplete H collar at anchor {anchor}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--anchor",
        action="append",
        type=Fraction,
        dest="anchors",
    )
    parser.add_argument(
        "--profile",
        choices=tuple(point_cache.PROFILE_SPECS),
        required=True,
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--width", type=Fraction, default=Fraction(2))
    parser.add_argument("--point-cache", type=Path, required=True)
    parser.add_argument("--h-cache", type=Path, default=h_cache.DEFAULT_CACHE)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    anchors = sorted(set(args.anchors or DEFAULT_ANCHORS))
    if args.width.denominator != 1 or not 0 < args.width <= 2:
        raise ValueError("scout width must be one or two")
    if any(anchor < 5700 or anchor > 38018 for anchor in anchors):
        raise ValueError("scout anchors must lie in 5700..38018")
    tasks = selected_point_tasks(anchors, args.profile)
    records = point_cache.build_cache(
        args.point_cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_points=None,
    )
    if len(records) != len(tasks):
        raise RuntimeError("selected exact-point cache is incomplete")

    flint.ctx.prec = PRECISION_BITS
    selected_h = load_selected_h_rows(args.h_cache, anchors, args.width)
    points = point_source(records)
    blocks = []
    for anchor in anchors:
        row = localized_seventh_formula_continuation_row(
            anchor,
            anchor + args.width,
            block_h_rows(selected_h, anchor, args.width),
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
    output = {
        "profile": args.profile,
        "profile_spec": point_cache.PROFILE_SPECS[args.profile],
        "point_rows": len(records),
        "blocks": blocks,
        "all_blocks_passed": all(row["passed"] for row in blocks),
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0 if output["all_blocks_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
