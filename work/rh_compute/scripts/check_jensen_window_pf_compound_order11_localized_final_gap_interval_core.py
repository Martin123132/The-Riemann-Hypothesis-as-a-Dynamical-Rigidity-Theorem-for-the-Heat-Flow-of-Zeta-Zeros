#!/usr/bin/env python3
"""Exercise the order-eleven localized final-gap interval core."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order11_localized_final_gap_interval_core as core  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_compact_certificate import (  # noqa: E402
    interval_from_text,
)


POINT_CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
POINT_ANCHOR = Fraction(2000)
INTERVAL_ANCHOR = Fraction(5701)


def load_point_source() -> dict[Fraction, tuple[list, dict]]:
    needed = {
        POINT_ANCHOR + shift
        for shift in range(-9, 10)
    }
    result = {}
    with POINT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            target = Fraction(row["target_t"])
            if target not in needed:
                continue
            derivatives = row["h_derivatives"]
            result[target] = (
                [
                    interval_from_text(derivatives[str(degree)]) / math.factorial(degree)
                    for degree in range(core.POINT_ORDER + 1)
                ],
                {
                    "target_t": str(target),
                    "mode_bracket": [row["mode_left"], row["mode_right"]],
                    "maximum_panel_error_upper": row["maximum_panel_error_upper"],
                    "maximum_tail_moment_upper": row["maximum_tail_moment_upper"],
                },
            )
    if set(result) != needed:
        raise RuntimeError("point audit source misses the exact nine-unit collar")
    return result


def load_h_rows() -> list[dict]:
    needed = {Fraction(value) for value in range(5692, 5711)}
    rows = []
    with h_source.DEFAULT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            left = Fraction(row["target_t_left"])
            if left not in needed:
                continue
            rows.append(
                {
                    "target_t_left": left,
                    "target_t_right": Fraction(row["target_t_right"]),
                    "H": {
                        degree: interval_from_text(row["h_derivatives"][str(degree)])
                        for degree in range(2, 25)
                    },
                }
            )
    rows.sort(key=lambda row: row["target_t_left"])
    if {row["target_t_left"] for row in rows} != needed:
        raise RuntimeError("localized audit source misses the H2-H24 collar")
    return rows


def main() -> int:
    flint.ctx.prec = core.PRECISION_BITS
    issues = []
    if core.REMAINDER_ORDER + 18 != 24:
        issues.append("eighth-stage derivative budget no longer ends at H24")
    if not core.analytic_x_floor_point(Fraction(1252)) > 0:
        issues.append("point X floor is not positive at its first index")
    if not core.analytic_x_floor_interval(Fraction(5701), Fraction(5702)) > 0:
        issues.append("interval X floor is not positive on the audit tile")

    point_source = load_point_source()
    point = core.point_eighth_hierarchy(
        POINT_ANCHOR,
        core.POINT_ORDER,
        support_h_rows=[],
        point_h_source=point_source,
    )
    for name in ("B", "J", "R", "S", "T", "U", "V", "W", "X"):
        if not point[name][0] > 0:
            issues.append(f"point hierarchy lost {name}")
    point_scaled = 2 * point["y"][2] * flint.arb(POINT_ANCHOR.numerator) ** 2
    if not point_scaled < core.CURVATURE_CONSTANT:
        issues.append("point y curvature exceeds the 6000 target")

    common = core.localized_eighth_component_remainder_bounds(
        INTERVAL_ANCHOR,
        INTERVAL_ANCHOR + 1,
        load_h_rows(),
        core.REMAINDER_ORDER,
    )
    for name, value in common["coordinates"].items():
        if not value > 0:
            issues.append(f"localized hierarchy lost {name}")
    stage_counts = common["stage_rows"]
    if "X" not in stage_counts or stage_counts["X"] < 1:
        issues.append("localized hierarchy produced no eighth-stage rows")

    if issues:
        for issue in issues:
            print(issue)
        print(f"order-eleven localized final-gap core: {len(issues)} issues")
        return 1
    print(
        "validated order-eleven localized final-gap core: "
        "H24 derivative budget, 9 positive point coordinates, "
        "9 positive localized coordinates, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
