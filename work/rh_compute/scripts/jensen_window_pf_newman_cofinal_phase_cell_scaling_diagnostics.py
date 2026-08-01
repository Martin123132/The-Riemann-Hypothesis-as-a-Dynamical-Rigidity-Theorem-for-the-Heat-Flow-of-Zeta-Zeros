#!/usr/bin/env python3
"""Extract reproducible Q207/Q208 phase-cell scaling diagnostics."""

from __future__ import annotations

from collections import Counter
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from statistics import median


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_cofinal_phase_cell_scaling_diagnostics"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"

Q207_STEM = (
    "jensen_window_pf_newman_theta_forward_six_term_"
    "finite_bridge_interval_certificate"
)
Q208_BOTTOM_STEM = (
    "jensen_window_pf_newman_q208_bottom_phase_cell_certificate"
)
Q208_TOP_STEM = (
    "jensen_window_pf_newman_q208_top_phase_cell_certificate"
)
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
Q207_CACHE = RESULT_DIR / f"{Q207_STEM}.jsonl"
Q207_RESULT = RESULT_DIR / f"{Q207_STEM}.json"
BOTTOM_CACHE = RESULT_DIR / f"{Q208_BOTTOM_STEM}.jsonl"
BOTTOM_RESULT = RESULT_DIR / f"{Q208_BOTTOM_STEM}.json"
TOP_CACHE = RESULT_DIR / f"{Q208_TOP_STEM}.jsonl"
TOP_RESULT = RESULT_DIR / f"{Q208_TOP_STEM}.json"

DATE = "2026-07-25"
SCALES = ("unit", "half_log", "full_log")
ARB_BALL = re.compile(r"^\[(.+?) \+/- (.+?)\]$")
getcontext().prec = 100


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def parse_ball(text: str) -> tuple[Decimal, Decimal]:
    match = ARB_BALL.match(text)
    if match:
        midpoint = Decimal(match.group(1))
        radius = Decimal(match.group(2))
        return midpoint - radius, midpoint + radius
    if text.startswith("[") and text.endswith("]") and "," in text:
        lower, upper = text[1:-1].split(",", maxsplit=1)
        return Decimal(lower.strip()), Decimal(upper.strip())
    point = Decimal(text)
    return point, point


def lower_abs(interval: tuple[Decimal, Decimal]) -> Decimal:
    lower, upper = interval
    if lower <= 0 <= upper:
        return Decimal(0)
    return min(abs(lower), abs(upper))


def upper_abs(interval: tuple[Decimal, Decimal]) -> Decimal:
    lower, upper = interval
    return max(abs(lower), abs(upper))


def midpoint_radius(
    interval: tuple[Decimal, Decimal],
) -> tuple[Decimal, Decimal]:
    lower, upper = interval
    return (lower + upper) / 2, (upper - lower) / 2


def q208_leaves(path: Path) -> list[dict]:
    leaves = [
        leaf
        for record in load_jsonl(path)
        for leaf in record["result"]["certified_leaves"]
    ]
    leaves.sort(key=lambda row: Fraction(row["x_low"]))
    return leaves


def q207_low_time_leaves() -> list[dict]:
    selected: list[dict] = []
    for record in load_jsonl(Q207_CACHE):
        candidates = [
            leaf
            for leaf in record["result"]["certified_leaves"]
            if leaf["t_low"] == "1/1035"
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                "Q207 panel does not have exactly one lowest-time cell"
            )
        selected.append(candidates[0])
    selected.sort(key=lambda row: Fraction(row["x_low"]))
    return selected


def branch_statistics(leaves: list[dict]) -> dict:
    sequence = "".join(
        "V" if leaf["branch"] == "value" else "D"
        for leaf in leaves
    )
    runs = [
        len(match.group(0))
        for match in re.finditer(r"V+|D+", sequence)
    ]
    counts = Counter(leaf["branch"] for leaf in leaves)
    return {
        "cell_count": len(leaves),
        "value_cells": counts["value"],
        "derivative_cells": counts["derivative"],
        "branch_transitions": sum(
            left != right
            for left, right in zip(sequence, sequence[1:])
        ),
        "run_count": len(runs),
        "maximum_run": max(runs),
        "median_run": float(median(runs)),
        "branch_sequence_sha256": sha256(sequence.encode()).hexdigest(),
    }


def jet_scale(name: str, x: float) -> float:
    if name == "unit" or x <= 0:
        return 1.0
    logarithm = math.log(x / (4 * math.pi))
    if name == "half_log":
        return max(1.0, logarithm / 2)
    if name == "full_log":
        return max(1.0, logarithm)
    raise ValueError(f"unknown scale: {name}")


def cell_conditioning(leaves: list[dict], scale_name: str) -> dict:
    rows: list[dict] = []
    for leaf in leaves:
        x_low = float(Fraction(leaf["x_low"]))
        x_high = float(Fraction(leaf["x_high"]))
        x_center = (x_low + x_high) / 2
        scale = jet_scale(scale_name, x_center)
        f_interval = parse_ball(leaf["full_f"])
        fp_interval = parse_ball(leaf["full_f_prime"])
        f_center, f_radius = midpoint_radius(f_interval)
        fp_center, fp_radius = midpoint_radius(fp_interval)
        center_norm = math.hypot(
            float(f_center),
            float(fp_center) / scale,
        )
        clearance = math.hypot(
            float(lower_abs(f_interval)),
            float(lower_abs(fp_interval)) / scale,
        )
        uncertainty = math.hypot(
            float(f_radius),
            float(fp_radius) / scale,
        )
        tail = math.hypot(
            float(upper_abs(parse_ball(leaf["tail_f_upper"]))),
            float(upper_abs(parse_ball(leaf["tail_f_prime_upper"])))
            / scale,
        )
        rows.append(
            {
                "x_center": x_center,
                "scale": scale,
                "clearance": clearance,
                "relative_clearance": (
                    clearance / center_norm if center_norm else 0.0
                ),
                "relative_uncertainty": (
                    uncertainty / center_norm
                    if center_norm
                    else math.inf
                ),
                "tail_domination": (
                    clearance / tail if tail else math.inf
                ),
            }
        )

    min_clearance = min(rows, key=lambda row: row["relative_clearance"])
    max_uncertainty = max(
        rows,
        key=lambda row: row["relative_uncertainty"],
    )
    min_tail = min(rows, key=lambda row: row["tail_domination"])
    return {
        "scale": scale_name,
        "minimum_relative_clearance": min_clearance[
            "relative_clearance"
        ],
        "minimum_relative_clearance_x": min_clearance["x_center"],
        "median_relative_clearance": median(
            row["relative_clearance"] for row in rows
        ),
        "maximum_relative_uncertainty": max_uncertainty[
            "relative_uncertainty"
        ],
        "maximum_relative_uncertainty_x": max_uncertainty["x_center"],
        "minimum_tail_domination_ratio": min_tail["tail_domination"],
        "minimum_tail_domination_x": min_tail["x_center"],
        "origin_free_cells": sum(
            row["clearance"] > 0 for row in rows
        ),
    }


def exact_cross(
    left: tuple[Fraction, Fraction],
    right: tuple[Fraction, Fraction],
) -> Fraction:
    return left[0] * right[1] - left[1] * right[0]


def exact_dot(
    left: tuple[Fraction, Fraction],
    right: tuple[Fraction, Fraction],
) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def exact_open_crossing(
    vertices: list[tuple[Fraction, Fraction]],
) -> int:
    total = 0
    for left, right in zip(vertices, vertices[1:]):
        determinant = exact_cross(left, right)
        if left == (0, 0) or right == (0, 0):
            raise RuntimeError("zero witness")
        if determinant == 0 and exact_dot(left, right) <= 0:
            raise RuntimeError("witness edge meets the origin")
        if left[1] <= 0 < right[1] and determinant > 0:
            total += 1
        elif right[1] <= 0 < left[1] and determinant < 0:
            total -= 1
    return total


def source_witnesses(result: dict) -> list[tuple[Fraction, Fraction]]:
    chain = result["summary"]["open_phase_chain"]
    return [
        (Fraction(point[0]), Fraction(point[1]))
        for point in chain["witnesses"]
    ]


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    index = min(
        len(ordered) - 1,
        max(0, math.ceil(probability * len(ordered)) - 1),
    )
    return ordered[index]


def witness_turning(
    witnesses: list[tuple[Fraction, Fraction]],
    scale_name: str,
) -> dict:
    scaled: list[tuple[float, float]] = []
    for index, point in enumerate(witnesses):
        x = index / 2
        scale = jet_scale(scale_name, x)
        scaled.append((float(point[0]), float(point[1]) / scale))
    increments = [
        math.atan2(
            left[0] * right[1] - left[1] * right[0],
            left[0] * right[0] + left[1] * right[1],
        )
        for left, right in zip(scaled, scaled[1:])
    ]
    absolute = [abs(value) for value in increments]
    return {
        "scale": scale_name,
        "principal_turn_sum_radians": sum(increments),
        "principal_turn_sum_turns": sum(increments) / (2 * math.pi),
        "maximum_absolute_panel_turn": max(absolute),
        "p95_absolute_panel_turn": quantile(absolute, 0.95),
        "panels_over_half_turn": sum(
            value >= math.pi for value in absolute
        ),
        "panels_over_quarter_turn": sum(
            value >= math.pi / 2 for value in absolute
        ),
    }


def circular_difference(left: float, right: float) -> float:
    return math.atan2(
        math.sin(left - right),
        math.cos(left - right),
    )


def witness_comparison(
    bottom: list[tuple[Fraction, Fraction]],
    top: list[tuple[Fraction, Fraction]],
    scale_name: str,
) -> dict:
    if len(bottom) != len(top):
        raise RuntimeError("bottom/top witness counts differ")
    differences: list[tuple[float, float]] = []
    for index, (bottom_point, top_point) in enumerate(zip(bottom, top)):
        x = index / 2
        scale = jet_scale(scale_name, x)
        bottom_angle = math.atan2(
            float(bottom_point[1]) / scale,
            float(bottom_point[0]),
        )
        top_angle = math.atan2(
            float(top_point[1]) / scale,
            float(top_point[0]),
        )
        differences.append(
            (
                x,
                abs(circular_difference(bottom_angle, top_angle)),
            )
        )
    maximum = max(differences, key=lambda row: row[1])
    values = [row[1] for row in differences]
    return {
        "scale": scale_name,
        "maximum_direction_difference": maximum[1],
        "maximum_direction_difference_x": maximum[0],
        "median_direction_difference": median(values),
        "p95_direction_difference": quantile(values, 0.95),
        "right_endpoint_direction_difference": differences[-1][1],
    }


def q207_q208_branch_comparison(
    q207: list[dict],
    bottom: list[dict],
) -> dict:
    q207_by_x = {leaf["x_low"]: leaf for leaf in q207}
    compared = [
        (q207_by_x[leaf["x_low"]], leaf)
        for leaf in bottom
        if leaf["x_low"] in q207_by_x
    ]
    disagreements = [
        right["x_low"]
        for left, right in compared
        if left["branch"] != right["branch"]
    ]
    return {
        "comparable_panels": len(compared),
        "branch_agreements": len(compared) - len(disagreements),
        "branch_disagreements": len(disagreements),
        "disagreement_x_low": disagreements,
        "q207_low_time_value_cells": sum(
            leaf["branch"] == "value" for leaf in q207
        ),
        "q207_low_time_derivative_cells": sum(
            leaf["branch"] == "derivative" for leaf in q207
        ),
        "q207_low_time_cell_t_high_values": sorted(
            {leaf["t_high"] for leaf in q207}
        ),
        "interpretation": (
            "The Q207 rows are two-dimensional cells beginning at "
            "t=1/1035, whereas Q208 rows are fixed-time phase cells. "
            "Branch agreement is a stability diagnostic, not a phase "
            "homotopy certificate."
        ),
    }


def q208_branch_comparison(
    bottom: list[dict],
    top: list[dict],
) -> dict:
    if [row["x_low"] for row in bottom] != [
        row["x_low"] for row in top
    ]:
        raise RuntimeError("bottom/top panel partitions differ")
    disagreements = [
        left["x_low"]
        for left, right in zip(bottom, top)
        if left["branch"] != right["branch"]
    ]
    return {
        "comparable_panels": len(bottom),
        "branch_agreements": len(bottom) - len(disagreements),
        "branch_disagreements": len(disagreements),
        "disagreement_x_low": disagreements,
    }


def shell_geometry() -> dict:
    rows: list[dict] = []
    for stage in (31, 207, 208, 209, 512, 1024, 10000):
        time_value = Fraction(1, 5 * stage)
        right = stage + 38
        product = time_value * right
        delta = time_value - Fraction(1, 5 * (stage + 1))
        logarithm = math.log(right / (4 * math.pi))
        rows.append(
            {
                "stage": stage,
                "bottom_time": str(time_value),
                "right_endpoint": right,
                "time_times_right": str(product),
                "time_step_to_successor": str(delta),
                "L": logarithm,
                "time_times_L": float(time_value) * logarithm,
                "dominant_saddle_tL_at_least_25": (
                    float(time_value) * logarithm >= 25
                ),
            }
        )
    return {
        "definition": "t_j=1/(5*j), R_j=j+38",
        "exact_identities": [
            "t_j*R_j=1/5+38/(5*j)",
            "t_j-t_(j+1)=1/(5*j*(j+1))",
            "R_(j+1)-R_j=1",
        ],
        "limits": [
            "t_j*R_j -> 1/5",
            "t_j*log(R_j/(4*pi)) -> 0",
        ],
        "rows": rows,
        "consequence": (
            "The cofinal bottom edge remains in the nonuniform "
            "t*L->0 layer. The proved dominant-saddle ray t*L>=25 "
            "cannot certify it."
        ),
    }


def render_note(result: dict) -> str:
    branch = result["branch_structure"]
    comparison = result["comparisons"]
    geometry = result["shell_geometry"]
    lines = [
        "# Newman Cofinal Phase-Cell Scaling Diagnostics",
        "",
        f"Date: {DATE}",
        "",
        "Status: reproducible finite diagnostics and exact shell-scaling",
        "identities. This is not a cofinal theorem and not a proof of",
        "`Lambda<=0` or RH.",
        "",
        "## Exact Shell Geometry",
        "",
        "```text",
        "t_j=1/(5*j), R_j=j+38",
        "t_j*R_j=1/5+38/(5*j) -> 1/5",
        "t_j-t_(j+1)=1/(5*j*(j+1))",
        "t_j*log(R_j/(4*pi)) -> 0",
        "```",
        "",
        geometry["consequence"],
        "",
        "## Branch Stability",
        "",
        "| comparison | agreements | panels | disagreements |",
        "|---|---:|---:|---:|",
        (
            "| Q207 lowest-time 2D cells vs Q208 bottom | "
            f"{comparison['q207_q208']['branch_agreements']} | "
            f"{comparison['q207_q208']['comparable_panels']} | "
            f"{comparison['q207_q208']['branch_disagreements']} |"
        ),
        (
            "| Q208 bottom vs Q208 top | "
            f"{comparison['bottom_top']['branch_agreements']} | "
            f"{comparison['bottom_top']['comparable_panels']} | "
            f"{comparison['bottom_top']['branch_disagreements']} |"
        ),
        "",
        "A branch records which coordinate supplied the stronger interval",
        "separation. Agreement is useful evidence about certificate",
        "conditioning, but it is not a homotopy or no-contact theorem.",
        "",
        "## Q208 Cells",
        "",
        "| edge | cells | value | derivative | transitions | max run |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for edge in ("bottom", "top"):
        row = branch[edge]
        lines.append(
            f"| {edge} | {row['cell_count']} | {row['value_cells']} | "
            f"{row['derivative_cells']} | {row['branch_transitions']} | "
            f"{row['maximum_run']} |"
        )
    lines.extend(
        [
            "",
            "Both exact witness chains have positive-ray crossing count",
            f"`{comparison['bottom_crossing']}`. The complete closed",
            "Q208 artifact, not this diagnostic, proves winding zero.",
            "",
            "## Scaling Tests",
            "",
            "The tested positive derivative scales are",
            "`1`, `max(1,L/2)`, and `max(1,L)`, where",
            "`L=log(x/(4*pi))`. Positive scaling preserves the exact",
            "zero set and winding, but the metrics below are floating",
            "conditioning diagnostics only.",
            "",
            "| edge | scale | min relative clearance | max relative uncertainty | min tail domination | max panel turn |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for edge in ("bottom", "top"):
        conditioning = {
            row["scale"]: row
            for row in result["conditioning"][edge]
        }
        turning = {
            row["scale"]: row
            for row in result["witness_turning"][edge]
        }
        for scale_name in SCALES:
            c_row = conditioning[scale_name]
            t_row = turning[scale_name]
            lines.append(
                f"| {edge} | {scale_name} | "
                f"{c_row['minimum_relative_clearance']:.6g} | "
                f"{c_row['maximum_relative_uncertainty']:.6g} | "
                f"{c_row['minimum_tail_domination_ratio']:.6g} | "
                f"{t_row['maximum_absolute_panel_turn']:.6g} |"
            )
    lines.extend(
        [
            "",
            "## Route Decision",
            "",
            "The half-unit phase cells remain well inside the exact",
            "convex-cell interface at Q208, and the six-term arithmetic",
            "tail is far smaller than the certified first-jet clearance.",
            "The hard quantity is retained first-jet separation in the",
            "small-time layer, not omitted-tail control.",
            "",
            "The next theorem should therefore be an adiabatic successor",
            "criterion: transfer the old bottom path across the",
            "`1/(5*j*(j+1))` time collar using heat-equation derivative",
            "bounds, then close the one-unit right strip by a uniform",
            "phase or derivative cone. Computing Q209 alone would not",
            "supply that theorem.",
            "",
            "## Proof Boundary",
            "",
            result["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    q207 = q207_low_time_leaves()
    bottom = q208_leaves(BOTTOM_CACHE)
    top = q208_leaves(TOP_CACHE)
    bottom_result = load_json(BOTTOM_RESULT)
    top_result = load_json(TOP_RESULT)
    bottom_witnesses = source_witnesses(bottom_result)
    top_witnesses = source_witnesses(top_result)
    bottom_crossing = exact_open_crossing(bottom_witnesses)
    top_crossing = exact_open_crossing(top_witnesses)
    if bottom_crossing != top_crossing:
        raise RuntimeError("bottom/top crossing counts differ")

    result = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "reproducible Q207/Q208 finite phase-cell scaling "
            "diagnostics and exact cofinal-shell geometry; no "
            "cofinal theorem"
        ),
        "builder_sha256": file_hash(Path(__file__)),
        "source_sha256": {
            "q207_cache": file_hash(Q207_CACHE),
            "q207_result": file_hash(Q207_RESULT),
            "q208_bottom_cache": file_hash(BOTTOM_CACHE),
            "q208_bottom_result": file_hash(BOTTOM_RESULT),
            "q208_top_cache": file_hash(TOP_CACHE),
            "q208_top_result": file_hash(TOP_RESULT),
        },
        "shell_geometry": shell_geometry(),
        "branch_structure": {
            "bottom": branch_statistics(bottom),
            "top": branch_statistics(top),
        },
        "comparisons": {
            "q207_q208": q207_q208_branch_comparison(q207, bottom),
            "bottom_top": q208_branch_comparison(bottom, top),
            "bottom_crossing": bottom_crossing,
            "top_crossing": top_crossing,
        },
        "conditioning": {
            "bottom": [
                cell_conditioning(bottom, scale_name)
                for scale_name in SCALES
            ],
            "top": [
                cell_conditioning(top, scale_name)
                for scale_name in SCALES
            ],
        },
        "witness_turning": {
            "bottom": [
                witness_turning(bottom_witnesses, scale_name)
                for scale_name in SCALES
            ],
            "top": [
                witness_turning(top_witnesses, scale_name)
                for scale_name in SCALES
            ],
        },
        "bottom_top_witness_direction": [
            witness_comparison(
                bottom_witnesses,
                top_witnesses,
                scale_name,
            )
            for scale_name in SCALES
        ],
        "route_decision": {
            "selected": (
                "adiabatic bottom-path transfer plus a uniform "
                "one-unit right-strip cone"
            ),
            "rejected_as_sufficient": [
                "blind Q209 phase-cell enumeration",
                "dominant-saddle t*L>=25 on the cofinal bottom edge",
                "branch agreement without an origin-free homotopy",
                "floating phase metrics without interval hypotheses",
            ],
            "next_exact_obligation": (
                "Bound the scaled first-jet time displacement across "
                "[t_(j+1),t_j] panel by panel, using "
                "partial_t H=-H_xx and partial_t H_x=-H_xxx, and "
                "prove it is strictly smaller than the old bottom "
                "phase-cell clearance. Separately prove a uniform "
                "origin-free cone on [R_j,R_(j+1)]x[t_(j+1),1/5]."
            ),
        },
        "proof_boundary": (
            "The exact identities describe the chosen exhaustion and "
            "the source-derived counts reproduce finite Q207/Q208 "
            "certificates. The logarithmic scaling, angular, and "
            "conditioning summaries are diagnostics. No uniform "
            "time-collar displacement bound, no all-j right-strip "
            "cone, no Q209 theorem, no cofinal boundary theorem, no "
            "Lambda<=0, and no RH proof is supplied."
        ),
    }
    RESULT_PATH.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE_PATH.write_text(render_note(result), encoding="utf-8")
    print(
        "wrote cofinal phase-cell scaling diagnostics: "
        f"{len(q207)} Q207 panels, {len(bottom)} Q208 bottom cells, "
        f"{len(top)} Q208 top cells, crossing={bottom_crossing}"
    )


if __name__ == "__main__":
    main()
