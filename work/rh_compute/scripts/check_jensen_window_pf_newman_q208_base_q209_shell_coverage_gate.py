#!/usr/bin/env python3
"""Independently validate the Q208-base/Q209 shell coverage gate."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import jensen_window_pf_newman_q208_base_q209_shell_coverage_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q208_base_q209_shell_coverage_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

EXPECTED_COVERAGE_IDS = [
    "q209_cov_01_inherited_p208",
    "q209_cov_02_compact_collar",
    "q209_cov_03_outer_collar",
    "q209_cov_04_new_right_strip",
    "q209_cov_05_high_time",
]
EXPECTED_BOUNDARY_IDS = [
    "q209_edge_01_bottom_compact",
    "q209_edge_02_bottom_outer",
    "q209_edge_03_bottom_new",
    "q209_edge_04_right_low",
    "q209_edge_05_right_high",
    "q209_edge_06_top",
    "q209_edge_07_left",
]
EXPECTED_GATE_IDS = [
    f"q209_gate_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        [
            "handoff_correction",
            "q208_base",
            "time_step",
            "shell_decomposition",
            "area_audit",
            "orientation",
            "inherited_coverage",
            "outer_gap",
            "strip_gap",
            "high_time",
            "time_noninheritance",
            "space_noninheritance",
            "conditional_q209",
            "finite_nonpromotion",
            "ray_route",
        ],
        start=1,
    )
]


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def parse_fraction(value: str) -> Fraction:
    return Fraction(value)


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "stored Q209 shell result", issues)
    if not NOTE.is_file():
        issues.append("missing rendered note")
    if not stored:
        return issues

    rebuilt = gate.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from deterministic reconstruction")
    if stored.get("kind") != STEM:
        issues.append("kind drifted")
    if stored.get("status") != (
        "exact Q208-base/Q209 successor-shell coverage gate "
        "with two explicit open regions"
    ):
        issues.append("status drifted")

    coverage = stored.get("coverage", [])
    boundary = stored.get("boundary", [])
    rows = stored.get("rows", [])
    if [row.get("id") for row in coverage] != EXPECTED_COVERAGE_IDS:
        issues.append("coverage ids or ordering drifted")
    if [row.get("id") for row in boundary] != EXPECTED_BOUNDARY_IDS:
        issues.append("boundary ids or ordering drifted")
    if [row.get("id") for row in rows] != EXPECTED_GATE_IDS:
        issues.append("gate ids or ordering drifted")

    # Reconstruct the linear successor arithmetic independently.
    t_208 = Fraction(1, 5 * 208)
    t_209 = Fraction(1, 5 * 209)
    delta = t_208 - t_209
    if t_208 != Fraction(1, 1040):
        issues.append("independent t_208 arithmetic failed")
    if t_209 != Fraction(1, 1045):
        issues.append("independent t_209 arithmetic failed")
    if delta != Fraction(1, 217360):
        issues.append("independent delta_208 arithmetic failed")
    if 208 + 38 != 246 or 209 + 38 != 247:
        issues.append("independent radii arithmetic failed")

    exact = stored.get("exact", {})
    schedule = exact.get("linear_schedule", {})
    if parse_fraction(schedule.get("delta_208", "0")) != delta:
        issues.append("stored delta_208 is wrong")
    if exact.get("rectangles", {}).get("Q_209") != (
        "[1/1045,1/4]x[0,247]"
    ):
        issues.append("stored Q209 rectangle is wrong")

    # Independently verify the refined low-shell area decomposition.
    low_top = Fraction(1, 5)
    inherited = (low_top - t_208) * 246
    compact_collar = delta * 38
    outer_collar = delta * (246 - 38)
    right_strip = (low_top - t_209) * (247 - 246)
    total = (low_top - t_209) * 247
    if inherited + compact_collar + outer_collar + right_strip != total:
        issues.append("independent low-shell area decomposition failed")
    shell = exact.get("low_shell_decomposition", {})
    expected_areas = {
        "area_P_208": inherited,
        "area_compact_collar": compact_collar,
        "area_outer_collar": outer_collar,
        "area_right_strip": right_strip,
        "area_P_209": total,
    }
    for key, expected in expected_areas.items():
        if parse_fraction(shell.get(key, "0")) != expected:
            issues.append(f"stored shell area failed: {key}")
    if not shell.get("area_sum_verified"):
        issues.append("stored area verification is false")

    # Reconstruct the oriented rectangle polygon in the standard (x,t) order.
    vertices = [
        (Fraction(0), t_209),
        (Fraction(247), t_209),
        (Fraction(247), Fraction(1, 4)),
        (Fraction(0), Fraction(1, 4)),
    ]
    twice_area = sum(
        x_0 * t_1 - t_0 * x_1
        for (x_0, t_0), (x_1, t_1) in zip(
            vertices, vertices[1:] + vertices[:1], strict=True
        )
    )
    if twice_area != 2 * 247 * (Fraction(1, 4) - t_209):
        issues.append("independent oriented-area magnitude failed")
    if twice_area <= 0:
        issues.append("Q209 orientation is not counterclockwise")
    orientation = exact.get("boundary_orientation", {})
    if parse_fraction(orientation.get("twice_signed_area", "0")) != twice_area:
        issues.append("stored oriented area failed")
    if orientation.get("coordinate_order") != "(x,t)":
        issues.append("coordinate orientation drifted")

    # Cell counts follow exactly from the stored half-unit Q208 phase mesh.
    if Fraction(246 - 38, 1) / Fraction(1, 2) != 416:
        issues.append("independent outer half-cell count failed")
    if Fraction(247 - 246, 1) / Fraction(1, 2) != 2:
        issues.append("independent new-strip half-column count failed")
    cell_counts = exact.get("phase_cell_counts", {})
    if cell_counts.get("outer_old_half_unit_cells") != 416:
        issues.append("stored outer half-cell count failed")
    if cell_counts.get("new_right_strip_half_unit_columns") != 2:
        issues.append("stored right-strip column count failed")

    # The coverage map must expose exactly two analytic shell obligations.
    covered = [
        row for row in coverage if row.get("status") == "covered_exact"
    ]
    open_regions = [
        row
        for row in coverage
        if row.get("status") == "open_xi_antecedent"
    ]
    if len(covered) != 3:
        issues.append("expected three covered coverage rows")
    if [row.get("id") for row in open_regions] != [
        "q209_cov_03_outer_collar",
        "q209_cov_04_new_right_strip",
    ]:
        issues.append("the two open shell regions drifted")
    if any(
        row.get("status") != "open_xi_antecedent"
        for row in boundary[1:4]
    ):
        issues.append("expected the three low open boundary traces")
    if any(
        row.get("status") != "covered_exact"
        for row in [boundary[0], *boundary[4:]]
    ):
        issues.append("covered boundary arcs drifted")

    # Independently audit the source domains that forbid automatic reuse.
    q208 = load_json(gate.SOURCES["q208_closed"], "Q208 source", issues)
    bottom = load_json(gate.SOURCES["q208_bottom"], "bottom source", issues)
    right = load_json(gate.SOURCES["q208_right"], "right source", issues)
    forward = load_json(
        gate.SOURCES["q207_q208_forward"], "forward source", issues
    )
    compact = load_json(
        gate.SOURCES["compact_core"], "compact source", issues
    )
    successor = load_json(
        gate.SOURCES["successor_lemma"], "successor source", issues
    )
    scaled = load_json(
        gate.SOURCES["scaled_successor"], "scaled source", issues
    )
    ray = load_json(gate.SOURCES["ray_aligned"], "ray source", issues)
    cofinal = load_json(
        gate.SOURCES["cofinal_target"], "cofinal source", issues
    )
    if q208.get("exact_winding") != 0:
        issues.append("Q208 source winding is not zero")
    if bottom.get("summary", {}).get("tasks_total") != 492:
        issues.append("Q208 bottom source cell count failed")
    if right.get("summary", {}).get("right_strip", {}).get("domain") != (
        "[1/1040,1/5]x[245,246]"
    ):
        issues.append("old right-strip source domain failed")
    if forward.get("contract", {}).get("time_step") != "1/215280":
        issues.append("old forward time-step source failed")
    if forward.get("contract", {}).get("transported_collar") != ["38", "245"]:
        issues.append("old transported-collar source failed")
    compact_exact = compact.get("exact", {})
    if "1/4<=x<=38" not in compact_exact.get("proved_rectangle", ""):
        issues.append("compact outer source failed")
    if "|x|<=1/4" not in compact_exact.get("origin_overlap", ""):
        issues.append("compact origin source failed")
    successor_exact = successor.get("exact", {})
    if "delta_j*M_(j,k)<d_(j,k)" not in successor_exact.get(
        "transport", {}
    ).get("strict_gate", ""):
        issues.append("successor transport source failed")
    if scaled.get("exact_summary", {}).get("all_j_xi_theorem") is not False:
        issues.append("scaled source all-j guard failed")
    if ray.get("summary", {}).get("new_strip_open_antecedents") != 0:
        issues.append("ray new-strip source failed")
    if ray.get("summary", {}).get("old_collar_open_antecedents") != 1:
        issues.append("ray old-collar source failed")
    cofinal_contract = cofinal.get("exact", {}).get("cofinal_contract", {})
    if "Q_207=" not in cofinal_contract.get("current_base", ""):
        issues.append("historical cofinal base source failed")

    guards = exact.get("noninheritance_guards", {})
    if "[1/1040,1/1035]" not in guards.get("time_enclosures", ""):
        issues.append("time noninheritance guard missing old interval")
    if "[1/1045,1/1040]" not in guards.get("time_enclosures", ""):
        issues.append("time noninheritance guard missing new interval")
    if "[245,246]" not in guards.get("space_translation", ""):
        issues.append("space noninheritance guard missing old strip")
    if "[246,247]" not in guards.get("space_translation", ""):
        issues.append("space noninheritance guard missing new strip")

    ray_route = exact.get("ray_aligned_route", {})
    for phrase in [
        "L(x)=log(x/(4*pi))",
        "No polygon",
        "Defining kappa as the quotient",
    ]:
        combined = json.dumps(ray_route, sort_keys=True)
        if phrase not in combined:
            issues.append(f"ray route missing phrase: {phrase}")

    summary = stored.get("summary", {})
    expected_summary = {
        "proved_base_index": 208,
        "first_unresolved_linear_stage": 209,
        "coverage_regions": 5,
        "covered_regions": 3,
        "open_shell_regions": 2,
        "boundary_arcs": 7,
        "open_boundary_arcs": 3,
        "outer_old_half_unit_cells": 416,
        "new_right_strip_half_unit_columns": 2,
        "q209_certified": False,
        "ray_aligned_new_strips_open": 0,
        "ray_aligned_old_collar_antecedents": 1,
        "gate_rows": 15,
    }
    if summary != expected_summary:
        issues.append("summary drifted")

    proof_boundary = stored.get("proof_boundary", "")
    for phrase in [
        "does not certify either open Q209 region",
        "prove a uniform successor theorem",
        "prove Lambda<=0",
        "prove RH",
        "Clay-prize conclusion",
    ]:
        if phrase not in proof_boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    if NOTE.is_file():
        note = NOTE.read_text(encoding="utf-8")
        for phrase in [
            "Q209 not certified",
            "Exactly two low-time shell regions remain open",
            "Neither antecedent is currently proved.",
            "Pi provenance:",
            "No polygon, curvature image",
            "A finite Q209 computation may calibrate or falsify",
        ]:
            if phrase not in note:
                issues.append(f"note missing phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Q208-base/Q209 shell coverage gate: 15 rows, "
        "5 coverage regions, 2 open shell regions, 7 oriented boundary "
        "arcs, 416 outer old half-cells, 2 new strip half-columns, "
        "Q209 not certified"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
