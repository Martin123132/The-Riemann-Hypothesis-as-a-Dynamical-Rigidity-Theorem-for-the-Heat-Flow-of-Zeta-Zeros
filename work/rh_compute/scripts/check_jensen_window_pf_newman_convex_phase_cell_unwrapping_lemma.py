#!/usr/bin/env python3
"""Validate the convex-cell phase-unwrapping and exact winding lemma."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma as lemma


EXPECTED_IDS = [
    "npcu_01_projection_separation",
    "npcu_02_rectangle_projection",
    "npcu_03_cell_to_chord",
    "npcu_04_witness_polygon",
    "npcu_05_exact_crossing_count",
    "npcu_06_phase_unwrapping",
    "npcu_07_interval_contract",
    "npcu_08_axis_safe_proxy",
    "npcu_09_q208_right_edge",
    "npcu_10_q208_handoff",
]


def independent_winding(vertices: list[lemma.Point]) -> int:
    """Independent exact winding implementation using the negative ray."""
    winding = 0
    for index, left in enumerate(vertices):
        right = vertices[(index + 1) % len(vertices)]
        determinant = lemma.cross(left, right)
        if determinant == 0 and lemma.dot(left, right) <= 0:
            raise ValueError("edge meets origin")
        if left[1] < 0 <= right[1] and determinant > 0:
            winding += 1
        elif right[1] < 0 <= left[1] and determinant < 0:
            winding -= 1
    return winding


def geometry_audit() -> list[str]:
    issues: list[str] = []
    fixtures: list[lemma.Rectangle] = [
        (Fraction(2), Fraction(5), Fraction(-3), Fraction(4)),
        (Fraction(-5), Fraction(-2), Fraction(-4), Fraction(6)),
        (Fraction(-3), Fraction(7), Fraction(2), Fraction(5)),
        (Fraction(-4), Fraction(8), Fraction(-6), Fraction(-1)),
        (Fraction(2), Fraction(5), Fraction(3), Fraction(7)),
    ]
    for index, rectangle in enumerate(fixtures):
        point = lemma.project_zero_to_rectangle(rectangle)
        norm_squared = lemma.dot(point, point)
        if norm_squared <= 0:
            issues.append(f"projection fixture {index} has zero distance")
            continue
        x_low, x_high, y_low, y_high = rectangle
        for vertex in (
            (x_low, y_low),
            (x_low, y_high),
            (x_high, y_low),
            (x_high, y_high),
        ):
            if lemma.dot(point, vertex) < norm_squared:
                issues.append(
                    f"projection separation failed at fixture {index}"
                )
    if not lemma.rectangle_contains_origin(
        (Fraction(-1), Fraction(2), Fraction(-3), Fraction(4))
    ):
        issues.append("origin-containing rectangle was not detected")
    if lemma.rectangle_contains_origin(
        (Fraction(1), Fraction(2), Fraction(-3), Fraction(4))
    ):
        issues.append("origin-free rectangle was rejected")
    return issues


def polygon_audit() -> list[str]:
    issues: list[str] = []
    cases: list[tuple[str, list[lemma.Point], int]] = [
        (
            "ccw",
            [
                (Fraction(1), Fraction(-1)),
                (Fraction(1), Fraction(1)),
                (Fraction(-1), Fraction(1)),
                (Fraction(-1), Fraction(-1)),
            ],
            1,
        ),
        (
            "cw",
            [
                (Fraction(1), Fraction(-1)),
                (Fraction(-1), Fraction(-1)),
                (Fraction(-1), Fraction(1)),
                (Fraction(1), Fraction(1)),
            ],
            -1,
        ),
        (
            "zero",
            [
                (Fraction(2), Fraction(1)),
                (Fraction(-2), Fraction(1)),
                (Fraction(0), Fraction(3)),
            ],
            0,
        ),
        (
            "axis_vertices",
            [
                (Fraction(2), Fraction(0)),
                (Fraction(0), Fraction(2)),
                (Fraction(-2), Fraction(0)),
                (Fraction(0), Fraction(-2)),
            ],
            1,
        ),
    ]
    for name, vertices, expected in cases:
        direct = lemma.polygon_winding(vertices)
        independent = independent_winding(vertices)
        if direct != expected or independent != expected:
            issues.append(
                f"{name} winding mismatch: {direct}, {independent}, "
                f"expected {expected}"
            )
        if name != "axis_vertices":
            cells = lemma.edge_rectangles(vertices)
            certified = lemma.certified_polygon_winding(cells, vertices)
            if certified != expected:
                issues.append(f"{name} certified winding mismatch")
        else:
            try:
                lemma.certified_polygon_winding(
                    lemma.edge_rectangles(vertices),
                    vertices,
                )
                issues.append(
                    "coarse axis-vertex rectangles containing zero were accepted"
                )
            except ValueError:
                pass

    try:
        lemma.polygon_winding(
            [
                (Fraction(1), Fraction(0)),
                (Fraction(-1), Fraction(0)),
                (Fraction(0), Fraction(1)),
            ]
        )
        issues.append("origin-crossing edge was accepted")
    except ValueError:
        pass

    good_vertices = cases[0][1]
    bad_cells = lemma.edge_rectangles(good_vertices)
    bad_cells[0] = (
        Fraction(-1),
        Fraction(1),
        Fraction(-1),
        Fraction(1),
    )
    try:
        lemma.certified_polygon_winding(bad_cells, good_vertices)
        issues.append("origin-containing cell was accepted")
    except ValueError:
        pass

    disconnected = lemma.edge_rectangles(good_vertices)
    disconnected[0] = (
        Fraction(2),
        Fraction(3),
        Fraction(-1),
        Fraction(1),
    )
    try:
        lemma.certified_polygon_winding(disconnected, good_vertices)
        issues.append("missing adjacent-cell witness was accepted")
    except ValueError:
        pass
    return issues


def independent_q208_audit() -> tuple[list[str], dict]:
    issues: list[str] = []
    lemma.q208.bridge.compact.flint.ctx.prec = 352
    source = json.loads(lemma.Q208_RESULT.read_text(encoding="utf-8"))
    right_records = [
        record
        for record in source.get("records", [])
        if record.get("task", {}).get("region")
        == "full_time_right_strip"
    ]
    strip_leaves = [
        leaf
        for record in right_records
        for leaf in record.get("result", {}).get("certified_leaves", [])
    ]
    edge_leaves = [
        leaf
        for record in right_records
        if record.get("task", {}).get("x_high") == "246"
        for leaf in record.get("result", {}).get("certified_leaves", [])
    ]
    arb = lemma.q208.bridge.compact.arb
    lowers = [
        arb(leaf["j_retained_prime_lower"]).lower()
        - arb(leaf["tail_derivative_upper"]).upper()
        for leaf in strip_leaves
    ]
    if len(strip_leaves) != 40:
        issues.append(f"Q208 strip count is {len(strip_leaves)}, not 40")
    if len(edge_leaves) != 20:
        issues.append(f"Q208 edge count is {len(edge_leaves)}, not 20")
    if not lowers or not all(lower > 0 for lower in lowers):
        issues.append("Q208 strip has a nonpositive derivative lower bound")

    x = arb(246)
    scale = 1 + 1 / x**4
    shear = 4 / x**5
    proxy_lowers = []
    for leaf in edge_leaves:
        j_box = arb(leaf["j_retained_lower"]).union(
            arb(leaf["j_retained_upper"])
        ) + arb(
            0,
            arb(leaf["tail_value_upper"]).upper().str(120),
        )
        j_prime_box = arb(leaf["j_retained_prime_lower"]).union(
            arb(leaf["j_retained_prime_upper"])
        ) + arb(
            0,
            arb(leaf["tail_derivative_upper"]).upper().str(120),
        )
        proxy_lowers.append(
            (scale * j_prime_box - shear * j_box).lower()
        )
    if not proxy_lowers or not all(lower > 0 for lower in proxy_lowers):
        issues.append("Q208 transformed proxy derivative is not positive")

    intervals = sorted(
        (
            Fraction(leaf["t_low"]),
            Fraction(leaf["t_high"]),
        )
        for leaf in edge_leaves
    )
    if not intervals:
        issues.append("Q208 right edge has no intervals")
    else:
        if intervals[0][0] != Fraction(1, 1040):
            issues.append("Q208 right-edge lower endpoint drifted")
        if intervals[-1][1] != Fraction(1, 5):
            issues.append("Q208 right-edge upper endpoint drifted")
        if any(
            left[1] != right[0]
            for left, right in zip(intervals, intervals[1:])
        ):
            issues.append("Q208 right-edge cells are not contiguous")
    return issues, {
        "strip_cells": len(strip_leaves),
        "edge_cells": len(edge_leaves),
        "all_derivative_positive": bool(lowers)
        and all(lower > 0 for lower in lowers),
        "all_proxy_derivative_positive": bool(proxy_lowers)
        and all(lower > 0 for lower in proxy_lowers),
    }


def validate(path: Path) -> list[str]:
    issues = geometry_audit() + polygon_audit()
    q208_issues, q208_independent = independent_q208_audit()
    issues.extend(q208_issues)

    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != lemma.STEM:
        issues.append("artifact kind drifted")
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order drifted")
    if artifact.get("exact") != lemma.build_exact():
        issues.append("exact theorem payload drifted")
    if artifact.get("projection_audits") != lemma.rectangle_projection_audits():
        issues.append("projection audit payload drifted")
    if artifact.get("polygon_audits") != lemma.exact_polygon_audits():
        issues.append("polygon audit payload drifted")

    roles = {
        "exact_convex_geometry_lemma": 1,
        "exact_interval_geometry": 1,
        "exact_homotopy_lemma": 1,
        "exact_homotopy_composition": 1,
        "exact_integer_algorithm": 1,
        "exact_phase_lift_equivalence": 1,
        "rigorous_numerical_contract": 1,
        "exact_orientation_preserving_proxy": 1,
        "rigorous_source_composition": 1,
        "open_theorem_target": 1,
    }
    actual_roles = {
        role: sum(row.get("role") == role for row in rows)
        for role in roles
    }
    if actual_roles != roles:
        issues.append(f"row roles drifted: {actual_roles}")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 1 or open_rows[0].get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("Q208 bottom target was promoted")

    right = artifact.get("exact", {}).get("right_edge_calibration", {})
    if right.get("right_strip_cells") != q208_independent["strip_cells"]:
        issues.append("stored Q208 strip count drifted")
    if right.get("right_edge_cells") != q208_independent["edge_cells"]:
        issues.append("stored Q208 edge count drifted")
    if not q208_independent["all_derivative_positive"]:
        issues.append("independent Q208 positivity audit failed")
    if not q208_independent["all_proxy_derivative_positive"]:
        issues.append("independent Q208 proxy positivity audit failed")
    if right.get("phase_half_plane") != (
        "Im(F_t(246)+i*F_t'(246))>0"
    ):
        issues.append("Q208 half-plane statement drifted")

    exact_text = json.dumps(artifact.get("exact", {}))
    for marker in (
        "p dot v >= ||p||^2>0",
        "relative endpoints",
        "rational witness polygon",
        "det(p,q)>0",
        "unique 2pi translate",
        "whole continuous path segment",
        "positive determinant",
        "t=1/1040",
        "closed witness-polygon winding",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    status = artifact.get("status", "")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "exact convex-cell",
        "Q208 right-edge calibration",
        "remain open",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "does not certify",
        "complete bottom edge",
        "closed Q208 boundary winding",
        "Lambda<=0",
        "RH",
        "Clay prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = lemma.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Convex Cell Lemma",
        "Exact Winding",
        "Q208 Right Edge",
        "Live Handoff",
        "No floating-point",
        "globally regular",
        "closed-boundary winding.",
        "Q208 remains open",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    issues = validate(lemma.DEFAULT_OUT)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    artifact = json.loads(lemma.DEFAULT_OUT.read_text(encoding="utf-8"))
    print(
        "validated Newman convex phase-cell unwrapping lemma: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        f"{len(artifact['projection_audits'])} projection audits, "
        f"{len(artifact['polygon_audits']) + 1} exact polygon tests, "
        "3 rejection guards, 20 Q208 right-edge cells, "
        "1 open bottom-edge target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
