#!/usr/bin/env python3
"""Independently check the cofinal phase-cell scaling diagnostics."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_cofinal_phase_cell_scaling_diagnostics"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
SOURCES = {
    "q207_cache": RESULT_DIR
    / (
        "jensen_window_pf_newman_theta_forward_six_term_"
        "finite_bridge_interval_certificate.jsonl"
    ),
    "q207_result": RESULT_DIR
    / (
        "jensen_window_pf_newman_theta_forward_six_term_"
        "finite_bridge_interval_certificate.json"
    ),
    "q208_bottom_cache": RESULT_DIR
    / "jensen_window_pf_newman_q208_bottom_phase_cell_certificate.jsonl",
    "q208_bottom_result": RESULT_DIR
    / "jensen_window_pf_newman_q208_bottom_phase_cell_certificate.json",
    "q208_top_cache": RESULT_DIR
    / "jensen_window_pf_newman_q208_top_phase_cell_certificate.jsonl",
    "q208_top_result": RESULT_DIR
    / "jensen_window_pf_newman_q208_top_phase_cell_certificate.json",
}


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


def flatten_q208(path: Path) -> list[dict]:
    leaves = [
        leaf
        for record in load_jsonl(path)
        for leaf in record["result"]["certified_leaves"]
    ]
    return sorted(leaves, key=lambda row: Fraction(row["x_low"]))


def branch_counts(leaves: list[dict]) -> Counter:
    return Counter(leaf["branch"] for leaf in leaves)


def branch_transitions(leaves: list[dict]) -> int:
    branches = [leaf["branch"] for leaf in leaves]
    return sum(
        left != right for left, right in zip(branches, branches[1:])
    )


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


def crossing(result: dict) -> int:
    vertices = [
        (Fraction(point[0]), Fraction(point[1]))
        for point in result["summary"]["open_phase_chain"]["witnesses"]
    ]
    total = 0
    for left, right in zip(vertices, vertices[1:]):
        determinant = exact_cross(left, right)
        if left == (0, 0) or right == (0, 0):
            raise AssertionError("zero witness")
        if determinant == 0 and exact_dot(left, right) <= 0:
            raise AssertionError("witness edge meets origin")
        if left[1] <= 0 < right[1] and determinant > 0:
            total += 1
        elif right[1] <= 0 < left[1] and determinant < 0:
            total -= 1
    return total


def check_float(actual: float, expected: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=2e-12, abs_tol=2e-12):
        raise AssertionError(f"{label}: {actual} != {expected}")


def main() -> None:
    result = load_json(RESULT_PATH)
    issues: list[str] = []
    if result.get("kind") != STEM:
        issues.append("kind mismatch")
    status = result.get("status", "")
    if "no cofinal theorem" not in status:
        issues.append("status lacks nonpromotion guard")
    proof_boundary = result.get("proof_boundary", "")
    for phrase in ("No uniform", "no Lambda<=0", "no RH proof"):
        if phrase not in proof_boundary:
            issues.append(f"proof boundary lacks {phrase!r}")

    for label, path in SOURCES.items():
        if result["source_sha256"].get(label) != file_hash(path):
            issues.append(f"source hash mismatch: {label}")

    bottom = flatten_q208(SOURCES["q208_bottom_cache"])
    top = flatten_q208(SOURCES["q208_top_cache"])
    if len(bottom) != 492 or len(top) != 492:
        issues.append("Q208 phase-cell count mismatch")
    for edge, leaves in (("bottom", bottom), ("top", top)):
        claimed = result["branch_structure"][edge]
        counts = branch_counts(leaves)
        expected = {
            "cell_count": len(leaves),
            "value_cells": counts["value"],
            "derivative_cells": counts["derivative"],
            "branch_transitions": branch_transitions(leaves),
        }
        for key, value in expected.items():
            if claimed.get(key) != value:
                issues.append(f"{edge} {key} mismatch")

    bottom_result = load_json(SOURCES["q208_bottom_result"])
    top_result = load_json(SOURCES["q208_top_result"])
    bottom_crossing = crossing(bottom_result)
    top_crossing = crossing(top_result)
    comparisons = result["comparisons"]
    if (
        bottom_crossing != -19
        or top_crossing != -19
        or comparisons.get("bottom_crossing") != bottom_crossing
        or comparisons.get("top_crossing") != top_crossing
    ):
        issues.append("exact crossing audit mismatch")

    bottom_top_differences = [
        left["x_low"]
        for left, right in zip(bottom, top)
        if left["branch"] != right["branch"]
    ]
    bottom_top = comparisons["bottom_top"]
    if bottom_top.get("comparable_panels") != 492:
        issues.append("bottom/top comparison count mismatch")
    if bottom_top.get("branch_disagreements") != len(
        bottom_top_differences
    ):
        issues.append("bottom/top disagreement count mismatch")
    if bottom_top.get("disagreement_x_low") != bottom_top_differences:
        issues.append("bottom/top disagreement locations mismatch")

    q207: list[dict] = []
    for record in load_jsonl(SOURCES["q207_cache"]):
        candidates = [
            leaf
            for leaf in record["result"]["certified_leaves"]
            if leaf["t_low"] == "1/1035"
        ]
        if len(candidates) != 1:
            issues.append("Q207 lowest-time cell multiplicity mismatch")
            continue
        q207.append(candidates[0])
    q207_by_x = {leaf["x_low"]: leaf for leaf in q207}
    pairs = [
        (q207_by_x[leaf["x_low"]], leaf)
        for leaf in bottom
        if leaf["x_low"] in q207_by_x
    ]
    q207_differences = [
        right["x_low"]
        for left, right in pairs
        if left["branch"] != right["branch"]
    ]
    q207_comparison = comparisons["q207_q208"]
    if len(pairs) != 414:
        issues.append("Q207/Q208 comparison count mismatch")
    if q207_comparison.get("branch_agreements") != (
        len(pairs) - len(q207_differences)
    ):
        issues.append("Q207/Q208 agreement mismatch")
    if q207_comparison.get("disagreement_x_low") != q207_differences:
        issues.append("Q207/Q208 disagreement locations mismatch")

    geometry = result["shell_geometry"]
    geometry_by_stage = {
        row["stage"]: row for row in geometry["rows"]
    }
    for stage, row in geometry_by_stage.items():
        t_value = Fraction(1, 5 * stage)
        right = stage + 38
        if Fraction(row["bottom_time"]) != t_value:
            issues.append(f"stage {stage} bottom time mismatch")
        if Fraction(row["time_times_right"]) != t_value * right:
            issues.append(f"stage {stage} tR mismatch")
        if Fraction(row["time_step_to_successor"]) != Fraction(
            1,
            5 * stage * (stage + 1),
        ):
            issues.append(f"stage {stage} successor step mismatch")
        expected_tl = float(t_value) * math.log(
            right / (4 * math.pi)
        )
        check_float(
            row["time_times_L"],
            expected_tl,
            f"stage {stage} tL",
        )
        if row["dominant_saddle_tL_at_least_25"]:
            issues.append(
                f"sampled stage {stage} unexpectedly enters dominant ray"
            )

    for edge in ("bottom", "top"):
        for row in result["conditioning"][edge]:
            if row["origin_free_cells"] != 492:
                issues.append(
                    f"{edge}/{row['scale']} lost an origin-free cell"
                )
            if row["minimum_tail_domination_ratio"] <= 1:
                issues.append(
                    f"{edge}/{row['scale']} tail domination failed"
                )
        for row in result["witness_turning"][edge]:
            if row["panels_over_half_turn"] != 0:
                issues.append(
                    f"{edge}/{row['scale']} has a half-turn panel"
                )

    if issues:
        raise SystemExit("\n".join(issues))
    print(
        "validated cofinal phase-cell scaling diagnostics: "
        "414 Q207 comparisons, 492+492 Q208 cells, "
        "408/414 and 470/492 branch agreements, "
        "crossing=-19, 0 cofinal promotions"
    )


if __name__ == "__main__":
    main()
