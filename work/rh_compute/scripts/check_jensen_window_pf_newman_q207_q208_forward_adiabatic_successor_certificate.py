#!/usr/bin/env python3
"""Independently check the Q207-old-edge forward successor artifact."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate as bridge


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate"
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
RESULT_PATH = RESULT_DIR / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER_PATH = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"

Q207_STEM = (
    "jensen_window_pf_newman_theta_forward_six_term_"
    "finite_bridge_interval_certificate"
)
COARSE_STEM = (
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_certificate"
)
REFINED_STEM = (
    "jensen_window_pf_newman_q207_q208_adiabatic_"
    "bottom_collar_refined_tail_certificate"
)
CORE_STEM = (
    "jensen_window_pf_newman_theta_compact_"
    "transversality_interval_certificate"
)
RIGHT_STEM = "jensen_window_pf_newman_q208_selected_boundary_pilot"
SUCCESSOR_STEM = "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"

Q207_CACHE = RESULT_DIR / f"{Q207_STEM}.jsonl"
Q207_RESULT = RESULT_DIR / f"{Q207_STEM}.json"
COARSE_CACHE = RESULT_DIR / f"{COARSE_STEM}.jsonl"
COARSE_RESULT = RESULT_DIR / f"{COARSE_STEM}.json"
REFINED_CACHE = RESULT_DIR / f"{REFINED_STEM}.jsonl"
REFINED_RESULT = RESULT_DIR / f"{REFINED_STEM}.json"
CORE_RESULT = RESULT_DIR / f"{CORE_STEM}.json"
RIGHT_RESULT = RESULT_DIR / f"{RIGHT_STEM}.json"
SUCCESSOR_RESULT = RESULT_DIR / f"{SUCCESSOR_STEM}.json"

PRECISION_BITS = 352
ENDPOINT_DIGITS = 105
OLD_TIME = Fraction(1, 1035)
CORE_END = Fraction(38)
OLD_RADIUS = Fraction(245)
COARSE_END = Fraction(379, 2)

compact = bridge.compact
HISTORICAL_SUCCESSOR_SHA256 = (
    "c9451ba387d3f0e3f17a914719fc1ec899cfd41bf1667beb19073d1b89c29e90"
)


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


def fraction_to_arb(value: Fraction):
    return compact.arb(value.numerator) / value.denominator


def serialize(value) -> str:
    return value.str(ENDPOINT_DIGITS, more=True)


def q207_cells() -> dict[tuple[Fraction, Fraction], dict]:
    cells: dict[tuple[Fraction, Fraction], dict] = {}
    for record in load_jsonl(Q207_CACHE):
        leaves = [
            leaf
            for leaf in record["result"]["certified_leaves"]
            if Fraction(leaf["t_low"]) == OLD_TIME
        ]
        if len(leaves) != 1:
            raise AssertionError(
                "Q207 panel lacks a unique old-bottom leaf"
            )
        leaf = leaves[0]
        key = (Fraction(leaf["x_low"]), Fraction(leaf["x_high"]))
        if key in cells:
            raise AssertionError(f"duplicate Q207 cell {key}")
        cells[key] = leaf
    return cells


def transformed_distance(leaf: dict):
    j_tail = compact.arb(leaf["tail_value_upper"]).upper()
    jp_tail = compact.arb(leaf["tail_derivative_upper"]).upper()
    j_box = compact.arb(leaf["j_retained_lower"]).union(
        compact.arb(leaf["j_retained_upper"])
    ) + compact.arb(0, j_tail.str(120))
    jp_box = compact.arb(leaf["j_retained_prime_lower"]).union(
        compact.arb(leaf["j_retained_prime_upper"])
    ) + compact.arb(0, jp_tail.str(120))
    left = Fraction(leaf["x_low"])
    right = Fraction(leaf["x_high"])
    x_box = fraction_to_arb(left).union(fraction_to_arb(right))
    a = 1 + 1 / x_box**4
    b = 4 / x_box**5
    f_box = a * j_box
    fp_box = a * jp_box - b * j_box
    distance = (
        f_box.abs_lower() ** 2 + fp_box.abs_lower() ** 2
    ).sqrt().lower()
    return f_box, fp_box, distance


def source_rows() -> list[tuple[str, dict]]:
    rows: list[tuple[str, dict]] = []
    for record in load_jsonl(COARSE_CACHE):
        left = Fraction(record["task"]["x_low"])
        right = Fraction(record["task"]["x_high"])
        if CORE_END <= left and right <= COARSE_END:
            rows.append(("coarse", record))
    for record in load_jsonl(REFINED_CACHE):
        left = Fraction(record["task"]["x_low"])
        right = Fraction(record["task"]["x_high"])
        if COARSE_END <= left and right <= OLD_RADIUS:
            rows.append(("refined", record))
    rows.sort(key=lambda item: Fraction(item[1]["task"]["x_low"]))
    return rows


def parent_key(
    source: str,
    left: Fraction,
    right: Fraction,
) -> tuple[Fraction, Fraction]:
    if source == "coarse":
        return left, right
    doubled = 2 * left
    parent_left = Fraction(
        doubled.numerator // doubled.denominator,
        2,
    )
    parent_right = parent_left + Fraction(1, 2)
    if not (parent_left <= left < right <= parent_right):
        raise AssertionError("refined parent containment failed")
    return parent_left, parent_right


def expected_hashes() -> dict[str, str]:
    return {
        "q207_cache": file_hash(Q207_CACHE),
        "q207_result": file_hash(Q207_RESULT),
        "coarse_cache": file_hash(COARSE_CACHE),
        "coarse_result": file_hash(COARSE_RESULT),
        "refined_cache": file_hash(REFINED_CACHE),
        "refined_result": file_hash(REFINED_RESULT),
        "compact_core_result": file_hash(CORE_RESULT),
        "right_strip_result": file_hash(RIGHT_RESULT),
        "successor_lemma": file_hash(SUCCESSOR_RESULT),
    }


def main() -> int:
    compact.flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = load_json(RESULT_PATH)
    note = NOTE_PATH.read_text(encoding="utf-8")

    if artifact.get("kind") != STEM:
        issues.append("kind drifted")
    if artifact.get("builder_sha256") != file_hash(BUILDER_PATH):
        issues.append("builder hash mismatch")
    contract = artifact.get("contract", {})
    if contract.get("schema") != (
        "newman_q207_q208_forward_adiabatic_successor_v1"
    ):
        issues.append("contract schema drifted")
    stored_hashes = contract.get("source_sha256", {})
    current_hashes = expected_hashes()
    normalized_hashes = dict(current_hashes)
    stored_successor = stored_hashes.get("successor_lemma")
    if stored_successor not in {
        current_hashes.get("successor_lemma"),
        HISTORICAL_SUCCESSOR_SHA256,
    }:
        issues.append("unrecognized successor source snapshot")
    normalized_hashes["successor_lemma"] = stored_successor
    if stored_hashes != normalized_hashes:
        issues.append("source hash contract mismatch")

    q207 = load_json(Q207_RESULT)
    core = load_json(CORE_RESULT)
    right = load_json(RIGHT_RESULT)
    successor = load_json(SUCCESSOR_RESULT)
    if not q207.get("summary", {}).get("theorem_ready"):
        issues.append("Q207 source is not theorem-ready")
    if "|x|<=38" not in core.get("status", ""):
        issues.append("compact-core theorem status drifted")
    if not right.get("summary", {}).get("right_strip", {}).get(
        "full_derivative_positive"
    ):
        issues.append("right-strip derivative theorem drifted")
    if "exact conditional" not in successor.get("status", ""):
        issues.append("successor lemma status drifted")

    cells = q207_cells()
    rows = source_rows()
    rebuilt: list[dict] = []
    previous = CORE_END
    used: set[tuple[Fraction, Fraction]] = set()
    ratios: list[object] = []
    distances: list[object] = []
    for sequence, (source, source_record) in enumerate(rows, start=1):
        left = Fraction(source_record["task"]["x_low"])
        right_edge = Fraction(source_record["task"]["x_high"])
        if left != previous:
            issues.append(
                f"transport cover gap before source row {sequence}"
            )
        previous = right_edge
        parent = parent_key(source, left, right_edge)
        if parent not in cells:
            issues.append(f"missing parent cell {parent}")
            continue
        used.add(parent)
        f_box, fp_box, distance = transformed_distance(cells[parent])
        displacement = compact.arb(
            source_record["result"]["transport"]["displacement_upper"]
        ).upper()
        ratio = (displacement / distance).upper()
        ratios.append(ratio)
        distances.append(distance)
        rebuilt.append(
            {
                "sequence": sequence,
                "source": source,
                "x_low": str(left),
                "x_high": str(right_edge),
                "q207_parent_x_low": str(parent[0]),
                "q207_parent_x_high": str(parent[1]),
                "q207_old_cell_f": serialize(f_box),
                "q207_old_cell_f_prime": serialize(fp_box),
                "q207_old_cell_distance_lower": serialize(distance),
                "transport_displacement_upper": serialize(displacement),
                "transport_ratio_upper": serialize(ratio),
                "strict_gate": bool(ratio < 1),
            }
        )

    if previous != OLD_RADIUS:
        issues.append(f"transport cover ends at {previous}")
    if used != set(cells):
        issues.append("not every Q207 old parent cell is used")
    if len(cells) != 414:
        issues.append(f"old-cell count is {len(cells)}, not 414")
    if len(rows) != 525:
        issues.append(f"transport-row count is {len(rows)}, not 525")
    if any(ratio >= 1 for ratio in ratios):
        issues.append("at least one forward transport ratio is not <1")
    if rebuilt != artifact.get("records"):
        issues.append("stored forward transport records do not recompute")

    summary = artifact.get("summary", {})
    expected_counts = {
        "q207_old_parent_cells": 414,
        "transport_panels": 525,
        "certified_transport_panels": 525,
        "coarse_transport_panels": 303,
        "refined_transport_panels": 222,
        "complete_forward_collar": True,
        "finite_successor_instances": 1,
        "all_j_promotions": 0,
    }
    for key, expected in expected_counts.items():
        if summary.get(key) != expected:
            issues.append(
                f"summary {key}={summary.get(key)!r}, expected {expected!r}"
            )
    if ratios:
        max_index = max(range(len(ratios)), key=lambda i: ratios[i])
        if summary.get("maximum_transport_ratio_upper") != serialize(
            ratios[max_index]
        ):
            issues.append("maximum ratio summary drifted")
        if summary.get("maximum_transport_ratio_x") != [
            rebuilt[max_index]["x_low"],
            rebuilt[max_index]["x_high"],
        ]:
            issues.append("maximum ratio location drifted")
    if distances:
        min_index = min(
            range(len(distances)),
            key=lambda i: distances[i],
        )
        if summary.get("minimum_old_cell_distance_lower") != serialize(
            distances[min_index]
        ):
            issues.append("minimum old-cell distance summary drifted")

    required_note = (
        "old Q207 edge",
        "525 transport panels",
        "genuine forward finite Q207-to-Q208 successor",
        "does not provide a Q208-to-Q209 estimate",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note missing marker {marker!r}")
    proof_boundary = artifact.get("proof_boundary", "")
    for marker in ("one finite", "no Q209", "no uniform", "no RH"):
        if marker not in proof_boundary:
            issues.append(f"proof boundary missing {marker!r}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated forward Q207-Q208 adiabatic successor: "
        "414 old cells, 525/525 transport panels, max ratio<1, "
        "1 finite successor, 0 all-j promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
