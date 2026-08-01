#!/usr/bin/env python3
"""Validate the second diagonal-shell two-block interval certificate."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate as target  # noqa: E402


EXPECTED_IDS = [
    "ntsdstbic_01_two_block_identity",
    "ntsdstbic_02_tail_moments",
    "ntsdstbic_03_direct_tail_bars",
    "ntsdstbic_04_cutoff_tail",
    "ntsdstbic_05_slab_partition",
    "ntsdstbic_06_full_kernel_separation",
    "ntsdstbic_07_core_composition",
    "ntsdstbic_08_second_stage_theorem",
    "ntsdstbic_09_winding_composition",
    "ntsdstbic_10_third_stage_handoff",
]


def validate(payload: dict, note: str) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != target.STEM:
        issues.append("artifact kind mismatch")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if any(
        row.get("readiness") != "ready_to_apply" for row in rows[:-1]
    ):
        issues.append("closed theorem row was demoted")
    if not rows or rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append("third-stage handoff was promoted")
    if payload.get("exact") != target.build_exact():
        issues.append("exact payload drifted")
    if payload.get("source_audit") != target.source_audit():
        issues.append("source audit drifted")

    certificate = payload.get("certificate", {})
    expected_summary = {
        "precision_bits": 160,
        "endpoint_serialization_digits": 50,
        "retained_theta_blocks": [1, 2],
        "integration_cutoff": "2",
        "cutoff_tail_radius": "1e-800",
        "tail_moment_zero_upper": "1/10000000000",
        "tail_moment_one_upper": "1/1000000000000",
        "time_lower": "1/10",
        "time_upper": "1/5",
        "x_lower": "38",
        "x_upper": "40",
        "initial_time_step": "1/50",
        "initial_x_step": "2/5",
        "initial_boxes": 25,
        "evaluated_boxes": 25,
        "certified_boxes": 25,
        "subdivisions": 0,
        "unresolved_boxes": 0,
        "maximum_depth": 0,
        "branch_counts": {"value": 0, "derivative": 25},
        "negative_value_boxes": 25,
    }
    for key, expected in expected_summary.items():
        if certificate.get(key) != expected:
            issues.append(
                f"certificate {key} drifted: "
                f"{certificate.get(key)!r} != {expected!r}"
            )
    if certificate.get("below_normal_priority_applied") is not True:
        issues.append("below-normal priority was not applied")
    if certificate.get("evaluated_boxes", 0) < 25:
        issues.append("too few boxes were evaluated")
    if certificate.get("certified_boxes", 0) < 25:
        issues.append("too few slab boxes were certified")
    if certificate.get("maximum_depth", -1) > target.MAX_DEPTH:
        issues.append("maximum depth exceeds the contract")
    try:
        minimum = target.compact.arb(
            certificate["minimum_certified_ratio_lower"]
        ).lower()
        if not minimum > 1000:
            issues.append("minimum normalized ratio is not above 1000")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid minimum ratio: {exc}")
    try:
        minimum_value = target.compact.arb(
            certificate["minimum_value_ratio_lower"]
        ).lower()
        if not minimum_value > 39000:
            issues.append(
                "minimum normalized value ratio is not above 39000"
            )
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid minimum value ratio: {exc}")

    records = certificate.get("records", [])
    for index, row in enumerate(records):
        if row.get("certified") is not True:
            issues.append(f"record {index} is not certified")
        try:
            t_low = Fraction(row["t_low"])
            t_high = Fraction(row["t_high"])
            x_low = Fraction(row["x_low"])
            x_high = Fraction(row["x_high"])
            ratio = target.compact.arb(
                row["certified_ratio_lower"]
            ).lower()
            j_lower = target.compact.arb(
                row["j_two_block_lower"]
            )
            j_upper = target.compact.arb(
                row["j_two_block_upper"]
            )
            j_prime_lower = target.compact.arb(
                row["j_two_block_prime_lower"]
            )
            j_prime_upper = target.compact.arb(
                row["j_two_block_prime_upper"]
            )
            value_tail = target.compact.arb(
                row["tail_value_upper"]
            ).upper()
            derivative_tail = target.compact.arb(
                row["tail_derivative_upper"]
            ).upper()
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"record {index} is invalid: {exc}")
            continue
        if not (
            target.TIME_LOWER <= t_low < t_high <= target.TIME_UPPER
        ):
            issues.append(f"record {index} leaves the time slab")
        if not target.X_LOWER <= x_low < x_high <= target.X_UPPER:
            issues.append(f"record {index} leaves the x slab")
        if not ratio > 1:
            issues.append(f"record {index} lost strict separation")
        if row.get("two_block_value_negative") is not True:
            issues.append(f"record {index} lost the negative value sign")
        if not j_lower.upper() <= j_upper.lower():
            issues.append(f"record {index} has reversed value endpoints")
        if not j_upper.upper() < -value_tail:
            issues.append(
                f"record {index} does not prove full-kernel negativity"
            )
        branch = row.get("branch")
        if branch == "value":
            value_separated = (
                j_lower.lower() > value_tail
                or j_upper.upper() < -value_tail
            )
            if not value_separated:
                issues.append(f"record {index} value branch is invalid")
        elif branch == "derivative":
            if not j_prime_lower.upper() <= j_prime_upper.lower():
                issues.append(
                    f"record {index} has reversed derivative endpoints"
                )
            derivative_separated = (
                j_prime_lower.lower() > derivative_tail
                or j_prime_upper.upper() < -derivative_tail
            )
            if not derivative_separated:
                issues.append(
                    f"record {index} derivative branch is invalid"
                )
        else:
            issues.append(f"record {index} has invalid branch")

    status = payload.get("status", "")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "closing the second",
        "j>=3",
        "Lambda<=0",
        "RH remain open",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "retains theta blocks n=1,2",
        "closes Q_2",
        "does not certify Q_3",
        "Lambda<=0",
        "prove RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    required_note = (
        "Two-Block Replacement",
        "Certified Slab",
        "Second Stage",
        "Proof Boundary",
        "n=2",
        "n>=3",
        "initial boxes=25",
        "minimum normalized value-sign ratio",
        "J_t(x)<0",
        "Q_2=[1/10,1/4]x[0,40]",
        "Stages `j>=3`",
        "not a proof of RH",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    try:
        payload = json.loads(target.DEFAULT_OUT.read_text(encoding="utf-8"))
        note = target.DEFAULT_NOTE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as exc:
        print(f"second-shell two-block certificate: source failure: {exc}")
        return 1
    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored certificate differs from full Arb replay")
    success = target.success_line(payload)
    if success not in note:
        issues.append("success line missing from note")
    if issues:
        print(f"second-shell two-block certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
