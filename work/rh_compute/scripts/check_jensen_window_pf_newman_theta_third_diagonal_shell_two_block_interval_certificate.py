#!/usr/bin/env python3
"""Validate the third diagonal-shell two-block interval certificate."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate as target  # noqa: E402


EXPECTED_IDS = [
    "nttdstbic_01_inherited_two_block_tail",
    "nttdstbic_02_slab_partition",
    "nttdstbic_03_full_kernel_sign",
    "nttdstbic_04_core_composition",
    "nttdstbic_05_third_stage_theorem",
    "nttdstbic_06_winding_composition",
    "nttdstbic_07_method_scope",
    "nttdstbic_08_fourth_stage_handoff",
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
        issues.append("fourth-stage handoff was promoted")
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
        "time_lower": "1/15",
        "time_upper": "1/5",
        "x_lower": "38",
        "x_upper": "41",
        "initial_time_step": "1/60",
        "initial_x_step": "1/2",
        "initial_boxes": 48,
        "evaluated_boxes": 48,
        "certified_boxes": 48,
        "subdivisions": 0,
        "unresolved_boxes": 0,
        "maximum_depth": 0,
        "branch_counts": {"value": 0, "derivative": 48},
        "negative_value_boxes": 48,
    }
    for key, expected in expected_summary.items():
        if certificate.get(key) != expected:
            issues.append(
                f"certificate {key} drifted: "
                f"{certificate.get(key)!r} != {expected!r}"
            )
    if certificate.get("below_normal_priority_applied") is not True:
        issues.append("below-normal priority was not applied")
    try:
        minimum = target.second.compact.arb(
            certificate["minimum_value_ratio_lower"]
        ).lower()
        if not minimum > 1000:
            issues.append("minimum value ratio is not above 1000")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid minimum value ratio: {exc}")

    records = certificate.get("records", [])
    if len(records) != 48:
        issues.append("record count mismatch")
    for index, row in enumerate(records):
        if row.get("certified") is not True:
            issues.append(f"record {index} is not certified")
        if row.get("two_block_value_negative") is not True:
            issues.append(f"record {index} lost full-kernel negativity")
        try:
            t_low = Fraction(row["t_low"])
            t_high = Fraction(row["t_high"])
            x_low = Fraction(row["x_low"])
            x_high = Fraction(row["x_high"])
            j_lower = target.second.compact.arb(
                row["j_two_block_lower"]
            )
            j_upper = target.second.compact.arb(
                row["j_two_block_upper"]
            )
            j_prime_lower = target.second.compact.arb(
                row["j_two_block_prime_lower"]
            )
            j_prime_upper = target.second.compact.arb(
                row["j_two_block_prime_upper"]
            )
            tail = target.second.compact.arb(
                row["tail_value_upper"]
            ).upper()
            derivative_tail = target.second.compact.arb(
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
        if not j_lower.upper() <= j_upper.lower():
            issues.append(f"record {index} has reversed value endpoints")
        if not j_upper.upper() < -tail:
            issues.append(f"record {index} does not prove H_t<0")
        if not j_prime_lower.upper() <= j_prime_upper.lower():
            issues.append(
                f"record {index} has reversed derivative endpoints"
            )
        derivative_separated = (
            j_prime_lower.lower() > derivative_tail
            or j_prime_upper.upper() < -derivative_tail
        )
        if row.get("branch") != "derivative":
            issues.append(f"record {index} changed certificate branch")
        elif not derivative_separated:
            issues.append(
                f"record {index} derivative branch is invalid"
            )

    status = payload.get("status", "")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "closing the third",
        "j>=4",
        "Lambda<=0",
        "RH remain open",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "closes Q_3",
        "does not certify Q_4",
        "Lambda<=0",
        "prove RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    required_note = (
        "Inherited Two-Block Theorem",
        "Certified Slab",
        "Third Stage",
        "Proof Boundary",
        "initial boxes=48",
        "H_t(x)<0",
        "Q_3=[1/15,1/4]x[0,41]",
        "Stages `j>=4`",
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
        print(f"third-shell two-block certificate: source failure: {exc}")
        return 1
    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored certificate differs from full Arb replay")
    success = target.success_line(payload)
    if success not in note:
        issues.append("success line missing from note")
    if issues:
        print(f"third-shell two-block certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
