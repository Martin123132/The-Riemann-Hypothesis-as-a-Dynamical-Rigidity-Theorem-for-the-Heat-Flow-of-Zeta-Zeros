#!/usr/bin/env python3
"""Validate the fourth diagonal-shell two-block interval certificate."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate as target  # noqa: E402


EXPECTED_IDS = [
    "ntfdstbic_01_inherited_two_block_tail",
    "ntfdstbic_02_slab_partition",
    "ntfdstbic_03_mixed_branch_transition",
    "ntfdstbic_04_full_kernel_separation",
    "ntfdstbic_05_core_composition",
    "ntfdstbic_06_fourth_stage_theorem",
    "ntfdstbic_07_winding_composition",
    "ntfdstbic_08_fifth_stage_handoff",
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
        issues.append("fifth-stage handoff was promoted")
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
        "time_lower": "1/20",
        "time_upper": "1/5",
        "x_lower": "38",
        "x_upper": "42",
        "initial_time_step": "1/100",
        "initial_x_step": "1/2",
        "initial_boxes": 120,
        "evaluated_boxes": 120,
        "certified_boxes": 120,
        "subdivisions": 0,
        "unresolved_boxes": 0,
        "maximum_depth": 0,
        "branch_counts": {"value": 0, "derivative": 120},
        "negative_value_boxes": 105,
        "value_separated_boxes": 105,
        "derivative_separated_boxes": 120,
        "derivative_only_boxes": 15,
    }
    for key, expected in expected_summary.items():
        if certificate.get(key) != expected:
            issues.append(
                f"certificate {key} drifted: "
                f"{certificate.get(key)!r} != {expected!r}"
            )
    if certificate.get("below_normal_priority_applied") is not True:
        issues.append("below-normal priority was not applied")
    arb = target.third.second.compact.arb
    try:
        minimum = arb(
            certificate["minimum_certified_ratio_lower"]
        ).lower()
        if not minimum > 88000:
            issues.append("minimum disjunction ratio is not above 88000")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid minimum ratio: {exc}")

    records = certificate.get("records", [])
    if len(records) != 120:
        issues.append("record count mismatch")
    transition: list[tuple[Fraction, Fraction]] = []
    negative_count = 0
    value_separated_count = 0
    for index, row in enumerate(records):
        if row.get("certified") is not True:
            issues.append(f"record {index} is not certified")
        try:
            t_low = Fraction(row["t_low"])
            t_high = Fraction(row["t_high"])
            x_low = Fraction(row["x_low"])
            x_high = Fraction(row["x_high"])
            ratio = arb(row["certified_ratio_lower"]).lower()
            value_ratio = arb(row["value_ratio_lower"]).lower()
            j_lower = arb(row["j_two_block_lower"])
            j_upper = arb(row["j_two_block_upper"])
            j_prime_lower = arb(row["j_two_block_prime_lower"])
            j_prime_upper = arb(row["j_two_block_prime_upper"])
            value_tail = arb(row["tail_value_upper"]).upper()
            derivative_tail = arb(
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
        if not ratio > 88000:
            issues.append(f"record {index} lost strict separation")
        if not j_lower.upper() <= j_upper.lower():
            issues.append(f"record {index} has reversed value endpoints")
        if not j_prime_lower.upper() <= j_prime_upper.lower():
            issues.append(
                f"record {index} has reversed derivative endpoints"
            )
        value_separated = (
            j_lower.lower() > value_tail
            or j_upper.upper() < -value_tail
        )
        derivative_separated = (
            j_prime_lower.lower() > derivative_tail
            or j_prime_upper.upper() < -derivative_tail
        )
        negative = j_upper.upper() < -value_tail
        if negative != (row.get("two_block_value_negative") is True):
            issues.append(f"record {index} value-sign flag is invalid")
        if value_separated != (value_ratio > 1):
            issues.append(f"record {index} value ratio is inconsistent")
        negative_count += int(negative)
        value_separated_count += int(value_separated)
        if row.get("branch") != "derivative":
            issues.append(f"record {index} changed certificate branch")
        elif not derivative_separated:
            issues.append(
                f"record {index} derivative branch is invalid"
            )
        if not value_separated:
            transition.append((t_low, t_high))
            if (x_low, x_high) != (Fraction(83, 2), Fraction(42)):
                issues.append(
                    f"record {index} moved the sign-transition cell"
                )
    if negative_count != 105:
        issues.append("independent negative-value count mismatch")
    if value_separated_count != 105:
        issues.append("independent value-separation count mismatch")
    expected_transition = [
        (
            target.TIME_LOWER + i * target.INITIAL_TIME_STEP,
            target.TIME_LOWER + (i + 1) * target.INITIAL_TIME_STEP,
        )
        for i in range(15)
    ]
    if sorted(transition) != expected_transition:
        issues.append("sign-transition time cover drifted")

    status = payload.get("status", "")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "closing the fourth",
        "j>=5",
        "Lambda<=0",
        "RH remain open",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "closes Q_4",
        "does not certify Q_5",
        "Lambda<=0",
        "prove RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    required_note = (
        "Inherited Two-Block Theorem",
        "Certified Slab",
        "Sign Transition",
        "Fourth Stage",
        "Proof Boundary",
        "initial boxes=120",
        "(H_t(x),H_t'(x))!=(0,0)",
        "Q_4=[1/20,1/4]x[0,42]",
        "Stages `j>=5`",
        "not a proof of RH",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    try:
        payload = json.loads(
            target.DEFAULT_OUT.read_text(encoding="utf-8")
        )
        note = target.DEFAULT_NOTE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as exc:
        print(f"fourth-shell two-block certificate: source failure: {exc}")
        return 1
    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored certificate differs from full Arb replay")
    success = target.success_line(payload)
    if success not in note:
        issues.append("success line missing from note")
    if issues:
        print(f"fourth-shell two-block certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
