#!/usr/bin/env python3
"""Validate the second diagonal-shell first-block route guard."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard as target  # noqa: E402


EXPECTED_IDS = [
    "ntsdsfbrg_01_inherited_disjunction",
    "ntsdsfbrg_02_transition_point",
    "ntsdsfbrg_03_strict_failed_points",
    "ntsdsfbrg_04_route_guard",
    "ntsdsfbrg_05_non_consequence",
    "ntsdsfbrg_06_second_stage_handoff",
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
        issues.append("closed route-guard row was demoted")
    if not rows or rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append("replacement handoff was promoted")
    if payload.get("exact") != target.build_exact():
        issues.append("exact payload drifted")
    if payload.get("source_audit") != target.source_audit():
        issues.append("source audit drifted")

    certificate = payload.get("certificate", {})
    expected_summary = {
        "precision_bits": 160,
        "integration_cutoff": "2",
        "tail_radius": "1e-800",
        "audit_time": "1/10",
        "audit_points": [
            "198/5",
            "397/10",
            "199/5",
            "399/10",
            "40",
        ],
        "evaluated_points": 5,
        "raw_disjunction_certified_points": 1,
        "raw_disjunction_rigorously_false_points": 4,
        "interval_undetermined_points": 0,
        "maximum_failed_ratio_point": "397/10",
        "maximum_failed_ratio_branch": "derivative",
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
        maximum = target.compact.arb(
            certificate["maximum_failed_ratio_upper"]
        ).upper()
        if not maximum < target.compact.arb(24) / 25:
            issues.append("maximum failed ratio is not below 24/25")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid maximum failed ratio: {exc}")

    records = certificate.get("records", [])
    if len(records) != len(target.AUDIT_POINTS):
        issues.append("point record count mismatch")
    for index, (row, expected_x) in enumerate(
        zip(records, target.AUDIT_POINTS)
    ):
        try:
            if Fraction(row["t"]) != target.AUDIT_TIME:
                issues.append(f"record {index} changed time")
            if Fraction(row["x"]) != expected_x:
                issues.append(f"record {index} changed point")
            value_lower = target.compact.arb(
                row["value_ratio_lower"]
            ).lower()
            value_upper = target.compact.arb(
                row["value_ratio_upper"]
            ).upper()
            derivative_lower = target.compact.arb(
                row["derivative_ratio_lower"]
            ).lower()
            derivative_upper = target.compact.arb(
                row["derivative_ratio_upper"]
            ).upper()
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"record {index} is invalid: {exc}")
            continue
        if not value_lower <= value_upper:
            issues.append(f"record {index} value interval inverted")
        if not derivative_lower <= derivative_upper:
            issues.append(f"record {index} derivative interval inverted")
        if index == 0:
            if row.get("classification") != "raw_disjunction_certified":
                issues.append("transition point lost certification")
            if not derivative_lower > 1:
                issues.append("transition derivative ratio is not above one")
        else:
            if (
                row.get("classification")
                != "raw_disjunction_rigorously_false"
            ):
                issues.append(f"record {index} lost strict route failure")
            if not (value_upper < 1 and derivative_upper < 1):
                issues.append(
                    f"record {index} does not reject both branches"
                )

    status = payload.get("status", "")
    boundary = payload.get("proof_boundary", "")
    for marker in ("route guard", "Q_2", "Lambda<=0", "RH remain open"):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "four exact points",
        "does not prove a common zero",
        "close Q_2",
        "Lambda<=0",
        "prove RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    required_note = (
        "Exact Point Audit",
        "Proof Boundary",
        "upper bounds",
        "397/10",
        "raw_disjunction_rigorously_false",
        "Failure of a sufficient bound does not imply",
        "maximum failed upper ratio",
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
        print(f"second diagonal-shell route guard: source failure: {exc}")
        return 1
    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored route guard differs from full Arb replay")
    success = target.success_line(payload)
    if success not in note:
        issues.append("success line missing from note")
    if issues:
        print(f"second diagonal-shell route guard: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
