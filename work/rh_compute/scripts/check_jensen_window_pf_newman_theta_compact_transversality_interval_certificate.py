#!/usr/bin/env python3
"""Validate the theta compact-transversality interval certificate."""

from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_compact_transversality_interval_certificate as target  # noqa: E402


RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "compact_transversality_interval_certificate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_"
    "compact_transversality_interval_certificate.md"
)

EXPECTED_IDS = [
    "ntctic_01_contact_implication",
    "ntctic_02_first_component_transform",
    "ntctic_03_retained_integral_tail",
    "ntctic_04_bivariate_taylor_enclosure",
    "ntctic_05_partition_integrity",
    "ntctic_06_compact_contact_exclusion",
    "ntctic_07_origin_composition",
    "ntctic_08_high_frequency_handoff",
]


def check_partition(records: list[dict], issues: list[str]) -> None:
    parsed: list[
        tuple[Fraction, Fraction, Fraction, Fraction]
    ] = []
    for index, row in enumerate(records):
        try:
            bounds = (
                Fraction(row["t_low"]),
                Fraction(row["t_high"]),
                Fraction(row["x_low"]),
                Fraction(row["x_high"]),
            )
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"record {index} has invalid rational bounds: {exc}")
            continue
        t_low, t_high, x_low, x_high = bounds
        if not (
            target.TIME_LOWER <= t_low < t_high <= target.TIME_UPPER
        ):
            issues.append(f"record {index} leaves the time rectangle")
        if not (
            target.X_LOWER <= x_low < x_high <= target.X_UPPER
        ):
            issues.append(f"record {index} leaves the frequency rectangle")
        parsed.append(bounds)

    time_boundaries = sorted(
        {bound for box in parsed for bound in box[:2]}
    )
    for time_low, time_high in zip(
        time_boundaries, time_boundaries[1:]
    ):
        if time_low == time_high:
            continue
        midpoint = (time_low + time_high) / 2
        intervals = sorted(
            (x_low, x_high)
            for t_low, t_high, x_low, x_high in parsed
            if t_low <= midpoint <= t_high
        )
        if not intervals:
            issues.append(
                f"partition has an empty time strip at {midpoint}"
            )
            continue
        cursor = target.X_LOWER
        for x_low, x_high in intervals:
            if x_low != cursor:
                issues.append(
                    "partition gap/overlap at "
                    f"t={midpoint}, expected x={cursor}, got {x_low}"
                )
                break
            cursor = x_high
        if cursor != target.X_UPPER:
            issues.append(
                f"partition time strip {midpoint} ends at {cursor}"
            )

    total_area = sum(
        (t_high - t_low) * (x_high - x_low)
        for t_low, t_high, x_low, x_high in parsed
    )
    expected_area = (
        (target.TIME_UPPER - target.TIME_LOWER)
        * (target.X_UPPER - target.X_LOWER)
    )
    if total_area != expected_area:
        issues.append(
            f"partition area mismatch: {total_area} != {expected_area}"
        )


def check_sources(issues: list[str]) -> None:
    theta_path = (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_"
        "theta_curvature_probability_operator_gate.json"
    )
    scout_path = (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_"
        "theta_compact_transversality_scout.json"
    )
    for path, phrases in (
        (
            theta_path,
            (
                "delta_3=1/54900000",
                "|J_(1,t)(x)|>B_0(t,x) or",
                "(J_t(x),J_t'(x))!=(0,0)",
            ),
        ),
        (
            scout_path,
            (
                "248/371925",
                "0<=t<=1/5 and |x|<=1/4",
                "61951",
            ),
        ),
    ):
        try:
            text = json.dumps(
                json.loads(path.read_text(encoding="utf-8"))
            )
        except (OSError, json.JSONDecodeError) as exc:
            issues.append(f"invalid source artifact {path}: {exc}")
            continue
        for phrase in phrases:
            if phrase not in text:
                issues.append(
                    f"source artifact {path.name} missing: {phrase}"
                )


def validate(payload: dict, note: str) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != (
        "jensen_window_pf_newman_theta_"
        "compact_transversality_interval_certificate"
    ):
        issues.append("artifact kind changed")
    status = payload.get("status", "")
    for phrase in (
        "rigorous Arb/Taylor contact exclusion for |x|<=38",
        "high-frequency transversality target",
        "Lambda<=0, and RH remain open",
    ):
        if phrase not in status:
            issues.append(f"status boundary missing: {phrase}")

    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order changed")
    for row in rows[:-1]:
        if row.get("readiness") != "ready_to_apply":
            issues.append(f"certified row was demoted: {row.get('id')}")
    if rows and rows[-1].get("readiness") != "open":
        issues.append("high-frequency handoff was promoted")

    exact = payload.get("exact", {})
    if exact.get("proved_rectangle") != (
        "0<=t<=1/5 and 1/4<=x<=38"
    ):
        issues.append("proved rectangle changed")
    tail = exact.get("tail_proof", {})
    if tail.get("radius_used") != "1e-800":
        issues.append("tail radius changed")
    for phrase in ("u>=2", "0<=m<=9"):
        if phrase not in tail.get("range", ""):
            issues.append(f"tail range missing: {phrase}")
    if "<10^-800" not in tail.get("moment_tail", ""):
        issues.append("analytic tail comparison missing")
    exp_eight_lower = sum(
        Fraction(8**k, math.factorial(k)) for k in range(8)
    )
    if not exp_eight_lower > 1000:
        issues.append("independent exp(8) lower audit failed")
    if not Fraction(22, 7) ** 2 < 10:
        issues.append("independent pi-squared upper audit failed")
    if not 20 * math.factorial(9) * 10**800 < 2**2979:
        issues.append("independent tail-radius integer audit failed")
    taylor = exact.get("taylor_enclosure", {})
    if taylor.get("time_order") != 2 or taylor.get("x_order") != 3:
        issues.append("Taylor orders changed")

    certificate = payload.get("certificate", {})
    expected_summary = {
        "precision_bits": 160,
        "integration_cutoff": "2",
        "tail_radius": "1e-800",
        "initial_time_step": "1/50",
        "initial_x_step": "1/5",
        "initial_boxes": 1890,
        "evaluated_boxes": 1910,
        "certified_boxes": 1900,
        "subdivisions": 10,
        "unresolved_boxes": 0,
        "maximum_depth": 1,
    }
    for key, expected in expected_summary.items():
        if certificate.get(key) != expected:
            issues.append(
                f"certificate {key} changed: "
                f"{certificate.get(key)!r} != {expected!r}"
            )
    if certificate.get("branch_counts") != {
        "value": 1115,
        "derivative": 785,
    }:
        issues.append("branch counts changed")
    if certificate.get("below_normal_priority_applied") is not True:
        issues.append("below-normal priority was not applied")

    try:
        minimum = target.arb(
            certificate["minimum_certified_ratio_lower"]
        ).lower()
        if not minimum > target.arb("1.75"):
            issues.append("minimum certified ratio lost its margin")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid minimum ratio: {exc}")

    records = certificate.get("records", [])
    if len(records) != certificate.get("certified_boxes"):
        issues.append("record count differs from certified-box count")
    for index, row in enumerate(records):
        if row.get("certified") is not True:
            issues.append(f"record {index} is not certified")
            continue
        if row.get("branch") not in {"value", "derivative"}:
            issues.append(f"record {index} has invalid branch")
        try:
            ratio = target.arb(row["certified_ratio_lower"]).lower()
            if not ratio > 1:
                issues.append(f"record {index} has non-strict ratio")
        except (KeyError, ValueError) as exc:
            issues.append(f"record {index} has invalid ratio: {exc}")
    check_partition(records, issues)
    check_sources(issues)

    required_note = (
        "# Newman Theta Compact-Transversality Interval Certificate",
        "## Certified Theorem",
        "0<=t<=1/5 and |x|<=38",
        "(H_t(x),H_t'(x))!=(0,0)",
        "248/371925",
        "## Interval Method",
        "1e-800",
        "## Partition",
        "certified leaf boxes=1900",
        "unresolved boxes=0",
        "## Proof Boundary",
        "does not imply `Lambda<=0` or RH",
    )
    lower_note = note.lower()
    for phrase in required_note:
        if phrase.lower() not in lower_note:
            issues.append(f"note missing phrase: {phrase}")
    return issues


def main() -> int:
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note = NOTE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as exc:
        print(f"compact transversality certificate: source failure: {exc}")
        return 1

    issues = validate(payload, note)
    rebuilt = target.build_payload(progress=False)
    if payload != rebuilt:
        issues.append("stored certificate differs from full Arb replay")
    success = target.success_line(payload)
    if success not in note:
        issues.append("success line missing from note")

    if issues:
        print(
            "compact transversality interval certificate: "
            f"{len(issues)} issues"
        )
        for issue in issues:
            print(f"- {issue}")
        return 1

    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
