#!/usr/bin/env python3
"""Validate the theta compact-transversality theorem/scout handoff."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_compact_transversality_scout as target  # noqa: E402


RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_compact_transversality_scout.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_compact_transversality_scout.md"
)

EXPECTED_IDS = [
    "ntcts_01_positive_phi_kernel",
    "ntcts_02_zeroth_moment_floor",
    "ntcts_03_second_moment_ceiling",
    "ntcts_04_near_origin_no_contact",
    "ntcts_05_component_contact_test",
    "ntcts_06_compact_grid_scout",
    "ntcts_07_quadrature_crosscheck",
    "ntcts_08_outer_scope_guard",
    "ntcts_09_dominant_saddle_composition",
    "ntcts_10_oscillatory_zeta_composition",
    "ntcts_11_critical_phase_target",
    "ntcts_12_compact_interval_upgrade",
]

SUCCESS = (
    "validated Newman theta compact-transversality scout: 12 rows, "
    "0 issues, 3 exact moment/kernel inequalities, "
    "1 uniform near-origin no-contact theorem, "
    "61951 compact diagnostic points, 0 compact grid failures, "
    "1 independent quadrature crosscheck, 1 outer nonpromotion guard, "
    "2 high-frequency composition theorems, 2 open certification targets"
)


def check_external_partition(issues: list[str]) -> None:
    sources = {
        "dominant": (
            REPO_ROOT
            / "work/rh_compute/results/"
            "jensen_window_pf_newman_polymath15_"
            "dominant_saddle_global_ray_certificate.json"
        ),
        "oscillatory": (
            REPO_ROOT
            / "work/rh_compute/results/"
            "jensen_window_pf_newman_polymath15_"
            "oscillatory_zeta_handoff_theorem.json"
        ),
        "critical": (
            REPO_ROOT
            / "work/rh_compute/results/"
            "jensen_window_pf_newman_polymath15_"
            "critical_transversality_target.json"
        ),
    }
    required = {
        "dominant": ("t*L>=25", "L_t(x)>0"),
        "oscillatory": (
            "4911678521/1933561194",
            "For every epsilon>0",
        ),
        "critical": (
            "0<tL<=c_*+o(1)",
            "compact no-double-zero certificates",
        ),
    }
    for name, path in sources.items():
        if not path.exists():
            issues.append(f"missing composition source: {path}")
            continue
        try:
            text = json.dumps(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError) as exc:
            issues.append(f"invalid composition source {path}: {exc}")
            continue
        for phrase in required[name]:
            if phrase not in text:
                issues.append(
                    f"{name} composition source missing phrase: {phrase}"
                )


def validate(payload: dict, note: str) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != (
        "jensen_window_pf_newman_theta_compact_transversality_scout"
    ):
        issues.append("artifact kind changed")
    status = payload.get("status", "")
    for phrase in (
        "exact uniform near-origin no-contact theorem",
        "non-rigorous compact first-block scout",
        "not a proof of Lambda<=0 or RH",
    ):
        if phrase not in status:
            issues.append(f"status boundary missing: {phrase}")

    rows = payload.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append("row ids or order changed")
    readiness = {row.get("id"): row.get("readiness") for row in rows}
    for row_id in EXPECTED_IDS[:5]:
        if readiness.get(row_id) != "ready_to_apply":
            issues.append(f"exact row is not ready to apply: {row_id}")
    if readiness.get("ntcts_06_compact_grid_scout") != "diagnostic_only":
        issues.append("compact point scout was promoted")
    if readiness.get("ntcts_07_quadrature_crosscheck") != "diagnostic_only":
        issues.append("quadrature agreement was promoted")
    if readiness.get("ntcts_08_outer_scope_guard") != "guard_validated":
        issues.append("outer scope guard readiness changed")
    for row_id in EXPECTED_IDS[-2:]:
        if readiness.get(row_id) != "open":
            issues.append(f"open target was promoted: {row_id}")

    exact = payload.get("exact", {})
    near = exact.get("near_origin_theorem", {})
    if near.get("margin") != "248/371925":
        issues.append("near-origin rational margin changed")
    if near.get("interval") != "0<=t<=1/5 and |x|<=1/4":
        issues.append("near-origin theorem interval changed")
    if Fraction(9, 2900) - Fraction(5, 2052) != Fraction(
        248, 371925
    ):
        issues.append("independent rational margin audit failed")
    if not Fraction(5, 2052) < Fraction(9, 2900):
        issues.append("near-origin positivity comparison failed")
    audit = exact.get("rational_audit", {})
    if not Fraction(audit.get("exp_one_25_upper", "1")) < Fraction(25, 24):
        issues.append("exp(1/25) upper audit failed")
    if not Fraction(audit.get("endpoint_exponent_upper", "1")) < Fraction(
        10, 3
    ):
        issues.append("endpoint exponent audit failed")
    if not Fraction(audit.get("exp_ten_thirds_upper", "30")) < 29:
        issues.append("exp(10/3) upper audit failed")
    if not Fraction(audit.get("exp_three_lower", "0")) > 20:
        issues.append("exp(3) lower audit failed")

    scout = payload.get("scout", {})
    method = scout.get("method", {})
    if "diagnostic only" not in method.get("status", ""):
        issues.append("quadrature method lost diagnostic boundary")
    compact = scout.get("compact_grid", {})
    expected_grid = {
        "time_min": 0.0,
        "time_max": 0.2,
        "time_step": 0.005,
        "time_rows": 41,
        "x_min": 0.25,
        "x_max": 38.0,
        "x_step": 0.025,
        "x_rows": 1511,
        "total_points": 61951,
    }
    for key, expected in expected_grid.items():
        if compact.get(key) != expected:
            issues.append(
                f"compact grid {key} changed: "
                f"{compact.get(key)!r} != {expected!r}"
            )
    fine_minimum = compact.get("fine_minimum", {})
    if not 1.84 < float(fine_minimum.get("ratio", 0)) < 1.86:
        issues.append("compact minimum ratio left the audited interval")
    if fine_minimum.get("t") != 0.2 or fine_minimum.get("x") != 38.0:
        issues.append("compact minimum location changed")
    if compact.get("fine_failures") != 0:
        issues.append("fine compact grid contains a failed point")
    if compact.get("coarse_failures") != 0:
        issues.append("coarse compact grid contains a failed point")
    difference = float(
        compact.get("maximum_coarse_fine_ratio_difference", 1)
    )
    if not 0 <= difference < 1e-7:
        issues.append("coarse/fine grid agreement failed")

    outer = scout.get("outer_scope_guard", {})
    if int(outer.get("failure_points", 0)) <= 0:
        issues.append("outer nonpromotion witness disappeared")
    first_failures = outer.get("first_failures", [])
    if len(first_failures) != 5:
        issues.append("outer time-row count changed")
    else:
        for row in first_failures:
            first_x = row.get("first_x")
            if first_x is None or not 39.5 <= float(first_x) <= 39.7:
                issues.append(
                    f"outer first-failure location changed: {row}"
                )
    if "cannot be promoted" not in outer.get("interpretation", ""):
        issues.append("outer nonpromotion wording missing")

    required_note = (
        "# Newman Theta Compact-Transversality Scout",
        "## Exact Near-Origin Theorem",
        "248/371925",
        "## Compact First-Block Scout",
        "61951",
        "not promoted to a compact theorem",
        "## Scope Guard",
        "## Global Proof Partition",
        "corrected phase transversality open",
        "## Live Handoff",
        "not a proof",
    )
    lower_note = note.lower()
    for phrase in required_note:
        if phrase.lower() not in lower_note:
            issues.append(f"note missing phrase: {phrase}")

    check_external_partition(issues)
    return issues


def main() -> int:
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note = NOTE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as exc:
        print(f"theta compact-transversality scout: source failure: {exc}")
        return 1

    rebuilt = target.build_payload()
    issues = validate(payload, note)
    if payload != rebuilt:
        issues.append("stored artifact differs from deterministic rebuild")
    if SUCCESS not in note:
        issues.append("success line missing from note")

    if issues:
        print(
            "theta compact-transversality scout: "
            f"{len(issues)} issues"
        )
        for issue in issues:
            print(f"- {issue}")
        return 1

    print(SUCCESS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
