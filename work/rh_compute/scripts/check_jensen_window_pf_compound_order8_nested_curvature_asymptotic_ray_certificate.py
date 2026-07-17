#!/usr/bin/env python3
"""Validate the order-eight nested-curvature asymptotic certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate import (  # noqa: E402
    B0_FLOOR,
    B1_MAGNITUDE_FLOOR,
    COLLAR_RADIUS,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    FIRST_COMPOSITION_CAP,
    INVERSE_T_CAP,
    MID_B_CAP,
    MID_CUMULANT_CAP,
    NESTED_COMPOSITION_CAP,
    NORMALIZED_CAP,
    PRECISION_BITS,
    SCALED_TARGET,
    SOURCE_FINITE_RAY,
    SOURCE_HIGH_CUMULANTS,
    SOURCE_MID_CUMULANTS,
    SOURCE_ORDER4_RAY,
    SOURCE_TOP_CUMULANTS,
    TOP_B_CAP,
    TOP_CUMULANT_CAP,
    ULTRA_B_CAP,
    ULTRA_CUMULANT_CAP,
    build_artifact,
)
from jensen_window_pf_compound_order8_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)


@dataclass(frozen=True)
class Issue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> Issue:
    return Issue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[Issue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[Issue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate"
    ):
        issues.append(issue("artifact", "kind", artifact.get("kind")))
    if "complete u>=20 ray" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "status", artifact.get("status")))

    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "inverse_t_cap": str(INVERSE_T_CAP),
        "collar_radius": COLLAR_RADIUS,
        "normalized_cap": str(NORMALIZED_CAP),
        "B0_floor": str(B0_FLOOR),
        "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
        "mid_cumulant_cap": MID_CUMULANT_CAP,
        "top_cumulant_cap": TOP_CUMULANT_CAP,
        "ultra_cumulant_cap": ULTRA_CUMULANT_CAP,
        "mid_B_cap": MID_B_CAP,
        "top_B_cap": TOP_B_CAP,
        "ultra_B_cap": ULTRA_B_CAP,
        "first_composition_cap": FIRST_COMPOSITION_CAP,
        "nested_composition_cap": NESTED_COMPOSITION_CAP,
        "scaled_target": SCALED_TARGET,
        "continuous_target": CURVATURE_CONSTANT,
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append(issue("parameters", "drift", artifact.get("parameters")))
    expected_summary = {
        "rows": 7,
        "ready_rows": 7,
        "normalized_H_box_theorems": 2,
        "stable_log_majorants": 1,
        "dimensionless_interval_theorems": 1,
        "asymptotic_curvature_theorems": 1,
        "global_first_summand_handoffs": 1,
        "open_rows": 0,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", key, summary.get(key)))

    sources = artifact.get("source_contract", {})
    if sources.get("normalized_H_boxes") != (
        "0<x_r<=1 for 2<=r<=8, x_2>=97/100, x_3>=24/25"
    ):
        issues.append(issue("sources", "H-boxes", sources))
    if sources.get("inverse_t_cap") != str(INVERSE_T_CAP):
        issues.append(issue("sources", "inverse-t", sources))
    for source in (
        SOURCE_ORDER4_RAY,
        SOURCE_MID_CUMULANTS,
        SOURCE_TOP_CUMULANTS,
        SOURCE_HIGH_CUMULANTS,
        SOURCE_FINITE_RAY,
    ):
        if not source.exists():
            issues.append(issue("sources", "missing", source))

    ratios = artifact.get("ratio_gates", {})
    if ratios.get("all_strict") is not True:
        issues.append(issue("ratios", "not-strict", ratios))
    if ratios.get("collar_radius") != COLLAR_RADIUS:
        issues.append(issue("ratios", "collar", ratios.get("collar_radius")))
    high_caps = ratios.get("high_normalized_caps", {})
    if set(high_caps) != {"9", "10", "11", "12", "13", "14"}:
        issues.append(issue("ratios", "orders", high_caps))
    for order in ("9", "10"):
        if high_caps.get(order, {}).get("target_cap") != MID_B_CAP:
            issues.append(issue("ratios", f"cap-{order}", high_caps.get(order)))
    for order in ("11", "12"):
        if high_caps.get(order, {}).get("target_cap") != TOP_B_CAP:
            issues.append(issue("ratios", f"cap-{order}", high_caps.get(order)))
    for order in ("13", "14"):
        if high_caps.get(order, {}).get("target_cap") != ULTRA_B_CAP:
            issues.append(issue("ratios", f"cap-{order}", high_caps.get(order)))

    interval = artifact.get("dimensionless_interval", {})
    decimal_fields = (
        "scaled_curvature_upper",
        "scaled_margin_lower",
        "j0_lower",
        "r0_lower",
        "s0_lower",
        "t0_lower",
        "u0_lower",
    )
    try:
        values = {key: Decimal(interval.get(key, "nan")) for key in decimal_fields}
    except Exception as exc:
        issues.append(issue("interval", "decimal", exc))
    else:
        if (
            not values["scaled_curvature_upper"].is_finite()
            or values["scaled_curvature_upper"] >= SCALED_TARGET
            or values["scaled_curvature_upper"] >= CURVATURE_CONSTANT
        ):
            issues.append(issue("interval", "scaled-upper", values))
        for key in decimal_fields[1:]:
            if values[key] <= 0:
                issues.append(issue("interval", key, values[key]))
    for key, count in (
        ("normalized_B_boxes", 13),
        ("normalized_J_boxes", 11),
        ("normalized_R_boxes", 9),
        ("normalized_S_boxes", 7),
        ("normalized_T_boxes", 5),
        ("normalized_U_boxes", 3),
    ):
        if len(interval.get(key, [])) != count:
            issues.append(issue("interval", key, len(interval.get(key, []))))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 7 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-unique", ids))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "not-ready", rows))
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if canonical != artifact:
            issues.append(issue("rebuild", "drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[Issue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous nested first-summand order-eight curvature theorem on",
        "This certificate is not a proof",
        "700001, r=11,12",
        "<1, r=13,14",
        "t+-6",
        "|b_11|,|b_12|<2",
        "R(w)=log((1-exp(-w))/w)",
        "exact partial-Bell identity",
        "0<=z<=10^-30",
        "t^2*s_1'' upper=",
        "<200<4000",
        "uniform interval theorem",
        "outputs/jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.md",
        "outputs/formal_core.md",
    )
    return [
        issue("note", "marker", marker)
        for marker in required
        if marker not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    for item in issues:
        print(f"ORDER8-NESTED-ASYMPTOTIC {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-eight nested asymptotic certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    interval = artifact["dimensionless_interval"]
    print(
        "validated order-eight nested curvature asymptotic certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['normalized_H_box_theorems']} normalized-H theorems, "
        f"{summary['dimensionless_interval_theorems']} dimensionless interval theorem, "
        f"{summary['asymptotic_curvature_theorems']} asymptotic curvature theorem, "
        f"{summary['open_rows']} open rows, "
        f"scaled upper {interval['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
