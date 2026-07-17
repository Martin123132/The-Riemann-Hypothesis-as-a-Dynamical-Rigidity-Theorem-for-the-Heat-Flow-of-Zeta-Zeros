#!/usr/bin/env python3
"""Validate the asymptotic nested-curvature certificate for order five."""

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

from jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate import (  # noqa: E402
    B0_FLOOR,
    B1_MAGNITUDE_FLOOR,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    FIRST_COMPOSITION_CAP,
    INVERSE_T_CAP,
    NESTED_COMPOSITION_CAP,
    NORMALIZED_CAP,
    SCALED_TARGET,
    SOURCE_FINITE_RAY,
    SOURCE_ORDER4_RAY,
    build_artifact,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)


@dataclass(frozen=True)
class RayIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> RayIssue:
    return RayIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[RayIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[RayIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "complete u>=20 ray" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_parameters = {
        "precision_bits": 256,
        "inverse_t_cap": str(INVERSE_T_CAP),
        "collar_radius": 3,
        "normalized_cap": str(NORMALIZED_CAP),
        "B0_floor": str(B0_FLOOR),
        "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
        "first_composition_cap": FIRST_COMPOSITION_CAP,
        "nested_composition_cap": NESTED_COMPOSITION_CAP,
        "scaled_target": SCALED_TARGET,
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append(issue("parameters", "drift", artifact.get("parameters")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 6,
        "ready_rows": 6,
        "normalized_H_box_theorems": 1,
        "stable_log_majorants": 1,
        "dimensionless_interval_theorems": 1,
        "asymptotic_curvature_theorems": 1,
        "open_rows": 0,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    sources = artifact.get("source_contract", {})
    if sources.get("normalized_H_boxes") != (
        "0<x_r<=1 for 2<=r<=8, x_2>=97/100, x_3>=24/25"
    ):
        issues.append(issue("sources", "bad-H-boxes", sources))
    if sources.get("inverse_t_cap") != str(INVERSE_T_CAP):
        issues.append(issue("sources", "bad-inverse-t-cap", sources))
    for source in (SOURCE_ORDER4_RAY, SOURCE_FINITE_RAY):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    ratios = artifact.get("ratio_gates", {})
    if ratios.get("all_strict") is not True:
        issues.append(issue("ratios", "not-strict", ratios))
    if ratios.get("upper_ratio_cap") != str(NORMALIZED_CAP):
        issues.append(issue("ratios", "bad-upper-cap", ratios))
    if ratios.get("lower_ratio_floor") != "999/1000":
        issues.append(issue("ratios", "bad-lower-floor", ratios))

    interval = artifact.get("dimensionless_interval", {})
    try:
        upper = Decimal(interval.get("scaled_curvature_upper", "nan"))
        margin = Decimal(interval.get("scaled_margin_lower", "nan"))
        j0 = Decimal(interval.get("j0_lower", "nan"))
        r0 = Decimal(interval.get("r0_lower", "nan"))
    except Exception as exc:
        issues.append(issue("interval", "bad-decimal", exc))
    else:
        if not upper.is_finite() or upper >= SCALED_TARGET:
            issues.append(issue("interval", "scaled-target-failed", upper))
        if margin <= 0 or j0 <= 0 or r0 <= 0:
            issues.append(issue("interval", "nonpositive-margin-or-floor", interval))
        if upper >= CURVATURE_CONSTANT:
            issues.append(issue("interval", "continuous-target-failed", upper))
    if interval.get("scaled_target") != SCALED_TARGET:
        issues.append(issue("interval", "bad-scaled-target", interval))
    if interval.get("first_composition_coefficient_cap") != FIRST_COMPOSITION_CAP:
        issues.append(issue("interval", "bad-first-cap", interval))
    if interval.get("nested_composition_coefficient_cap") != NESTED_COMPOSITION_CAP:
        issues.append(issue("interval", "bad-nested-cap", interval))
    if len(interval.get("normalized_B_boxes", [])) != 7:
        issues.append(issue("interval", "bad-B-box-count", interval))
    if len(interval.get("normalized_J_boxes", [])) != 5:
        issues.append(issue("interval", "bad-J-box-count", interval))
    if len(interval.get("normalized_R_boxes", [])) != 3:
        issues.append(issue("interval", "bad-R-box-count", interval))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 6 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "non-ready-row", rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[RayIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous nested first-summand order-five curvature theorem on",
        "This is not by itself a proof",
        "0<x_r<=1 for 2<=r<=8",
        "1/t<=10^-30",
        "within `t+-3`",
        "b_0 in [969/1000,1001/1000]",
        "b_1 in [-1001/1000,-959/1000]",
        "R_0(w)=log((1-exp(-w))/w)",
        "exact partial-Bell identity",
        "C=2",
        "C=100",
        "0<=z<=10^-30",
        "t^2*q_1''(t)<10<60 for every mode u>=20",
        "uniform interval theorem",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md",
        "outputs/formal_core.md",
    )
    issues: list[RayIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "full order-five entry is complete",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER5-NESTED-ASYMPTOTIC {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five nested curvature asymptotic certificate: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    interval = artifact["dimensionless_interval"]
    print(
        "validated order-five nested curvature asymptotic certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['normalized_H_box_theorems']} normalized-H theorem, "
        f"{summary['stable_log_majorants']} stable-log majorant, "
        f"{summary['dimensionless_interval_theorems']} dimensionless interval theorem, "
        f"{summary['asymptotic_curvature_theorems']} asymptotic curvature theorem, "
        f"{summary['open_rows']} open rows, "
        f"scaled upper {interval['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
