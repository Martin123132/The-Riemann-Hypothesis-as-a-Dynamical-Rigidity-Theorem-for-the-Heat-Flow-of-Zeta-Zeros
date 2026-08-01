#!/usr/bin/env python3
"""Validate the order-twelve nested-curvature asymptotic-ray certificate."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate as generator  # noqa: E402


def main() -> int:
    artifact = json.loads(generator.DEFAULT_OUT.read_text(encoding="utf-8"))
    issues = []
    if artifact.get("kind") != "jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate":
        issues.append("bad artifact kind")
    rebuilt = generator.build_artifact()
    if artifact != rebuilt:
        issues.append("rebuilt artifact differs")
    interval = artifact.get("dimensionless_interval", {})
    if Decimal(interval.get("scaled_curvature_upper", "Infinity")) >= Decimal(8000):
        issues.append("curvature margin failed")
    if Fraction(interval.get("dimensionless_D9_floor", "0")) <= Fraction(499, 100):
        issues.append("D9 floor failed")
    for name in ("J", "R", "S", "T", "U", "V", "W", "X", "D9"):
        if Decimal(interval.get(f"{name}0_lower", "-Infinity")) <= 0:
            issues.append(f"nonpositive {name} floor")
    expected_summary = {
        "rows": 3,
        "ready_rows": 3,
        "asymptotic_ray_theorems": 1,
        "open_finite_ray_ranges": 0,
        "open_lower_compact_ranges": 1,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("summary changed")
    note = generator.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "rigorous first-summand order-twelve theorem on `u>=20`",
        "t^2*v_1''(t)<8000",
        "dimensionless D9 floor=",
        "finite ray `2001/1000<=u<=20`",
        "compact handoff remain open",
    ):
        if marker not in note:
            issues.append(f"note misses marker: {marker}")
    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print(
        "validated order-twelve asymptotic ray: 0 issues, scaled upper "
        f"{interval['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
