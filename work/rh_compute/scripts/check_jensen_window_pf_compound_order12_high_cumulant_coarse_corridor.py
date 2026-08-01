#!/usr/bin/env python3
"""Validate the order-twelve high-cumulant coarse corridor."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_high_cumulant_coarse_corridor as generator  # noqa: E402


def main() -> int:
    artifact = json.loads(generator.DEFAULT_OUT.read_text(encoding="utf-8"))
    issues = []
    if artifact.get("kind") != "jensen_window_pf_compound_order12_high_cumulant_coarse_corridor":
        issues.append("bad artifact kind")
    rebuilt = generator.build_artifact()
    if artifact != rebuilt:
        issues.append("rebuilt artifact differs")
    exact = artifact.get("exact", {})
    if exact.get("formal_orders") != [21, 22]:
        issues.append("formal orders changed")
    if exact.get("cauchy_factor") != 462:
        issues.append("Cauchy factor changed")
    if Fraction(exact.get("finite_scaled_residual_bound", "1")) != Fraction(133, 1875):
        issues.append("finite residual changed")
    if exact.get("ray_scaled_residual_bound") != "7693/150000/u":
        issues.append("ray residual changed")
    expected_summary = {
        "rows": 4,
        "ready_rows": 4,
        "formal_orders": 2,
        "formal_terms": 0,
        "cauchy_extensions": 1,
        "global_coarse_corridors": 2,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("summary changed")
    note = generator.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "twenty-first- and twenty-second-cumulant",
        "scaled kappa_21^[10]=scaled kappa_22^[10]=0",
        "22*21=462",
        "finite residual < 133/1875",
        "ray residual < 7693/150000/u",
        "No numerical",
    ):
        if marker not in note:
            issues.append(f"note misses marker: {marker}")
    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print(
        "validated order-twelve high-cumulant corridor: "
        "0 formal terms, 2 exact corridors, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
