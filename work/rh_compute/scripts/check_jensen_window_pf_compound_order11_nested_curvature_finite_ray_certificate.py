#!/usr/bin/env python3
"""Validate the order-eleven nested-curvature finite-ray certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate as certificate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=certificate.DEFAULT_CACHE)
    parser.add_argument("--artifact", type=Path, default=certificate.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=certificate.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    tasks = certificate.ray_tasks()
    try:
        records = certificate.load_cache(args.cache, tasks)
        actual = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = certificate.build_artifact(records)
    except (OSError, json.JSONDecodeError, RuntimeError, KeyError, ValueError) as exc:
        print(f"order-eleven finite-ray certificate: source failure: {exc}")
        return 1
    if len(records) != 17999:
        issues.append(f"expected 17999 blocks, found {len(records)}")
    if actual != expected:
        issues.append("artifact differs from deterministic cache/source rebuild")
    finite = actual.get("finite_ray", {})
    if (
        actual.get("kind")
        != "jensen_window_pf_compound_order11_nested_curvature_finite_ray_certificate"
        or actual.get("theorem") != certificate.THEOREM
        or finite.get("mode_range") != ["2001/1000", "20"]
        or finite.get("blocks") != 17999
        or finite.get("all_blocks_passed") is not True
    ):
        issues.append("identity or complete finite-ray summary changed")
    try:
        if Decimal(finite["largest_scaled_curvature_upper"]) >= Decimal(6000):
            issues.append("curvature cap failed")
        if Decimal(finite["smallest_margin_lower"]) <= 0:
            issues.append("global margin is nonpositive")
        if Decimal(finite["weakest_X_lower"]) <= 0:
            issues.append("global X floor is nonpositive")
    except (KeyError, ValueError):
        issues.append("numeric summary parse failure")
    for index in (0, len(tasks) // 2, len(tasks) - 1):
        rebuilt = certificate.ray_task(tasks[index])
        if rebuilt != records[index]:
            issues.append(f"independent block rebuild mismatch at index {index}")
    if not args.note.exists():
        issues.append("missing note")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in ("17999", "2001/1000<=u<=20", "This is not a proof of RH"):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-eleven finite-ray certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-eleven finite ray: 17999 blocks, "
        "3 exact rebuilds, 0 issues, largest scaled upper "
        f"{finite['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
