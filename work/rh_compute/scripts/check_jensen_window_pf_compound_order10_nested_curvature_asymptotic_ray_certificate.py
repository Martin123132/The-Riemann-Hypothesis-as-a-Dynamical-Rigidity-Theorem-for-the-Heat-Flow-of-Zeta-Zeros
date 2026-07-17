#!/usr/bin/env python3
"""Validate the order-ten nested-curvature asymptotic-ray certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class RayIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> RayIssue:
    return RayIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[RayIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[RayIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    interval = artifact.get("dimensionless_interval", {})
    if Decimal(interval.get("scaled_curvature_upper", "Infinity")) >= 1000:
        issues.append(
            issue("interval", "curvature-margin", interval.get("scaled_curvature_upper"))
        )
    if Fraction(interval.get("dimensionless_W_floor", "0")) <= Fraction(399, 100):
        issues.append(
            issue("interval", "bad-W-floor", interval.get("dimensionless_W_floor"))
        )
    for name in ("J", "R", "S", "T", "U", "V", "W"):
        if Decimal(interval.get(f"{name}0_lower", "-Infinity")) <= 0:
            issues.append(issue("interval", f"nonpositive-{name}", interval.get(f"{name}0_lower")))
    summary = artifact.get("summary", {})
    expected = {
        "rows": 4,
        "ready_rows": 4,
        "asymptotic_ray_theorems": 1,
        "global_above_mode_2001_1000_compositions": 1,
        "open_compact_ranges": 1,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        canonical = build_artifact()
    except Exception as exc:
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
        "Status: rigorous first-summand order-ten theorem on `u>=20`",
        "This is not a proof",
        "t^2*z_1''(t)<1000<5500",
        "dimensionless W floor=",
        "u>=2001/1000",
        "compact handoff remains open",
    )
    return [
        issue("note", "missing-marker", marker)
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
    if issues:
        for finding in issues:
            print(
                f"ORDER10-ASYMPTOTIC {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-ten asymptotic ray: {len(issues)} issues")
        return 1
    print(
        "validated order-ten asymptotic ray: 0 issues, scaled upper "
        f"{artifact['dimensionless_interval']['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
