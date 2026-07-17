#!/usr/bin/env python3
"""Validate the order-nine asymptotic-ray curvature certificate."""

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

import jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate as ray9  # noqa: E402


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
    issues = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    interval = artifact.get("dimensionless_interval", {})
    try:
        scaled = Decimal(interval.get("scaled_curvature_upper", "Infinity"))
        margin = Decimal(interval.get("scaled_margin_lower", "-Infinity"))
    except Exception as exc:
        issues.append(issue("interval", "bad-decimal", exc))
    else:
        if scaled >= Decimal("500"):
            issues.append(issue("interval", "scaled-failure", scaled))
        if margin <= 0:
            issues.append(issue("interval", "nonpositive-margin", margin))
    for name in ("J", "R", "S", "T", "U", "V"):
        try:
            lower = Decimal(interval.get(f"{name}0_lower", "-Infinity"))
        except Exception as exc:
            issues.append(issue("interval", f"bad-{name}", exc))
        else:
            if lower <= 0:
                issues.append(issue("interval", f"nonpositive-{name}", lower))
    gates = artifact.get("ratio_gates", {})
    if gates.get("all_strict") is not True or gates.get("collar_radius") != 7:
        issues.append(issue("gates", "bad-collar-contract", gates.get("collar_radius")))
    if len(gates.get("high_normalized_caps", {})) != 8:
        issues.append(issue("gates", "bad-high-cap-count", None))
    expected_summary = {
        "rows": 4,
        "ready_rows": 4,
        "asymptotic_ray_theorems": 1,
        "global_above_5700_compositions": 1,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    return artifact, issues


def validate_rebuild(artifact: dict) -> list[RayIssue]:
    try:
        canonical = ray9.build_artifact()
    except Exception as exc:
        return [issue("rebuild", "exception", exc)]
    if canonical != artifact:
        return [issue("rebuild", "artifact-drift", ray9.DEFAULT_OUT)]
    return []


def validate_note(path: Path) -> list[RayIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    markers = (
        "Status: rigorous first-summand order-nine theorem on `u>=20`",
        "This is not a proof",
        "t^2*w_1''(t)<500<4200 for every mode u>=20",
        "B,J,R,S,T,U,V>0",
        "t>=5700",
        "1250<=t<=5700",
    )
    issues = []
    for marker in markers:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in ("therefore rh", "pf-infinity follows", "entry is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ray9.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=ray9.DEFAULT_NOTE)
    parser.add_argument("--skip-rebuild", action="store_true")
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if artifact and not args.skip_rebuild:
        issues.extend(validate_rebuild(artifact))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-ASYMPTOTIC {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine asymptotic-ray certificate: {len(issues)} issues")
        return 1
    interval = artifact["dimensionless_interval"]
    print(
        "validated order-nine asymptotic-ray certificate: 0 issues, "
        f"scaled upper {interval['scaled_curvature_upper']}<500<4200, "
        "6 positive stable coordinates, 1 global above-5700 composition"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
