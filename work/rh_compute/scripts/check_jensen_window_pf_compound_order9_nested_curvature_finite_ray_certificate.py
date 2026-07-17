#!/usr/bin/env python3
"""Validate the order-nine finite-ray curvature certificate."""

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

import jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate as finite9  # noqa: E402


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
        "jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    finite = artifact.get("finite_ray", {})
    if finite.get("mode_range") != ["2", "20"]:
        issues.append(issue("finite", "bad-mode-range", finite.get("mode_range")))
    if finite.get("ray_blocks") != 17999 or finite.get("all_blocks_passed") is not True:
        issues.append(issue("finite", "bad-block-contract", finite.get("ray_blocks")))
    if finite.get("theorem") != (
        "w_1''(t)<=4200/t^2 for every saddle mode 2<=u<=20"
    ):
        issues.append(issue("finite", "bad-theorem", finite.get("theorem")))
    try:
        largest = Decimal(finite.get("largest_scaled_curvature_upper", "Infinity"))
        weakest = Decimal(finite.get("weakest_V_lower", "-Infinity"))
    except Exception as exc:
        issues.append(issue("finite", "bad-decimal", exc))
    else:
        if largest >= Decimal("4200"):
            issues.append(issue("finite", "curvature-failure", largest))
        if weakest <= 0:
            issues.append(issue("finite", "nonpositive-V", weakest))
    summary = artifact.get("summary", {})
    expected = {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "initial_collar_blocks": 1,
        "finite_ray_blocks": 17999,
        "finite_ray_theorems": 1,
        "open_asymptotic_rays": 1,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    return artifact, issues


def validate_rebuild(artifact: dict) -> list[RayIssue]:
    issues = []
    try:
        extension = finite9.load_cache(
            finite9.DEFAULT_EXTENSION_CACHE,
            finite9.extension_tasks(),
        )
        ray = finite9.load_cache(finite9.DEFAULT_RAY_CACHE, finite9.ray_tasks())
        canonical = finite9.build_artifact(extension, ray)
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", finite9.DEFAULT_OUT))
    return issues


def validate_note(path: Path) -> list[RayIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    markers = (
        "Status: rigorous first-summand order-nine theorem on `2<=u<=20`",
        "This is not a proof",
        "H2-H16",
        "`t+-7`",
        "`1/1000`",
        "w_1''(t)<=4200/t^2 for every saddle mode 2<=u<=20",
        "u>=20",
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
    parser.add_argument("--artifact", type=Path, default=finite9.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=finite9.DEFAULT_NOTE)
    parser.add_argument("--skip-rebuild", action="store_true")
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if artifact and not args.skip_rebuild:
        issues.extend(validate_rebuild(artifact))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-FINITE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine finite-ray certificate: {len(issues)} issues")
        return 1
    finite = artifact["finite_ray"]
    print(
        "validated order-nine finite-ray certificate: 17999 blocks, 0 issues, "
        f"largest scaled upper {finite['largest_scaled_curvature_upper']}, "
        f"weakest V lower {finite['weakest_V_lower']}, 1 open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
