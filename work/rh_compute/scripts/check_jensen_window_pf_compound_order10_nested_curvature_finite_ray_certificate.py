#!/usr/bin/env python3
"""Validate the order-ten nested-curvature finite-ray certificate."""

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

from jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_RAY_CACHE,
    build_artifact,
    order9,
    ray_tasks,
    sha256,
)


@dataclass(frozen=True)
class RayIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> RayIssue:
    return RayIssue(section, code, str(detail))


def validate_artifact(
    artifact_path: Path,
    cache_path: Path,
) -> tuple[dict, list[RayIssue]]:
    if not artifact_path.exists():
        return {}, [issue("artifact", "missing-file", artifact_path)]
    if not cache_path.exists():
        return {}, [issue("cache", "missing-file", cache_path)]
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    issues: list[RayIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    finite = artifact.get("finite_ray", {})
    if finite.get("mode_range") != ["2001/1000", "20"]:
        issues.append(issue("finite", "bad-mode-range", finite.get("mode_range")))
    if finite.get("theorem") != (
        "z_1''(t)<=5500/t^2 for every saddle mode 2001/1000<=u<=20"
    ):
        issues.append(issue("finite", "bad-theorem", finite.get("theorem")))
    tasks = ray_tasks()
    try:
        records = order9.order7_ray.load_cache(cache_path, tasks)
    except Exception as exc:
        return artifact, [*issues, issue("cache", "load-exception", exc)]
    if len(records) != len(tasks):
        issues.append(issue("cache", "incomplete", f"{len(records)}/{len(tasks)}"))
    for index, record in enumerate(records):
        if record.get("passed") is not True:
            issues.append(issue("cache", "failed-row", index))
            break
        if Decimal(record.get("scaled_curvature_upper", "Infinity")) >= 5500:
            issues.append(
                issue(
                    "cache",
                    "curvature-margin",
                    f"{index}:{record.get('scaled_curvature_upper')}",
                )
            )
            break
        for name in ("J", "R", "S", "T", "U", "V", "W"):
            if Decimal(record.get(f"{name}_lower", "-Infinity")) <= 0:
                issues.append(issue("cache", f"nonpositive-{name}", index))
                break
        if issues:
            break
    cache = artifact.get("cache", {})
    if cache.get("ray_rows") != len(tasks):
        issues.append(issue("artifact", "bad-ray-rows", cache.get("ray_rows")))
    actual_hash = sha256(cache_path)
    if cache.get("ray_sha256") != actual_hash:
        issues.append(issue("artifact", "bad-ray-hash", cache.get("ray_sha256")))
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 4,
        "ready_rows": 2,
        "open_rows": 2,
        "initial_collar_blocks": 0,
        "finite_ray_blocks": len(tasks),
        "finite_ray_theorems": 1,
        "open_asymptotic_rays": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        canonical = build_artifact(records)
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", artifact_path))
    return artifact, issues


def validate_note(path: Path) -> list[RayIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous first-summand order-ten theorem",
        "`2001/1000<=u<=20`",
        "This is not a proof",
        "H2-H18 corridors",
        "ray blocks=17999",
        "compact handoff",
        "u>=20",
    )
    return [
        issue("note", "missing-marker", marker)
        for marker in required
        if marker not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_RAY_CACHE)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact, args.cache)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER10-FINITE-RAY {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-ten finite ray: {len(issues)} issues")
        return 1
    finite = artifact["finite_ray"]
    print(
        "validated order-ten finite ray: "
        f"{finite['ray_blocks']} blocks, 0 issues, scaled upper "
        f"{finite['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
