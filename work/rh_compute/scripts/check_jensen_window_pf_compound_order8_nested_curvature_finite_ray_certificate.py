#!/usr/bin/env python3
"""Validate the order-eight nested-curvature finite-ray certificate."""

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

from jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate import (  # noqa: E402
    CURVATURE_CONSTANT,
    DEFAULT_EXTENSION_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_RAY_CACHE,
    build_artifact,
    extension_tasks,
    load_cache,
    ray_task,
    ray_tasks,
    sha256,
)


@dataclass(frozen=True)
class Issue:
    section: str
    code: str
    detail: str


def validate(artifact_path: Path, note_path: Path) -> tuple[dict, list[Issue]]:
    if not artifact_path.exists():
        return {}, [Issue("artifact", "missing", str(artifact_path))]
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    issues: list[Issue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate"
    ):
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))

    extension_task_rows = extension_tasks()
    ray_task_rows = ray_tasks()
    extension = load_cache(DEFAULT_EXTENSION_CACHE, extension_task_rows)
    ray = load_cache(DEFAULT_RAY_CACHE, ray_task_rows)
    expected = {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "initial_collar_blocks": 1,
        "finite_ray_blocks": len(ray_task_rows),
        "finite_ray_theorems": 1,
        "open_asymptotic_rays": 1,
    }
    summary = artifact.get("summary", {})
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(Issue("summary", key, str(summary.get(key))))
    if len(extension) != len(extension_task_rows):
        issues.append(
            Issue(
                "cache",
                "extension-rows",
                f"{len(extension)}/{len(extension_task_rows)}",
            )
        )
    if len(ray) != len(ray_task_rows):
        issues.append(Issue("cache", "ray-rows", f"{len(ray)}/{len(ray_task_rows)}"))
    cache = artifact.get("cache", {})
    if cache.get("extension_sha256") != sha256(DEFAULT_EXTENSION_CACHE):
        issues.append(Issue("cache", "extension-hash", str(DEFAULT_EXTENSION_CACHE)))
    if cache.get("ray_sha256") != sha256(DEFAULT_RAY_CACHE):
        issues.append(Issue("cache", "ray-hash", str(DEFAULT_RAY_CACHE)))

    finite = artifact.get("finite_ray", {})
    if finite.get("mode_range") != ["2", "20"]:
        issues.append(Issue("finite", "range", str(finite.get("mode_range"))))
    if finite.get("all_blocks_passed") is not True:
        issues.append(Issue("finite", "coverage", "not all blocks passed"))
    if finite.get("ray_blocks") != len(ray):
        issues.append(Issue("finite", "block-count", str(finite.get("ray_blocks"))))
    if Decimal(
        str(finite.get("largest_scaled_curvature_upper", CURVATURE_CONSTANT))
    ) >= CURVATURE_CONSTANT:
        issues.append(
            Issue(
                "finite",
                "scaled-upper",
                str(finite.get("largest_scaled_curvature_upper")),
            )
        )
    if Decimal(str(finite.get("weakest_U_lower", "0"))) <= 0:
        issues.append(Issue("finite", "U-lower", str(finite.get("weakest_U_lower"))))
    collar = finite.get("initial_collar", {})
    for side in ("left_t_collar", "right_t_collar"):
        ball = collar.get("diagnostics", {}).get(side, "")
        if not ball.startswith("["):
            issues.append(Issue("collar", side, str(ball)))
    for index, row in enumerate(ray):
        if row.get("passed") is not True:
            issues.append(Issue("ray", "failed", str(index)))
        if Decimal(str(row.get("U_lower", "0"))) <= 0:
            issues.append(Issue("ray", "U-lower", str(index)))
        if Decimal(
            str(row.get("scaled_curvature_upper", CURVATURE_CONSTANT))
        ) >= CURVATURE_CONSTANT:
            issues.append(Issue("ray", "scaled-upper", str(index)))
        caps = row.get("diagnostics", {}).get("cumulant_caps", {})
        if caps.get("13") != 1 or caps.get("14") != 1:
            issues.append(Issue("ray", "ultra-cap", str(index)))
    for index, (previous, current) in enumerate(zip(ray, ray[1:])):
        if previous.get("mode_right") != current.get("mode_left"):
            issues.append(Issue("ray", "mode-gap", str(index)))

    for index in (0, len(ray_task_rows) // 2, len(ray_task_rows) - 1):
        if ray_task(ray_task_rows[index]) != ray[index]:
            issues.append(Issue("recompute", "representative-drift", str(index)))
    try:
        canonical = build_artifact(extension, ray)
    except Exception as exc:  # pragma: no cover
        issues.append(Issue("rebuild", "exception", repr(exc)))
    else:
        if canonical != artifact:
            issues.append(Issue("rebuild", "drift", str(artifact_path)))

    if not note_path.exists():
        issues.append(Issue("note", "missing", str(note_path)))
    else:
        text = note_path.read_text(encoding="utf-8")
        for marker in (
            "Status: rigorous first-summand order-eight",
            "2<=u<=20",
            "H^(2),...,H^(14)",
            "t+-6",
            "width `1/1000`",
            "This is not a proof",
            "outputs/formal_core.md",
        ):
            if marker not in text:
                issues.append(Issue("note", "marker", marker))
    return artifact, issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate(args.artifact, args.note)
    for item in issues:
        print(f"ORDER8-NESTED-FINITE-RAY {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-eight nested finite-ray certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight nested finite-ray certificate: "
        f"{summary['finite_ray_blocks']} ray blocks, 0 issues, "
        f"{summary['finite_ray_theorems']} theorem, "
        f"{summary['open_asymptotic_rays']} open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
