#!/usr/bin/env python3
"""Validate the order-six nested-curvature finite-ray certificate."""

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

from jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_EXTENSION_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_RAY_CACHE,
    build_artifact,
    extension_tasks,
    load_cache,
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
    if artifact.get("kind") != "jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate":
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    extension = load_cache(DEFAULT_EXTENSION_CACHE, extension_tasks())
    ray = load_cache(DEFAULT_RAY_CACHE, ray_tasks())
    expected = {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "initial_collar_blocks": 1,
        "finite_ray_blocks": len(ray_tasks()),
        "finite_ray_theorems": 1,
        "open_asymptotic_rays": 1,
    }
    for key, value in expected.items():
        if artifact.get("summary", {}).get(key) != value:
            issues.append(Issue("summary", key, str(artifact.get("summary", {}).get(key))))
    cache = artifact.get("cache", {})
    if cache.get("extension_sha256") != sha256(DEFAULT_EXTENSION_CACHE):
        issues.append(Issue("cache", "extension-hash", str(DEFAULT_EXTENSION_CACHE)))
    if cache.get("ray_sha256") != sha256(DEFAULT_RAY_CACHE):
        issues.append(Issue("cache", "ray-hash", str(DEFAULT_RAY_CACHE)))
    finite = artifact.get("finite_ray", {})
    if finite.get("mode_range") != ["2", "20"] or not finite.get("all_blocks_passed"):
        issues.append(Issue("finite", "coverage", str(finite)))
    if Decimal(str(finite.get("largest_scaled_curvature_upper", "200"))) >= 200:
        issues.append(Issue("finite", "scaled-upper", str(finite.get("largest_scaled_curvature_upper"))))
    if Decimal(str(finite.get("weakest_S_lower", "0"))) <= 0:
        issues.append(Issue("finite", "S-lower", str(finite.get("weakest_S_lower"))))
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
            "Status: rigorous first-summand order-six",
            "2<=u<=20",
            "t+-4",
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
        print(f"ORDER6-NESTED-FINITE-RAY {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-six nested finite-ray certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-six nested finite-ray certificate: "
        f"{summary['finite_ray_blocks']} ray blocks, 0 issues, "
        f"{summary['finite_ray_theorems']} theorem, "
        f"{summary['open_asymptotic_rays']} open asymptotic ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
