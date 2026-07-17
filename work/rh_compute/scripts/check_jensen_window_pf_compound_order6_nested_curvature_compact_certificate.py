#!/usr/bin/env python3
"""Validate the order-six nested-curvature compact certificate."""

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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
from jensen_window_pf_compound_order6_nested_curvature_compact_certificate import (  # noqa: E402
    DEFAULT_EXTENSION_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_RIGHT_COLLAR,
    append_right_collar,
    build_artifact,
    build_right_collar,
    deterministic_tasks,
    load_extension_cache,
    sha256,
)


@dataclass(frozen=True)
class Issue:
    section: str
    code: str
    detail: str


def validate(
    artifact_path: Path,
    note_path: Path,
    cache_path: Path,
    right_collar_path: Path,
) -> tuple[dict, list[Issue]]:
    if not artifact_path.exists():
        return {}, [Issue("artifact", "missing", str(artifact_path))]
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    issues: list[Issue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order6_nested_curvature_compact_certificate":
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    summary = artifact.get("summary", {})
    for key, value in {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "all_blocks_passed": True,
        "compact_curvature_theorems": 1,
        "open_ray_targets": 1,
    }.items():
        if summary.get(key) != value:
            issues.append(Issue("summary", key, str(summary.get(key))))
    tasks = deterministic_tasks()
    extension = load_extension_cache(cache_path, tasks)
    if len(extension) != len(tasks):
        issues.append(Issue("cache", "row-count", f"{len(extension)}/{len(tasks)}"))
        return artifact, issues
    contract = artifact.get("source_contract", {})
    if contract.get("extension_cache_sha256") != sha256(cache_path):
        issues.append(Issue("cache", "extension-hash", str(cache_path)))
    if contract.get("base_cache_sha256") != sha256(order4_compact.DEFAULT_CACHE):
        issues.append(Issue("cache", "base-hash", str(order4_compact.DEFAULT_CACHE)))
    if contract.get("right_collar_sha256") != sha256(right_collar_path):
        issues.append(Issue("cache", "right-collar-hash", str(right_collar_path)))
    compact = artifact.get("compact", {})
    if compact.get("certified_t_range") != "321<=t<=V'(2)":
        issues.append(Issue("compact", "range", str(compact.get("certified_t_range"))))
    if Decimal(str(compact.get("largest_scaled_curvature_upper", "200"))) >= 200:
        issues.append(Issue("compact", "scaled-upper", str(compact.get("largest_scaled_curvature_upper"))))
    if Decimal(str(compact.get("weakest_S_lower", "0"))) <= 0:
        issues.append(Issue("compact", "S-lower", str(compact.get("weakest_S_lower"))))
    if not compact.get("all_blocks_passed"):
        issues.append(Issue("compact", "blocks", "not all passed"))
    try:
        base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
        collar = build_right_collar(right_collar_path)
        base, extension = append_right_collar(base, extension, collar)
        canonical = build_artifact(
            base, extension, cache_path, right_collar_path
        )
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
            "H^(9),H^(10)",
            "t+-4",
            "p_1''(t)<=200/t^2",
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
    parser.add_argument("--extension-cache", type=Path, default=DEFAULT_EXTENSION_CACHE)
    parser.add_argument("--right-collar", type=Path, default=DEFAULT_RIGHT_COLLAR)
    args = parser.parse_args()
    artifact, issues = validate(
        args.artifact, args.note, args.extension_cache, args.right_collar
    )
    for item in issues:
        print(f"ORDER6-NESTED-COMPACT {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-six nested compact certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-six nested compact certificate: "
        f"{summary['accepted_blocks']} blocks, 0 issues, "
        f"{summary['compact_curvature_theorems']} compact theorem, "
        f"{summary['open_ray_targets']} open ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
