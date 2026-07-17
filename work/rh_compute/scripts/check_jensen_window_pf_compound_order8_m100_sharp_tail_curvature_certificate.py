#!/usr/bin/env python3
"""Validate the sharp order-eight tail curvature certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7_compact  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_certificate as order8_compact  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as h13_h14_cache  # noqa: E402
from jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class SharpIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> SharpIssue:
    return SharpIssue(section, code, str(detail))


def rebuild_artifact(right_collar: Path) -> dict:
    tasks = order8_compact.deterministic_tasks()
    base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
    high = order6_compact.load_extension_cache(
        order6_compact.DEFAULT_EXTENSION_CACHE, tasks
    )
    top = order7_compact.load_extension_cache(
        order7_compact.DEFAULT_EXTENSION_CACHE, tasks
    )
    ultra = h13_h14_cache.load_cache(h13_h14_cache.DEFAULT_CACHE, tasks)
    collar = order8_compact.build_right_collar(right_collar)
    base, high, top, ultra = order8_compact.append_right_collar(
        base, high, top, ultra, collar
    )
    return build_artifact(
        base,
        high,
        top,
        ultra,
        h13_h14_cache.DEFAULT_CACHE,
        right_collar,
    )


def validate_artifact(
    path: Path, right_collar: Path
) -> tuple[dict, list[SharpIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[SharpIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "sharp order-eight tail" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 6,
        "ready_rows": 6,
        "open_rows": 0,
        "compact_blocks": 95,
        "sharp_compact_theorems": 1,
        "global_sharp_curvature_theorems": 1,
        "sharp_full_ceiling_theorems": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    compact = artifact.get("compact", {})
    if compact.get("all_blocks_passed") is not True:
        issues.append(issue("compact", "failed-blocks", compact))
    try:
        largest = Decimal(compact.get("largest_scaled_curvature_upper", "nan"))
    except Exception as exc:
        issues.append(issue("compact", "bad-largest", exc))
    else:
        if not (Decimal("2498") < largest < Decimal("2500")):
            issues.append(issue("compact", "bad-largest", largest))
    exact = artifact.get("exact", {})
    required = {
        "global_first_curvature": "s_1''(t)<=2500/t^2 for every real t>=1249",
        "tent_transfer": "W_k^(1)<=2500*[-log(1-1/k^2)]<2501/k^2, k>=1250",
        "sharp_full_ceiling": "W_k<2501/k^2+190/k^2=2691/k^2, k>=1250",
    }
    for key, expected in required.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    rows = artifact.get("rows", [])
    if len(rows) != 6 or len({row.get("id") for row in rows}) != 6:
        issues.append(issue("rows", "bad-rows", rows))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "bad-readiness", rows))
    try:
        canonical = rebuild_artifact(right_collar)
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[SharpIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous order-eight sharp tail theorem",
        "This is not a proof of order-nine entry",
        "all 95 sharp compact blocks pass",
        "<2500",
        "s_1''(t)<=2500/t^2 for every real t>=1249",
        "W_k<2501/k^2+190/k^2=2691/k^2, k>=1250",
        "outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-nine entry is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument(
        "--right-collar", type=Path, default=order8_compact.DEFAULT_RIGHT_COLLAR
    )
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact, args.right_collar)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER8-SHARP-TAIL {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-eight sharp tail curvature: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight sharp tail curvature: "
        f"{summary['compact_blocks']} compact blocks, 0 issues, "
        f"{summary['global_sharp_curvature_theorems']} global sharp theorem, "
        f"{summary['sharp_full_ceiling_theorems']} sharp full-kernel ceiling"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
