#!/usr/bin/env python3
"""Validate the order-eight compact nested-curvature certificate."""

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

from jensen_window_pf_compound_order8_nested_curvature_compact_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
    build_right_collar,
    deterministic_tasks,
    append_right_collar,
)
import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7_compact  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_certificate as compact8  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as h13_h14  # noqa: E402


@dataclass(frozen=True)
class CompactIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CompactIssue:
    return CompactIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[CompactIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[CompactIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order8_nested_curvature_compact_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "999<=t<=V'(2)" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    compact = artifact.get("compact", {})
    if compact.get("all_blocks_passed") is not True:
        issues.append(issue("compact", "not-all-passed", compact.get("all_blocks_passed")))
    blocks = compact.get("blocks", [])
    if compact.get("block_count") != len(blocks) or not blocks:
        issues.append(issue("compact", "bad-block-count", compact.get("block_count")))
    if compact.get("theorem") != "s_1''(t)<=4000/t^2 for every real 999<=t<=V'(2)":
        issues.append(issue("compact", "bad-theorem", compact.get("theorem")))
    try:
        largest = Decimal(compact.get("largest_scaled_curvature_upper", "Infinity"))
        weakest = Decimal(compact.get("weakest_U_lower", "-Infinity"))
    except Exception as exc:
        issues.append(issue("compact", "bad-decimal", exc))
    else:
        if largest >= Decimal("4000"):
            issues.append(issue("compact", "curvature-failure", largest))
        if weakest <= 0:
            issues.append(issue("compact", "nonpositive-U", weakest))
    for block in blocks:
        if block.get("passed") is not True:
            issues.append(issue("compact", "failed-block", block.get("central_mode")))
            break
    for previous, current in zip(blocks, blocks[1:]):
        if previous.get("central_mode", [None, None])[1] != current.get("central_mode", [None])[0]:
            issues.append(issue("compact", "mode-gap", current.get("central_mode")))
            break
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "compact_blocks": len(blocks),
        "compact_curvature_theorems": 1,
        "open_finite_rays": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    return artifact, issues


def validate_rebuild(artifact: dict) -> list[CompactIssue]:
    issues = []
    tasks = deterministic_tasks()
    try:
        base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
        high = order6_compact.load_extension_cache(order6_compact.DEFAULT_EXTENSION_CACHE, tasks)
        top = order7_compact.load_extension_cache(order7_compact.DEFAULT_EXTENSION_CACHE, tasks)
        ultra = h13_h14.load_cache(h13_h14.DEFAULT_CACHE, tasks)
        collar = build_right_collar(compact8.DEFAULT_RIGHT_COLLAR)
        base, high, top, ultra = append_right_collar(base, high, top, ultra, collar)
        canonical = build_artifact(
            base,
            high,
            top,
            ultra,
            h13_h14.DEFAULT_CACHE,
            compact8.DEFAULT_RIGHT_COLLAR,
        )
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", DEFAULT_OUT))
    return issues


def validate_note(path: Path) -> list[CompactIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous first-summand curvature theorem on `999<=t<=V'(2)`",
        "This is not a proof",
        "H2-H14",
        "strict t+-6",
        "B,J,R,S,T,U>0",
        "s_1''(t)<=4000/t^2 for every real 999<=t<=V'(2)",
        "saddle mode 2<=u<=20",
        "outputs/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-eight entry is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--skip-rebuild", action="store_true")
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if artifact and not args.skip_rebuild:
        issues.extend(validate_rebuild(artifact))
    if issues:
        for finding in issues:
            print(
                f"ORDER8-COMPACT {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-eight nested compact certificate: {len(issues)} issues")
        return 1
    compact = artifact["compact"]
    print(
        "validated order-eight nested compact certificate: "
        f"{compact['block_count']} blocks, 0 issues, "
        f"largest scaled upper {compact['largest_scaled_curvature_upper']}, "
        f"weakest U lower {compact['weakest_U_lower']}, 1 open finite ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
