#!/usr/bin/env python3
"""Validate the order-nine compact nested-curvature certificate."""

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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as order8  # noqa: E402
import jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_cache as order9  # noqa: E402
import jensen_window_pf_compound_order9_nested_curvature_compact_certificate as compact9  # noqa: E402


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
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order9_nested_curvature_compact_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "5700<=t<=V'(2)" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    compact = artifact.get("compact", {})
    blocks = compact.get("blocks", [])
    if compact.get("all_blocks_passed") is not True:
        issues.append(issue("compact", "not-all-passed", None))
    if compact.get("block_count") != len(blocks) or not blocks:
        issues.append(issue("compact", "bad-block-count", len(blocks)))
    if compact.get("theorem") != (
        "w_1''(t)<=4200/t^2 for every real 5700<=t<=V'(2)"
    ):
        issues.append(issue("compact", "bad-theorem", compact.get("theorem")))
    try:
        largest = Decimal(compact.get("largest_scaled_curvature_upper", "Infinity"))
        weakest = Decimal(compact.get("weakest_V_lower", "-Infinity"))
    except Exception as exc:
        issues.append(issue("compact", "bad-decimal", exc))
    else:
        if largest >= Decimal("4200"):
            issues.append(issue("compact", "curvature-failure", largest))
        if weakest <= 0:
            issues.append(issue("compact", "nonpositive-V", weakest))
    for block in blocks:
        if block.get("passed") is not True:
            issues.append(issue("compact", "failed-block", block.get("central_mode")))
            break
    for previous, current in zip(blocks, blocks[1:]):
        if previous.get("central_mode", [None, None])[1] != current.get(
            "central_mode", [None]
        )[0]:
            issues.append(issue("compact", "mode-gap", current.get("central_mode")))
            break
    if blocks and blocks[-1].get("central_mode", [None, None])[1] != "2":
        issues.append(issue("compact", "bad-right-end", blocks[-1].get("central_mode")))
    expected_summary = {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "compact_blocks": len(blocks),
        "compact_curvature_theorems": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    return artifact, issues


def validate_rebuild(artifact: dict) -> list[CompactIssue]:
    issues = []
    try:
        tasks = compact9.deterministic_tasks()
        caches = [
            order4.load_cache(order4.DEFAULT_CACHE, tasks),
            order6.load_extension_cache(order6.DEFAULT_EXTENSION_CACHE, tasks),
            order7.load_extension_cache(order7.DEFAULT_EXTENSION_CACHE, tasks),
            order8.load_cache(order8.DEFAULT_CACHE, tasks),
            order9.load_cache(order9.DEFAULT_CACHE, tasks),
        ]
        collar = compact9.build_right_collar(compact9.DEFAULT_RIGHT_COLLAR)
        caches = compact9.append_right_collar(caches, collar)
        canonical = compact9.build_artifact(caches, compact9.DEFAULT_RIGHT_COLLAR)
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", compact9.DEFAULT_OUT))
    return issues


def validate_note(path: Path) -> list[CompactIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous first-summand curvature theorem on `5700<=t<=V'(2)`",
        "This is not a proof",
        "H2-H16",
        "`t+-7`",
        "B,J,R,S,T,U,V>0",
        "w_1''(t)<=4200/t^2 for every real 5700<=t<=V'(2)",
        "1250<=t<=5700",
        "u>=2",
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
    parser.add_argument("--artifact", type=Path, default=compact9.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=compact9.DEFAULT_NOTE)
    parser.add_argument("--skip-rebuild", action="store_true")
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if artifact and not args.skip_rebuild:
        issues.extend(validate_rebuild(artifact))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-COMPACT {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine nested compact certificate: {len(issues)} issues")
        return 1
    compact = artifact["compact"]
    print(
        "validated order-nine nested compact certificate: "
        f"{compact['block_count']} blocks, 0 issues, largest scaled upper "
        f"{compact['largest_scaled_curvature_upper']}, weakest V lower "
        f"{compact['weakest_V_lower']}, 2 open handoff ranges"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
