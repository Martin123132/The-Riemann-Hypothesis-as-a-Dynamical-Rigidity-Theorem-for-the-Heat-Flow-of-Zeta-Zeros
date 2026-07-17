#!/usr/bin/env python3
"""Validate the order-seven nested-curvature compact certificate."""

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
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
from jensen_window_pf_compound_order7_nested_curvature_compact_certificate import (  # noqa: E402
    DEFAULT_EXTENSION_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_RIGHT_COLLAR,
    SOURCE_BRIDGE,
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
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order7_nested_curvature_compact_certificate"
    ):
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))

    summary = artifact.get("summary", {})
    for key, value in {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "compact_curvature_theorems": 1,
        "open_finite_rays": 1,
    }.items():
        if summary.get(key) != value:
            issues.append(Issue("summary", key, str(summary.get(key))))

    tasks = deterministic_tasks()
    extension = load_extension_cache(cache_path, tasks)
    if len(extension) != len(tasks):
        issues.append(Issue("cache", "h11-h12-row-count", f"{len(extension)}/{len(tasks)}"))
        return artifact, issues

    contract = artifact.get("source_contract", {})
    expected_contract = {
        "base_cache_rows": len(tasks),
        "h9_h10_cache_rows": len(tasks),
        "h11_h12_cache_rows": len(tasks),
        "base_cache_sha256": sha256(order4_compact.DEFAULT_CACHE),
        "h9_h10_cache_sha256": sha256(order6_compact.DEFAULT_EXTENSION_CACHE),
        "h11_h12_cache_sha256": sha256(cache_path),
        "right_collar_sha256": sha256(right_collar_path),
        "compact_bridge_sha256": sha256(SOURCE_BRIDGE),
    }
    for key, value in expected_contract.items():
        if contract.get(key) != value:
            issues.append(Issue("source-contract", key, str(contract.get(key))))

    bridge = json.loads(SOURCE_BRIDGE.read_text(encoding="utf-8"))
    if bridge.get("summary", {}).get("all_blocks_passed") is not True:
        issues.append(Issue("source-contract", "compact-bridge", str(SOURCE_BRIDGE)))

    compact = artifact.get("compact", {})
    blocks = compact.get("blocks", [])
    if compact.get("theorem") != (
        "r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2)"
    ):
        issues.append(Issue("compact", "theorem", str(compact.get("theorem"))))
    if compact.get("mode_range", [None, None])[-1] != "2":
        issues.append(Issue("compact", "mode-end", str(compact.get("mode_range"))))
    if compact.get("all_blocks_passed") is not True:
        issues.append(Issue("compact", "blocks", "not all passed"))
    if compact.get("block_count") != len(blocks):
        issues.append(
            Issue(
                "compact",
                "block-count",
                f"{compact.get('block_count')}/{len(blocks)}",
            )
        )
    if summary.get("compact_blocks") != len(blocks):
        issues.append(
            Issue(
                "summary",
                "compact-blocks",
                f"{summary.get('compact_blocks')}/{len(blocks)}",
            )
        )
    if Decimal(str(compact.get("largest_scaled_curvature_upper", "600"))) >= 600:
        issues.append(
            Issue(
                "compact",
                "scaled-upper",
                str(compact.get("largest_scaled_curvature_upper")),
            )
        )
    if Decimal(str(compact.get("weakest_T_lower", "0"))) <= 0:
        issues.append(Issue("compact", "T-lower", str(compact.get("weakest_T_lower"))))

    positive_fields = (
        "B_lower",
        "J_lower",
        "R_lower",
        "S_lower",
        "T_lower",
        "d_lower",
        "target_lower",
        "margin_lower",
    )
    for index, block in enumerate(blocks):
        if block.get("passed") is not True:
            issues.append(Issue("blocks", "failed", str(index)))
        for field in positive_fields:
            if Decimal(str(block.get(field, "0"))) <= 0:
                issues.append(Issue("blocks", field, str(index)))
        if Decimal(str(block.get("scaled_curvature_upper", "600"))) >= 600:
            issues.append(Issue("blocks", "scaled-upper", str(index)))
        diagnostics = block.get("cover_diagnostics", {})
        for side in ("left_t_collar_ball", "right_t_collar_ball"):
            try:
                collar_width = order4_compact.flint.arb(
                    str(diagnostics.get(side, ""))
                )
            except Exception:
                issues.append(Issue("blocks", side, str(index)))
            else:
                if not bool(collar_width > 5):
                    issues.append(Issue("blocks", side, str(index)))
    for index, (previous, current) in enumerate(zip(blocks, blocks[1:])):
        if previous.get("central_mode", [None, None])[1] != current.get(
            "central_mode", [None, None]
        )[0]:
            issues.append(Issue("blocks", "mode-gap", str(index)))
        if previous.get("central_tile_index_last", -2) + 1 != current.get(
            "central_tile_index_first", -1
        ):
            issues.append(Issue("blocks", "tile-gap", str(index)))

    try:
        base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
        high = order6_compact.load_extension_cache(
            order6_compact.DEFAULT_EXTENSION_CACHE,
            tasks,
        )
        if len(base) != len(tasks) or len(high) != len(tasks):
            raise RuntimeError("existing compact source caches are incomplete")
        collar = build_right_collar(right_collar_path)
        base, high, extension = append_right_collar(
            base,
            high,
            extension,
            collar,
        )
        canonical = build_artifact(
            base,
            high,
            extension,
            cache_path,
            right_collar_path,
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
            "Status: rigorous first-summand order-seven",
            "H^(2),...,H^(12)",
            "t+-5",
            "r_1''(t)<=600/t^2",
            "finite or asymptotic",
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
        args.artifact,
        args.note,
        args.extension_cache,
        args.right_collar,
    )
    for item in issues:
        print(f"ORDER7-NESTED-COMPACT {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-seven nested compact certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-seven nested compact certificate: "
        f"{summary['compact_blocks']} blocks, 0 issues, "
        f"{summary['compact_curvature_theorems']} compact theorem, "
        f"{summary['open_finite_rays']} open finite ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
