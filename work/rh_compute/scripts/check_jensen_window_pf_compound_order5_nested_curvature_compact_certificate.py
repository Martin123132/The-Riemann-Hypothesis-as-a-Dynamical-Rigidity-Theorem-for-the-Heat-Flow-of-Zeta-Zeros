#!/usr/bin/env python3
"""Validate the compact nested-curvature certificate for order five."""

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

from jensen_window_pf_compound_order5_nested_curvature_compact_certificate import (  # noqa: E402
    COLLAR_T,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    INITIAL_CENTRAL_TILE_COUNT,
    SOURCE_BRIDGE,
    SOURCE_COMPACT,
    TARGET_T,
    build_artifact,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)


@dataclass(frozen=True)
class CompactIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CompactIssue:
    return CompactIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[CompactIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[CompactIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_nested_curvature_compact_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "one open analytic ray" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_parameters = {
        "precision_bits": 256,
        "target_t": TARGET_T,
        "curvature_constant": CURVATURE_CONSTANT,
        "collar_t": COLLAR_T,
        "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append(issue("parameters", "drift", artifact.get("parameters")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 5,
        "ready_rows": 4,
        "open_rows": 1,
        "accepted_blocks": 36,
        "positive_stable_layers": 2,
        "compact_curvature_theorems": 1,
        "open_rays": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    contract = artifact.get("source_contract", {})
    if contract.get("cache_rows") != 107452:
        issues.append(issue("sources", "bad-cache-row-count", contract))
    if contract.get("cache_sha256") != (
        "d721a22738543dd2f62181a31732b26666d13eb3f8c4f1c8c46ead3e84ada4cf"
    ):
        issues.append(issue("sources", "bad-cache-hash", contract))
    for source in (SOURCE_COMPACT, SOURCE_BRIDGE):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    compact = artifact.get("compact", {})
    if compact.get("certified_t_range") != "320<=t<=V'(2)":
        issues.append(issue("compact", "bad-range", compact))
    if compact.get("target_t") != TARGET_T:
        issues.append(issue("compact", "bad-target", compact))
    finite_summary = compact.get("summary", {})
    if finite_summary.get("all_blocks_passed") is not True:
        issues.append(issue("compact", "blocks-not-passed", finite_summary))
    if finite_summary.get("accepted_blocks") != summary.get("accepted_blocks"):
        issues.append(issue("compact", "block-count-mismatch", finite_summary))
    for key in (
        "worst_margin_lower",
        "weakest_J_lower",
        "weakest_R_lower",
    ):
        try:
            value = Decimal(finite_summary.get(key, "nan"))
        except Exception as exc:
            issues.append(issue("compact", f"bad-{key}", exc))
        else:
            if not value.is_finite() or value <= 0:
                issues.append(issue("compact", f"nonpositive-{key}", value))
    try:
        largest = Decimal(
            finite_summary.get("largest_scaled_curvature_upper", "nan")
        )
    except Exception as exc:
        issues.append(issue("compact", "bad-largest-scaled", exc))
    else:
        if not largest.is_finite() or largest >= CURVATURE_CONSTANT:
            issues.append(issue("compact", "scaled-target-failed", largest))

    selected = compact.get("selected", [])
    if len(selected) < 5:
        issues.append(issue("compact", "too-few-selected-blocks", len(selected)))
    for row in selected:
        for key in ("J_lower", "R_lower", "margin_lower"):
            if Decimal(row.get(key, "nan")) <= 0:
                issues.append(issue("compact", f"selected-{key}", row))
        if Decimal(row.get("scaled_curvature_upper", "nan")) >= CURVATURE_CONSTANT:
            issues.append(issue("compact", "selected-scaled-target", row))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 5 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "u>=2" not in open_rows[0].get("formula", ""):
        issues.append(issue("rows", "bad-open-ray", open_rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[CompactIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous compact first-summand order-five curvature theorem",
        "This is not a proof of full order-five",
        "common three-unit",
        "B^(r)=integral_[-1,1]",
        "J^(r)=2*B^(r)-integral_[-1,1]",
        "R^(r)=3*B^(r)-integral_[-1,1]",
        "107452",
        "Adaptive assembly accepts `36` central blocks.",
        "320<=t<=V'(2)",
        "q_1''(t)<=60/t^2 for every 320<=t<=V'(2)",
        "q_1''(t)<=60/t^2 for every mode u>=2",
        "only unproved part",
        "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
    )
    issues: list[CompactIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "the analytic ray is proved",
        "full order-five entry is complete",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER5-NESTED-COMPACT {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five nested curvature compact certificate: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    finite = artifact["compact"]["summary"]
    print(
        "validated order-five nested curvature compact certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['accepted_blocks']} interval blocks, "
        f"{summary['positive_stable_layers']} positive stable layers, "
        f"{summary['compact_curvature_theorems']} compact curvature theorem, "
        f"{summary['open_rays']} open ray, "
        f"largest scaled upper {finite['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
