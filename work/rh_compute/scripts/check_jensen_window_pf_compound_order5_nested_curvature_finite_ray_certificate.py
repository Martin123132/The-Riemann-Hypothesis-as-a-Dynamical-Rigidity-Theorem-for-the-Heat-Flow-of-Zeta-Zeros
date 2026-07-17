#!/usr/bin/env python3
"""Validate the finite-ray nested-curvature certificate for order five."""

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

from jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_EXTENSION_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_RAY_CACHE,
    SOURCE_COMPACT,
    SOURCE_EXACT_CORRIDOR,
    build_artifact,
    extension_tasks,
    load_cache,
    ray_tasks,
    sha256,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)


@dataclass(frozen=True)
class RayIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> RayIssue:
    return RayIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_cache(
    path: Path,
    tasks: list,
    expected_hash: str,
    section: str,
) -> list[RayIssue]:
    issues: list[RayIssue] = []
    if not path.exists():
        return [issue(section, "missing-cache", path)]
    actual_hash = sha256(path)
    if actual_hash != expected_hash:
        issues.append(issue(section, "hash-drift", actual_hash))
    try:
        records = load_cache(path, tasks)
    except Exception as exc:
        issues.append(issue(section, "cache-contract", exc))
        return issues
    if len(records) != len(tasks):
        issues.append(issue(section, "incomplete-cache", len(records)))
    for row in records:
        if row.get("passed") is not True:
            issues.append(issue(section, "failed-row", row.get("index")))
            break
    return issues


def validate_artifact(path: Path) -> tuple[dict, list[RayIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[RayIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "one open asymptotic ray" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 5,
        "ready_rows": 4,
        "open_rows": 1,
        "extension_tiles": 100,
        "initial_collar_blocks": 1,
        "exact_corridor_blocks": 1850,
        "finite_ray_theorems": 1,
        "open_asymptotic_rays": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        largest = Decimal(summary.get("largest_scaled_curvature_upper", "nan"))
    except Exception as exc:
        issues.append(issue("summary", "bad-largest-scaled", exc))
    else:
        if not largest.is_finite() or largest >= CURVATURE_CONSTANT:
            issues.append(issue("summary", "scaled-target-failed", largest))

    extension = artifact.get("extension", {})
    ray = artifact.get("finite_ray", {})
    if extension.get("rows") != len(extension_tasks()):
        issues.append(issue("extension", "bad-row-count", extension))
    if ray.get("block_count") != len(ray_tasks()):
        issues.append(issue("ray", "bad-block-count", ray))
    if ray.get("mode_range") != ["2001/1000", "20"]:
        issues.append(issue("ray", "bad-mode-range", ray))
    if ray.get("all_blocks_passed") is not True:
        issues.append(issue("ray", "blocks-not-passed", ray))
    for key in ("weakest_margin_lower", "weakest_J_lower", "weakest_R_lower"):
        try:
            value = Decimal(ray.get(key, "nan"))
        except Exception as exc:
            issues.append(issue("ray", f"bad-{key}", exc))
        else:
            if not value.is_finite() or value <= 0:
                issues.append(issue("ray", f"nonpositive-{key}", value))

    collar = artifact.get("initial_collar", {})
    if collar.get("central_mode") != ["2", "2001/1000"]:
        issues.append(issue("collar", "bad-mode-range", collar))
    if collar.get("passed") is not True:
        issues.append(issue("collar", "not-passed", collar))
    for key in ("margin_lower", "J_lower", "R_lower"):
        if Decimal(collar.get(key, "nan")) <= 0:
            issues.append(issue("collar", f"nonpositive-{key}", collar))

    envelope = extension.get("envelope", {})
    if envelope.get("strict") is not True or envelope.get("cap") != "1/50000":
        issues.append(issue("extension", "bad-envelope", envelope))

    issues.extend(
        validate_cache(
            DEFAULT_EXTENSION_CACHE,
            extension_tasks(),
            "3bd7b260926c48682930ad03f2233a2a79c05019fb8e9a04f4c581f3a916be19",
            "extension-cache",
        )
    )
    issues.extend(
        validate_cache(
            DEFAULT_RAY_CACHE,
            ray_tasks(),
            "ce365b6a6bbc527f5952579e865d8765d5368267a363c238f693d0a52061a660",
            "ray-cache",
        )
    )
    for source in (SOURCE_COMPACT, SOURCE_EXACT_CORRIDOR):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 5 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "u>=20" not in open_rows[0].get("formula", ""):
        issues.append(issue("rows", "bad-open-ray", open_rows))

    try:
        canonical = build_artifact(workers=1, overwrite=False)
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        canonical.get("parameters", {})["workers"] = artifact.get(
            "parameters", {}
        ).get("workers")
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[RayIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous first-summand nested-curvature theorem on",
        "This is not a proof of",
        "One hundred additional",
        "abs(V^(8)/V''^4)<1/50000 on 2<=u<=2.002",
        "q_1''(t)<=60/t^2 for every 2<=u<=2.001",
        "`1850` blocks from `u=2.001` to `u=20`",
        "q_1''(t)<=60/t^2 for every mode 2<=u<=20",
        "sole remaining part",
        "q_1''(t)<=60/t^2 for every mode u>=20",
        "0<x_r<=1",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
        "outputs/formal_core.md",
    )
    issues: list[RayIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "the asymptotic ray is proved",
        "complete order-five entry is complete",
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
                f"ORDER5-NESTED-FINITE-RAY {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five nested curvature finite-ray certificate: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    print(
        "validated order-five nested curvature finite-ray certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['extension_tiles']} extension tiles, "
        f"{summary['exact_corridor_blocks']} exact-corridor blocks, "
        f"{summary['finite_ray_theorems']} finite-ray theorem, "
        f"{summary['open_asymptotic_rays']} open asymptotic ray, "
        f"largest scaled upper {summary['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
