#!/usr/bin/env python3
"""Validate the exact-point order-nine H0-H8 half-grid cache."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order9_shifted_point_h0_h8_cache as point9  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


@dataclass(frozen=True)
class PointCacheIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> PointCacheIssue:
    return PointCacheIssue(section, code, str(detail))


def load_records(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def validate_manifest(
    path: Path,
    cache_path: Path,
    records: list[dict],
) -> tuple[dict, list[PointCacheIssue]]:
    if not path.exists():
        return {}, [issue("manifest", "missing-file", path)]
    manifest = json.loads(path.read_text(encoding="utf-8"))
    issues: list[PointCacheIssue] = []
    if manifest.get("kind") != (
        "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache"
    ):
        issues.append(issue("manifest", "bad-kind", manifest.get("kind")))
    expected_parameters = {
        "start_t": str(point9.DEFAULT_START_T),
        "end_t": str(point9.DEFAULT_END_T),
        "step_t": str(point9.DEFAULT_STEP_T),
        "precision_bits": point9.PRECISION_BITS,
        "mode_bisections": point9.MODE_BISECTIONS,
        "panels": point9.PANELS,
        "window_y": point9.WINDOW_Y,
        "taylor_order": point9.TAYLOR_ORDER,
        "max_moment": point9.MAX_MOMENT,
    }
    if manifest.get("parameters") != expected_parameters:
        issues.append(
            issue("manifest", "bad-parameters", manifest.get("parameters"))
        )
    cache = manifest.get("cache", {})
    expected_path = cache_path.resolve().relative_to(point9.REPO_ROOT).as_posix()
    if cache.get("path") != expected_path:
        issues.append(issue("manifest", "bad-cache-path", cache.get("path")))
    if cache.get("row_count") != len(records):
        issues.append(issue("manifest", "bad-row-count", cache.get("row_count")))
    if cache.get("all_rows_passed") is not True:
        issues.append(issue("manifest", "not-all-passed", None))
    if cache.get("h_derivative_orders") != [0, 8]:
        issues.append(
            issue(
                "manifest",
                "bad-derivative-orders",
                cache.get("h_derivative_orders"),
            )
        )
    if point9.sha256(cache_path) != cache.get("sha256"):
        issues.append(issue("manifest", "cache-hash-mismatch", cache_path))
    if records:
        maximum_panel = max(
            (record["maximum_panel_error_upper"] for record in records),
            key=float,
        )
        maximum_tail = max(
            (record["maximum_tail_moment_upper"] for record in records),
            key=float,
        )
        diagnostics = manifest.get("diagnostics", {})
        if diagnostics.get("maximum_panel_error_upper") != maximum_panel:
            issues.append(issue("manifest", "bad-panel-maximum", maximum_panel))
        if diagnostics.get("maximum_tail_moment_upper") != maximum_tail:
            issues.append(issue("manifest", "bad-tail-maximum", maximum_tail))
    return manifest, issues


def validate_rows(
    records: list[dict],
    expected: list[tuple],
) -> list[PointCacheIssue]:
    issues: list[PointCacheIssue] = []
    if len(records) != len(expected):
        return [
            issue(
                "rows",
                "bad-count",
                f"{len(records)} records for {len(expected)} tasks",
            )
        ]
    for record, task in zip(records, expected):
        index, target, mode_left, mode_right = task
        derivatives = record.get("h_derivatives", {})
        if (
            record.get("kind") != "order9_shifted_point_h0_h8_jet"
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("mode_left") != str(mode_left)
            or record.get("mode_right") != str(mode_right)
            or record.get("passed") is not True
            or set(derivatives) != {str(order) for order in range(9)}
        ):
            issues.append(issue("rows", "contract-mismatch", index))
            break
        mode_ball = (
            (arb_rational(mode_left) + arb_rational(mode_right)) / 2
            + flint.arb(
                0,
                (arb_rational(mode_right) - arb_rational(mode_left)) / 2,
            )
        )
        target_ball = potential_jet_arb(mode_ball, 1)[1]
        if not bool(target_ball.contains(arb_rational(target))):
            issues.append(issue("rows", "target-not-contained", target))
            break
        try:
            if Decimal(record["maximum_panel_error_upper"]) <= 0:
                raise ValueError("panel error is not positive")
            if Decimal(record["maximum_tail_moment_upper"]) <= 0:
                raise ValueError("tail moment is not positive")
            if Decimal(record["minimum_tail_slope_lower"]) <= 0:
                raise ValueError("tail slope is not positive")
        except Exception as exc:
            issues.append(issue("rows", "bad-diagnostic", f"{index}: {exc}"))
            break
    return issues


def validate_rebuild(
    records: list[dict],
    expected: list[tuple],
    *,
    workers: int,
) -> list[PointCacheIssue]:
    issues: list[PointCacheIssue] = []
    if workers <= 1:
        point9.initialize_worker()
        rebuilt = map(point9.point_task, expected)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=point9.initialize_worker,
        )
        rebuilt = executor.map(point9.point_task, expected, chunksize=1)
    try:
        for index, (actual, stored) in enumerate(zip(rebuilt, records)):
            if actual != stored:
                issues.append(issue("rebuild", "row-drift", index))
                break
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=point9.DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=point9.DEFAULT_MANIFEST)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--workers", type=int, default=point9.DEFAULT_WORKERS)
    args = parser.parse_args()

    flint.ctx.prec = point9.PRECISION_BITS
    issues: list[PointCacheIssue] = []
    if not args.cache.exists():
        issues.append(issue("cache", "missing-file", args.cache))
        records = []
    else:
        try:
            records = load_records(args.cache)
        except Exception as exc:
            issues.append(issue("cache", "parse-error", exc))
            records = []
    try:
        expected = point9.deterministic_tasks(
            point9.DEFAULT_START_T,
            point9.DEFAULT_END_T,
            point9.DEFAULT_STEP_T,
        )
    except Exception as exc:
        issues.append(issue("tasks", "rebuild-error", exc))
        expected = []
    if records:
        _, manifest_issues = validate_manifest(
            args.manifest,
            args.cache,
            records,
        )
        issues.extend(manifest_issues)
    issues.extend(validate_rows(records, expected))
    if args.rebuild and not issues:
        issues.extend(
            validate_rebuild(
                records,
                expected,
                workers=max(1, args.workers),
            )
        )
    if issues:
        for finding in issues:
            print(
                f"ORDER9-POINT-H {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine exact-point H0-H8 cache: {len(issues)} issues")
        return 1
    rebuild = " with full quadrature rebuild" if args.rebuild else ""
    print(
        f"validated order-nine exact-point H0-H8 cache{rebuild}: "
        f"{len(records)} rows, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
