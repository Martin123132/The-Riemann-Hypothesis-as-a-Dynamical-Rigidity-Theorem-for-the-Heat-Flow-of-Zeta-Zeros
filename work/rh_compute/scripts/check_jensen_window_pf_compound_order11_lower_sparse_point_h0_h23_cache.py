#!/usr/bin/env python3
"""Independently validate the order-eleven step-two exact H0-H23 cache."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache as source  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


LOWER_BRIDGE_START = Fraction(1252)
LOWER_BRIDGE_END = Fraction(5700)
SHIFT_RADIUS = 9
HALF_GRID_STEP = Fraction(1, 2)


def load_records(path: Path) -> list[dict]:
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise RuntimeError(f"invalid JSONL row {line_number}") from exc
    return records


def validate_rows(
    records: list[dict],
    expected: list[tuple[int, Fraction]],
) -> list[str]:
    issues = []
    if len(records) > len(expected):
        return [f"expected at most {len(expected)} rows, found {len(records)}"]
    derivative_keys = {str(order) for order in range(source.MAX_MOMENT + 1)}
    flint.ctx.prec = source.PRECISION_BITS
    for record, task in zip(records, expected):
        index, target = task
        prefix = f"row {index}"
        if (
            record.get("kind") != "order11_lower_sparse_point_h0_h23_jet"
            or record.get("contract_id") != source.ROW_CONTRACT
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != derivative_keys
        ):
            issues.append(f"{prefix}: identity or derivative contract changed")
            continue
        try:
            mode_left = Fraction(record["mode_left"])
            mode_right = Fraction(record["mode_right"])
            target_arb = arb_rational(target)
            if not mode_left < mode_right:
                issues.append(f"{prefix}: empty mode bracket")
            elif not bool(
                potential_jet_arb(arb_rational(mode_left), 1)[1]
                < target_arb
                < potential_jet_arb(arb_rational(mode_right), 1)[1]
            ):
                issues.append(f"{prefix}: mode bracket misses target")
            if not bool(flint.arb(record["target_t_ball"]).contains(target_arb)):
                issues.append(f"{prefix}: serialized target ball misses target")
            for value in record["h_derivatives"].values():
                compact.interval_from_text(value)
            for key in (
                "maximum_panel_error_upper",
                "maximum_tail_moment_upper",
                "minimum_tail_slope_lower",
            ):
                value = float(record[key])
                if not math.isfinite(value) or value < 0:
                    issues.append(f"{prefix}: invalid {key}")
        except (KeyError, ValueError, TypeError) as exc:
            issues.append(f"{prefix}: parse failure: {exc}")
        if len(issues) >= 20:
            issues.append("validation stopped after 20 issues")
            break
    return issues


def validate_manifest(
    manifest: dict,
    cache_path: Path,
    records: list[dict],
    expected: list[tuple[int, Fraction]],
) -> list[str]:
    expected_parameters = {
        "start_t": str(source.DEFAULT_START_T),
        "end_t": str(source.DEFAULT_END_T),
        "step_t": str(source.DEFAULT_STEP_T),
        "precision_bits": source.PRECISION_BITS,
        "mode_bisections": source.MODE_BISECTIONS,
        "panels": source.PANELS,
        "window_y": source.WINDOW_Y,
        "taylor_order": source.TAYLOR_ORDER,
        "max_moment": source.MAX_MOMENT,
        "row_contract": source.ROW_CONTRACT,
    }
    cache = manifest.get("cache", {})
    expected_path = cache_path.resolve().relative_to(source.REPO_ROOT).as_posix()
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache"
        or manifest.get("parameters") != expected_parameters
        or manifest.get("generator") != source.GENERATOR_PATH
    ):
        return ["manifest identity, parameters, or generator changed"]
    if (
        len(records) != len(expected)
        or cache.get("path") != expected_path
        or cache.get("sha256") != source.sha256(cache_path)
        or cache.get("row_count") != len(expected)
        or cache.get("all_rows_passed") is not True
        or cache.get("h_derivative_orders") != [0, source.MAX_MOMENT]
    ):
        return ["manifest cache contract changed"]
    return []


def validate_propagation_geometry() -> list[str]:
    anchors = sorted(
        target
        for _, target in source.deterministic_tasks(
            source.DEFAULT_START_T,
            source.DEFAULT_END_T,
            source.DEFAULT_STEP_T,
        )
    )
    target_start = LOWER_BRIDGE_START - SHIFT_RADIUS
    target_end = LOWER_BRIDGE_END + SHIFT_RADIUS
    target = target_start
    while target <= target_end:
        index = bisect_left(anchors, target)
        candidates = anchors[max(0, index - 1) : min(len(anchors), index + 1)]
        distance = min(abs(target - anchor) for anchor in candidates)
        if distance > 1:
            return [
                f"half-grid target {target} is {distance} from its nearest anchor"
            ]
        target += HALF_GRID_STEP
    if min(anchors) != target_start + 1 or max(anchors) != target_end - 1:
        return ["sparse endpoint collar no longer has exact distance-one coverage"]
    return []


def validate_rebuild_samples(
    records: list[dict],
    expected: list[tuple[int, Fraction]],
) -> list[str]:
    if not records:
        return ["cannot rebuild samples from an empty cache"]
    indices = sorted({0, len(records) // 2, len(records) - 1})
    source.initialize_worker()
    issues = []
    for index in indices:
        rebuilt = source.exact_task(expected[index])
        if rebuilt != records[index]:
            issues.append(f"exact sample rebuild differs at row {index}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=source.DEFAULT_MANIFEST)
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument("--rebuild-samples", action="store_true")
    args = parser.parse_args()

    expected = source.deterministic_tasks(
        source.DEFAULT_START_T,
        source.DEFAULT_END_T,
        source.DEFAULT_STEP_T,
    )
    issues = validate_propagation_geometry()
    if not args.cache.exists():
        records = []
        issues.append(f"missing cache: {args.cache}")
    else:
        try:
            records = load_records(args.cache)
        except RuntimeError as exc:
            records = []
            issues.append(str(exc))
    if records:
        issues.extend(validate_rows(records, expected))
    if args.require_complete and len(records) != len(expected):
        issues.append(f"expected {len(expected)} rows, found {len(records)}")
    if args.manifest.exists():
        try:
            manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
            issues.extend(validate_manifest(manifest, args.cache, records, expected))
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            issues.append(f"manifest parse failure: {exc}")
    elif args.require_complete:
        issues.append(f"missing manifest: {args.manifest}")
    if args.rebuild_samples and not issues:
        issues.extend(validate_rebuild_samples(records, expected))

    if issues:
        print(f"order-eleven lower sparse exact H0-H23 cache: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    status = "complete" if len(records) == len(expected) else "valid resumable prefix"
    suffix = " with three exact rebuilds" if args.rebuild_samples else ""
    print(
        "validated order-eleven lower sparse exact H0-H23 cache "
        f"({status}){suffix}: {len(records)}/{len(expected)} rows, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
