#!/usr/bin/env python3
"""Validate the sparse exact order-ten H0-H23 cache."""

from __future__ import annotations

import argparse
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
import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


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


def expected_parameters() -> dict:
    return {
        "start_t": str(source.START_T),
        "end_t": str(source.END_T),
        "step_t": str(source.STEP_T),
        "max_moment": source.MAX_MOMENT,
        "window_y": source.WINDOW_Y,
        "taylor_order": source.TAYLOR_ORDER,
        "profile_partition": {
            "medium": [str(source.START_T), str(source.MEDIUM_LAST_T)],
            "far": [str(source.FAR_FIRST_T), str(source.END_T)],
        },
        "profiles": {
            profile_name: {
                **source.profiles.PROFILE_SPECS[profile_name],
                "row_contract": source.row_contract(profile_name),
            }
            for profile_name in ("medium", "far")
        },
    }


def validate_manifest(
    manifest: dict,
    cache_path: Path,
    records: list[dict],
) -> list[str]:
    issues = []
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache"
        or manifest.get("parameters") != expected_parameters()
        or manifest.get("generator") != source.GENERATOR_PATH
    ):
        issues.append("manifest identity, parameters, or generator changed")
    cache = manifest.get("cache", {})
    profile_counts = {
        name: sum(record.get("profile") == name for record in records)
        for name in ("medium", "far")
    }
    expected_path = cache_path.resolve().relative_to(source.REPO_ROOT).as_posix()
    if (
        cache.get("path") != expected_path
        or cache.get("sha256") != source.sha256(cache_path)
        or cache.get("row_count") != len(records)
        or cache.get("profile_row_counts") != profile_counts
        or cache.get("all_rows_passed") is not True
        or cache.get("h_derivative_orders") != [0, source.MAX_MOMENT]
    ):
        issues.append("manifest cache contract changed")
    return issues


def validate_rows(
    records: list[dict],
    expected: list[tuple[int, Fraction, str]],
) -> list[str]:
    issues = []
    if len(records) != len(expected):
        return [f"expected {len(expected)} rows, found {len(records)}"]
    derivative_keys = {str(order) for order in range(source.MAX_MOMENT + 1)}
    for record, task in zip(records, expected):
        index, target, profile_name = task
        profile = source.profiles.PROFILE_SPECS[profile_name]
        flint.ctx.prec = profile["precision_bits"]
        prefix = f"row {index}"
        if (
            record.get("kind") != "order10_compact_sparse_point_h0_h23_jet"
            or record.get("contract_id") != source.row_contract(profile_name)
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("profile") != profile_name
            or record.get("passed") is not True
            or set(record.get("h_derivatives", {})) != derivative_keys
        ):
            issues.append(f"{prefix}: identity or derivative contract changed")
            continue
        try:
            mode_left = Fraction(record["mode_left"])
            mode_right = Fraction(record["mode_right"])
            target_arb = arb_rational(target)
            if not bool(
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


def validate_rebuild_samples(
    records: list[dict],
    expected: list[tuple[int, Fraction, str]],
) -> list[str]:
    issues = []
    indices = sorted({0, 1, 1519, 3038, 3039, len(expected) - 2, len(expected) - 1})
    source.initialize_worker()
    for index in indices:
        rebuilt = source.exact_task(expected[index])
        if rebuilt != records[index]:
            issues.append(f"exact sample rebuild differs at row {index}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=source.DEFAULT_MANIFEST)
    parser.add_argument("--rebuild-samples", action="store_true")
    args = parser.parse_args()

    issues = []
    expected = source.deterministic_tasks()
    if not args.cache.exists():
        issues.append(f"missing cache: {args.cache}")
        records = []
    else:
        try:
            records = load_records(args.cache)
        except RuntimeError as exc:
            issues.append(str(exc))
            records = []
    if records:
        issues.extend(validate_rows(records, expected))
    if not args.manifest.exists():
        issues.append(f"missing manifest: {args.manifest}")
    elif records:
        try:
            manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
            issues.extend(validate_manifest(manifest, args.cache, records))
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            issues.append(f"manifest parse failure: {exc}")
    if args.rebuild_samples and records and not issues:
        issues.extend(validate_rebuild_samples(records, expected))

    if issues:
        print(f"order-ten sparse exact H0-H23 cache: {len(issues)} issues")
        for item in issues:
            print(f"- {item}")
        return 1
    suffix = " with seven exact rebuilds" if args.rebuild_samples else ""
    print(
        "validated order-ten sparse exact H0-H23 cache"
        f"{suffix}: {len(records)} step-eight points, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
