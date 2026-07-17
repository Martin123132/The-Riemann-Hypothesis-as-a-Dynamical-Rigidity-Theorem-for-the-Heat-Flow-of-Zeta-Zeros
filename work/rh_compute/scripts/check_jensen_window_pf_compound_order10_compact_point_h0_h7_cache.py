#!/usr/bin/env python3
"""Validate an order-ten compact exact-point H0-H7 cache segment."""

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
import jensen_window_pf_compound_order10_compact_point_h0_h7_cache as source  # noqa: E402
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


def expected_parameters(
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
    profile_name: str,
) -> dict:
    profile = source.PROFILE_SPECS[profile_name]
    return {
        "start_t": str(start_t),
        "end_t": str(end_t),
        "step_t": str(step_t),
        "initial_mode_bracket": [str(source.MODE_LEFT), str(source.MODE_RIGHT)],
        "profile": profile_name,
        "precision_bits": profile["precision_bits"],
        "mode_bisections": profile["mode_bisections"],
        "panels": profile["panels"],
        "window_y": source.WINDOW_Y,
        "taylor_order": source.TAYLOR_ORDER,
        "max_moment": source.MAX_MOMENT,
        "row_contract": profile["row_contract"],
    }


def validate_manifest(
    manifest: dict,
    cache_path: Path,
    records: list[dict],
    *,
    start_t: Fraction,
    end_t: Fraction,
    step_t: Fraction,
    profile_name: str,
) -> list[str]:
    issues = []
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_point_h0_h7_cache"
    ):
        issues.append("manifest kind changed")
    if manifest.get("parameters") != expected_parameters(
        start_t, end_t, step_t, profile_name
    ):
        issues.append("manifest parameters changed")
    cache = manifest.get("cache", {})
    expected_path = cache_path.resolve().relative_to(source.REPO_ROOT).as_posix()
    if cache.get("path") != expected_path:
        issues.append("manifest cache path changed")
    if cache.get("row_count") != len(records):
        issues.append("manifest row count does not match JSONL")
    if cache.get("all_rows_passed") is not True:
        issues.append("manifest does not mark every row passed")
    if cache.get("h_derivative_orders") != [0, source.MAX_MOMENT]:
        issues.append("manifest derivative range changed")
    if cache.get("sha256") != source.sha256(cache_path):
        issues.append("manifest cache hash mismatch")
    if manifest.get("generator") != source.GENERATOR_PATH:
        issues.append("manifest generator changed")
    return issues


def validate_rows(
    records: list[dict],
    expected: list[tuple[int, Fraction, str]],
) -> list[str]:
    issues = []
    if len(records) != len(expected):
        return [f"expected {len(expected)} rows, found {len(records)}"]
    derivative_keys = {str(order) for order in range(source.MAX_MOMENT + 1)}
    profile_name = expected[0][2]
    profile = source.PROFILE_SPECS[profile_name]
    flint.ctx.prec = profile["precision_bits"]
    for record, task in zip(records, expected):
        index, target, task_profile = task
        prefix = f"row {index}"
        if (
            task_profile != profile_name
            or record.get("kind") != "order10_compact_point_h0_h7_jet"
            or record.get("profile") != profile_name
            or record.get("contract_id") != profile["row_contract"]
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("passed") is not True
        ):
            issues.append(f"{prefix}: identity or pass contract changed")
            continue
        derivatives = record.get("h_derivatives", {})
        if set(derivatives) != derivative_keys:
            issues.append(f"{prefix}: derivative orders changed")
            continue
        try:
            mode_left = Fraction(record["mode_left"])
            mode_right = Fraction(record["mode_right"])
            if not source.MODE_LEFT <= mode_left < mode_right <= source.MODE_RIGHT:
                issues.append(f"{prefix}: invalid mode interval")
                continue
            target_arb = arb_rational(target)
            if not bool(
                potential_jet_arb(arb_rational(mode_left), 1)[1]
                < target_arb
                < potential_jet_arb(arb_rational(mode_right), 1)[1]
            ):
                issues.append(f"{prefix}: mode bracket does not contain target")
            target_ball = flint.arb(record["target_t_ball"])
            if not bool(target_ball.contains(target_arb)):
                issues.append(f"{prefix}: serialized target ball misses target")
            for value in derivatives.values():
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
    sample_indices = sorted(
        {0, len(expected) // 4, len(expected) // 2, 3 * len(expected) // 4, len(expected) - 1}
    )
    source.initialize_worker(expected[0][2])
    for index in sample_indices:
        rebuilt = source.point_task(expected[index])
        if rebuilt != records[index]:
            issues.append(f"sample rebuild differs at row {index}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=source.DEFAULT_MANIFEST)
    parser.add_argument("--start-t", type=Fraction, default=source.DEFAULT_START_T)
    parser.add_argument("--end-t", type=Fraction, default=source.DEFAULT_END_T)
    parser.add_argument("--step-t", type=Fraction, default=source.DEFAULT_STEP_T)
    parser.add_argument(
        "--profile",
        choices=tuple(source.PROFILE_SPECS),
        default=source.DEFAULT_PROFILE,
    )
    parser.add_argument("--rebuild-samples", action="store_true")
    args = parser.parse_args()

    expected = source.deterministic_tasks(
        args.start_t,
        args.end_t,
        args.step_t,
        args.profile,
    )
    issues = []
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
            issues.extend(
                validate_manifest(
                    manifest,
                    args.cache,
                    records,
                    start_t=args.start_t,
                    end_t=args.end_t,
                    step_t=args.step_t,
                    profile_name=args.profile,
                )
            )
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            issues.append(f"manifest parse failure: {exc}")
    if args.rebuild_samples and records and not issues:
        issues.extend(validate_rebuild_samples(records, expected))

    if issues:
        print(f"order-ten compact point H0-H7 cache: {len(issues)} issues")
        for item in issues:
            print(f"- {item}")
        return 1
    suffix = " with five exact rebuilds" if args.rebuild_samples else ""
    print(
        "validated order-ten compact point H0-H7 cache"
        f"{suffix}: {len(records)} {args.profile}-profile points, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
