#!/usr/bin/env python3
"""Validate the order-ten compact H2-H24 unit-tile cache."""

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
import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as source  # noqa: E402
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
        "start_t": str(source.DEFAULT_START_T),
        "end_t": str(source.DEFAULT_END_T),
        "tile_width_t": str(source.DEFAULT_TILE_WIDTH_T),
        "initial_mode_bracket": [str(source.MODE_LEFT), str(source.MODE_RIGHT)],
        "mode_bisections": source.MODE_BISECTIONS,
        "max_moment": source.MAX_MOMENT,
        "precision_bits": source.PRECISION_BITS,
        "panels": source.PANELS,
        "window_y": source.WINDOW_Y,
        "eighth_envelope": str(compact.EIGHTH_ENVELOPE),
        "row_contract": source.ROW_CONTRACT,
    }


def validate_manifest(
    manifest: dict,
    cache_path: Path,
    records: list[dict],
) -> list[str]:
    issues = []
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_h2_h24_unit_cache"
    ):
        issues.append("manifest kind changed")
    if manifest.get("parameters") != expected_parameters():
        issues.append("manifest parameters changed")
    cache = manifest.get("cache", {})
    expected_path = cache_path.resolve().relative_to(source.REPO_ROOT).as_posix()
    if cache.get("path") != expected_path:
        issues.append("manifest cache path changed")
    if cache.get("row_count") != len(records):
        issues.append("manifest row count does not match JSONL")
    if cache.get("all_rows_passed") is not True:
        issues.append("manifest does not mark every row passed")
    if cache.get("h_derivative_orders") != [2, source.MAX_MOMENT]:
        issues.append("manifest derivative range changed")
    if cache.get("sha256") != source.sha256(cache_path):
        issues.append("manifest cache hash mismatch")
    if manifest.get("generator") != source.GENERATOR_PATH:
        issues.append("manifest generator changed")
    return issues


def validate_rows(
    records: list[dict],
    expected: list[tuple[int, Fraction, Fraction]],
) -> list[str]:
    issues = []
    if len(records) != len(expected):
        return [f"expected {len(expected)} rows, found {len(records)}"]
    derivative_keys = {str(order) for order in range(2, source.MAX_MOMENT + 1)}
    flint.ctx.prec = source.PRECISION_BITS
    for record, task in zip(records, expected):
        index, t_left, t_right = task
        prefix = f"row {index}"
        if (
            record.get("kind") != "order10_compact_h2_h24_unit_tile"
            or record.get("contract_id") != source.ROW_CONTRACT
            or record.get("index") != index
            or record.get("target_t_left") != str(t_left)
            or record.get("target_t_right") != str(t_right)
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
            if not bool(
                potential_jet_arb(arb_rational(mode_left), 1)[1]
                < arb_rational(t_left)
            ):
                issues.append(f"{prefix}: left mode does not under-cover tile")
            if not bool(
                potential_jet_arb(arb_rational(mode_right), 1)[1]
                > arb_rational(t_right)
            ):
                issues.append(f"{prefix}: right mode does not over-cover tile")
            actual_left = compact.interval_from_text(record["actual_t_left"])
            actual_right = compact.interval_from_text(record["actual_t_right"])
            if not bool(actual_left < arb_rational(t_left)):
                issues.append(f"{prefix}: serialized left image is not below tile")
            if not bool(actual_right > arb_rational(t_right)):
                issues.append(f"{prefix}: serialized right image is not above tile")
            for value in derivatives.values():
                compact.interval_from_text(value)
            for key in (
                "normalizer_lower",
                "minimum_tail_slope_lower",
                "maximum_tail_moment_upper",
                "maximum_simpson_error_upper",
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
    expected: list[tuple[int, Fraction, Fraction]],
) -> list[str]:
    issues = []
    sample_indices = sorted(
        {0, len(expected) // 4, len(expected) // 2, 3 * len(expected) // 4, len(expected) - 1}
    )
    source.initialize_worker()
    for index in sample_indices:
        rebuilt = source.tile_task(expected[index])
        if rebuilt != records[index]:
            issues.append(f"sample rebuild differs at row {index}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=source.DEFAULT_MANIFEST)
    parser.add_argument("--rebuild-samples", action="store_true")
    args = parser.parse_args()

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
    expected = source.deterministic_tasks(
        source.DEFAULT_START_T,
        source.DEFAULT_END_T,
        source.DEFAULT_TILE_WIDTH_T,
    )
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
        print(f"order-ten compact H2-H24 cache: {len(issues)} issues")
        for item in issues:
            print(f"- {item}")
        return 1
    suffix = " with five exact rebuilds" if args.rebuild_samples else ""
    print(
        "validated order-ten compact H2-H24 cache"
        f"{suffix}: {len(records)} contiguous unit tiles, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
