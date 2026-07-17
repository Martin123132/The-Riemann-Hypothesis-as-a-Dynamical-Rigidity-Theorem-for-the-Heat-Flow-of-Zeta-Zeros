#!/usr/bin/env python3
"""Validate the propagated order-ten compact H0-H7 integer-grid cache."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache as source  # noqa: E402
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


def decimal_nonnegative(value: object) -> bool:
    try:
        number = Decimal(str(value))
    except (InvalidOperation, ValueError):
        return False
    return number.is_finite() and number >= 0


def validate_rows(
    records: list[dict],
    expected: list[tuple[int, Fraction]],
    contract_id: str,
) -> list[str]:
    issues = []
    if len(records) != len(expected):
        return [f"expected {len(expected)} rows, found {len(records)}"]
    keys = {str(order) for order in range(8)}
    flint.ctx.prec = source.PRECISION_BITS
    for record, task in zip(records, expected):
        index, target = task
        base = source.sparse_base(target)
        displacement = int(target - base)
        profile_name = source.exact_source.profile_for_target(base)
        prefix = f"row {index}"
        if (
            record.get("kind")
            != "order10_compact_propagated_point_h0_h7_jet"
            or record.get("source_contract_id") != contract_id
            or record.get("index") != index
            or record.get("target_t") != str(target)
            or record.get("sparse_base_t") != str(base)
            or record.get("sparse_displacement") != displacement
            or record.get("base_profile") != profile_name
            or record.get("passed") is not True
            or set(record.get("h_taylor_coefficients", {})) != keys
            or set(record.get("remainder_coefficient_upper", {})) != keys
        ):
            issues.append(f"{prefix}: identity or propagation contract changed")
            continue
        try:
            mode_left = Fraction(record["base_mode_left"])
            mode_right = Fraction(record["base_mode_right"])
            base_arb = arb_rational(base)
            if not bool(
                potential_jet_arb(arb_rational(mode_left), 1)[1]
                < base_arb
                < potential_jet_arb(arb_rational(mode_right), 1)[1]
            ):
                issues.append(f"{prefix}: exact-base mode bracket misses base")
            for value in record["h_taylor_coefficients"].values():
                compact.interval_from_text(value)
            for value in record["remainder_coefficient_upper"].values():
                if not decimal_nonnegative(value):
                    issues.append(f"{prefix}: invalid remainder upper")
                    break
            if displacement == 0 and any(
                Decimal(value) != 0
                for value in record["remainder_coefficient_upper"].values()
            ):
                issues.append(f"{prefix}: exact lattice row has nonzero remainder")
        except (KeyError, ValueError, TypeError, InvalidOperation) as exc:
            issues.append(f"{prefix}: parse failure: {exc}")
        if len(issues) >= 20:
            issues.append("validation stopped after 20 issues")
            break
    return issues


def validate_manifest(
    manifest: dict,
    records: list[dict],
    contract: dict,
    cache_path: Path,
) -> list[str]:
    issues = []
    expected_parameters = {
        "target_start_t": str(source.TARGET_START),
        "target_end_t": str(source.TARGET_END),
        "target_step_t": str(source.TARGET_STEP),
        "output_orders": [0, source.OUTPUT_MAX_MOMENT],
        "remainder_order": source.REMAINDER_MOMENT,
    }
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache"
        or manifest.get("parameters") != expected_parameters
        or manifest.get("source_contract") != contract
        or manifest.get("generator") != source.GENERATOR_PATH
    ):
        issues.append("manifest identity, parameters, or source contract changed")
    maxima = {
        str(order): max(
            (record["remainder_coefficient_upper"][str(order)] for record in records),
            key=float,
        )
        for order in range(8)
    }
    cache = manifest.get("cache", {})
    if (
        cache.get("path") != source.relative(cache_path)
        or cache.get("sha256") != source.sha256(cache_path)
        or cache.get("row_count") != len(records)
        or cache.get("all_rows_passed") is not True
        or cache.get("coefficient_orders") != [0, 7]
        or cache.get("maximum_remainder_coefficient_upper") != maxima
    ):
        issues.append("manifest cache contract or remainder maxima changed")
    return issues


def rebuild_indices(expected: list[tuple[int, Fraction]]) -> list[int]:
    targets = {
        5692,
        5699,
        5700,
        10000,
        29995,
        29996,
        30003,
        30004,
        38020,
        38026,
    }
    return sorted(int(Fraction(target) - source.TARGET_START) for target in targets)


def validate_rebuild(
    records: list[dict],
    expected: list[tuple[int, Fraction]],
    contract_id: str,
    *,
    rebuild_all: bool,
) -> list[str]:
    issues = []
    source.initialize_worker(contract_id)
    indices = range(len(expected)) if rebuild_all else rebuild_indices(expected)
    for index in indices:
        rebuilt = source.propagate_task(expected[index])
        if rebuilt != records[index]:
            issues.append(f"deterministic rebuild differs at row {index}")
            if len(issues) >= 20:
                break
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=source.DEFAULT_CACHE)
    parser.add_argument("--manifest", type=Path, default=source.DEFAULT_MANIFEST)
    parser.add_argument("--rebuild-samples", action="store_true")
    parser.add_argument("--rebuild-all", action="store_true")
    args = parser.parse_args()

    try:
        contract = source.validate_sources()
    except (RuntimeError, OSError, json.JSONDecodeError) as exc:
        print(f"order-ten propagated H0-H7 cache: invalid sources: {exc}")
        return 1
    expected = source.deterministic_tasks()
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
        issues.extend(validate_rows(records, expected, contract["id"]))
    if not args.manifest.exists():
        issues.append(f"missing manifest: {args.manifest}")
    elif records:
        try:
            manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
            issues.extend(
                validate_manifest(
                    manifest,
                    records,
                    contract,
                    args.cache,
                )
            )
        except (json.JSONDecodeError, OSError, KeyError, RuntimeError) as exc:
            issues.append(f"manifest parse failure: {exc}")
    if (args.rebuild_samples or args.rebuild_all) and records and not issues:
        issues.extend(
            validate_rebuild(
                records,
                expected,
                contract["id"],
                rebuild_all=args.rebuild_all,
            )
        )

    if issues:
        print(f"order-ten propagated H0-H7 cache: {len(issues)} issues")
        for item in issues:
            print(f"- {item}")
        return 1
    if args.rebuild_all:
        suffix = " with a full deterministic rebuild"
    elif args.rebuild_samples:
        suffix = " with ten deterministic sample rebuilds"
    else:
        suffix = ""
    print(
        "validated order-ten propagated H0-H7 cache"
        f"{suffix}: {len(records)} integer-grid points, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
