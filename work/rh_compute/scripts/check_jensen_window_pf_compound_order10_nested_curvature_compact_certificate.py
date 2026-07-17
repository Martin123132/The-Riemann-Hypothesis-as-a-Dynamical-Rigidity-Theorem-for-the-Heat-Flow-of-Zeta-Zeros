#!/usr/bin/env python3
"""Validate the compact order-ten first-summand curvature certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation, localcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_nested_curvature_compact_certificate as certificate  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
    arb_upper_text,
)


EXPECTED_THEOREM = "z_1''(t)<=4200/t^2 for every real 5700<=t<=38020"
EXPECTED_FORMULA = (
    "z''=2*w''-s''+phi(W)*W''-chi(W)*(W')^2 "
    "<=2*w''-s''+phi(W)*max(W'',0)"
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


def decimal_value(value: object, label: str) -> Decimal:
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise RuntimeError(f"invalid decimal {label}: {value}") from exc
    if not result.is_finite():
        raise RuntimeError(f"nonfinite decimal {label}: {value}")
    return result


def validate_rows(
    records: list[dict],
    expected: list[tuple],
    contract_id: str,
) -> tuple[list[str], dict | None]:
    issues = []
    if len(records) != len(expected):
        return [f"expected {len(expected)} blocks, found {len(records)}"], None
    largest = None
    weakest_w = None
    smallest_margin = None
    with localcontext() as context:
        context.prec = 100
        cap = Decimal(certificate.CURVATURE_CONSTANT)
        theorem_cap = Decimal(certificate.THEOREM_CURVATURE_CONSTANT)
        for record, task in zip(records, expected):
            index, anchor, right = task
            prefix = f"block {index}"
            if (
                record.get("kind") != "order10_compact_curvature_block"
                or record.get("source_contract_id") != contract_id
                or record.get("index") != index
                or record.get("anchor") != str(anchor)
                or record.get("expansion_anchor") != str(anchor)
                or record.get("right") != str(right)
                or record.get("width") != str(right - anchor)
                or record.get("passed") is not True
                or record.get("proof_formula") != EXPECTED_FORMULA
            ):
                issues.append(f"{prefix}: identity, geometry, or formula changed")
                continue
            profiles = record.get("point_profiles")
            if (
                not isinstance(profiles, list)
                or not profiles
                or not set(profiles) <= {"medium", "far"}
            ):
                issues.append(f"{prefix}: invalid point-profile provenance")
                continue
            try:
                scaled = decimal_value(
                    record["scaled_curvature_upper"],
                    f"{prefix} scaled upper",
                )
                margin = decimal_value(
                    record["curvature_margin_lower"],
                    f"{prefix} margin",
                )
                w_lower = decimal_value(record["W_lower"], f"{prefix} W lower")
                if scaled < 0 or scaled >= cap:
                    issues.append(f"{prefix}: scaled upper is outside [0,5500)")
                if scaled >= theorem_cap:
                    issues.append(f"{prefix}: scaled upper is not below 4200")
                if margin <= 0:
                    issues.append(f"{prefix}: curvature margin is not positive")
                if w_lower <= 0:
                    issues.append(f"{prefix}: W floor is not positive")
                if margin > cap - scaled:
                    issues.append(f"{prefix}: stated lower margin exceeds cap-upper")
                if largest is None or scaled > largest[0]:
                    largest = (scaled, record)
                if weakest_w is None or w_lower < weakest_w[0]:
                    weakest_w = (w_lower, record)
                if smallest_margin is None or margin < smallest_margin[0]:
                    smallest_margin = (margin, record)
            except (KeyError, RuntimeError) as exc:
                issues.append(str(exc))
            if len(issues) >= 20:
                issues.append("validation stopped after 20 issues")
                break
    if issues or largest is None or weakest_w is None or smallest_margin is None:
        return issues, None
    return issues, {
        "largest": largest,
        "weakest_w": weakest_w,
        "smallest_margin": smallest_margin,
    }


def validate_artifact(
    artifact: dict,
    records: list[dict],
    extrema: dict,
    expected_contract: dict,
    block_path: Path,
) -> list[str]:
    issues = []
    if (
        artifact.get("kind")
        != "jensen_window_pf_compound_order10_nested_curvature_compact_certificate"
        or artifact.get("status")
        != "rigorous order-ten first-summand curvature theorem on 5700<=t<=38020"
        or artifact.get("theorem") != EXPECTED_THEOREM
        or artifact.get("generator") != certificate.GENERATOR_PATH
        or artifact.get("checker") != certificate.CHECKER_PATH
    ):
        issues.append("artifact identity, status, or theorem changed")
    if artifact.get("source_contract") != expected_contract:
        issues.append("artifact source contract does not match live sources")
    cache = artifact.get("block_cache", {})
    if (
        cache.get("path") != certificate.relative(block_path)
        or cache.get("sha256") != certificate.sha256(block_path)
        or cache.get("row_count") != len(records)
        or cache.get("all_rows_passed") is not True
        or cache.get("source_contract_id") != expected_contract["id"]
    ):
        issues.append("artifact block-cache contract changed")
    rows = artifact.get("rows", [])
    if (
        len(rows) != 5
        or any(row.get("readiness") != "ready_to_apply" for row in rows)
        or [row.get("id") for row in rows]
        != [
            "co10nccc_01_sources",
            "co10nccc_02_partition",
            "co10nccc_03_formula",
            "co10nccc_04_blocks",
            "co10nccc_05_transition",
        ]
    ):
        issues.append("artifact proof ledger changed")
    summary = artifact.get("summary", {})
    largest = extrema["largest"][1]
    weakest_w = extrema["weakest_w"][1]
    smallest_margin = extrema["smallest_margin"][1]
    theorem_margin = str(
        Decimal(certificate.THEOREM_CURVATURE_CONSTANT)
        - Decimal(largest["scaled_curvature_upper"])
    )
    if (
        summary.get("blocks") != 18310
        or summary.get("unit_blocks") != 4300
        or summary.get("double_blocks") != 14010
        or summary.get("largest_scaled_curvature_upper")
        != largest["scaled_curvature_upper"]
        or summary.get("smallest_curvature_margin_lower")
        != smallest_margin["curvature_margin_lower"]
        or summary.get("theorem_margin_lower") != theorem_margin
        or summary.get("weakest_W_lower") != weakest_w["W_lower"]
        or summary.get("compact_first_summand_theorems") != 1
        or summary.get("full_half_line_theorems") != 0
        or summary.get("full_kernel_theorems") != 0
        or summary.get("rh_claims") != 0
    ):
        issues.append("artifact summary does not match block extrema")
    flint.ctx.prec = certificate.PRECISION_BITS
    transition = potential_jet_arb(
        arb_rational(certificate.SADDLE_RAY_START),
        1,
    )[1]
    if not bool(transition < arb_rational(certificate.COMPACT_END)):
        issues.append("live saddle transition no longer lies below compact endpoint")
    if summary.get("saddle_transition_upper") != arb_upper_text(transition):
        issues.append("artifact saddle-transition enclosure changed")
    return issues


def validate_rebuild_samples(
    records: list[dict],
    expected: list[tuple],
    contract_id: str,
) -> list[str]:
    issues = []
    indices = sorted(
        {
            0,
            1,
            4299,
            4300,
            12100,
            len(expected) - 2,
            len(expected) - 1,
        }
    )
    certificate.initialize_worker(
        str(certificate.H_CACHE),
        str(certificate.POINT_CACHE),
        contract_id,
    )
    for index in indices:
        rebuilt = certificate.block_task(expected[index])
        if rebuilt != records[index]:
            issues.append(f"exact sample rebuild differs at block {index}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--block-cache",
        type=Path,
        default=certificate.DEFAULT_BLOCK_CACHE,
    )
    parser.add_argument("--artifact", type=Path, default=certificate.DEFAULT_OUT)
    parser.add_argument("--rebuild-samples", action="store_true")
    args = parser.parse_args()

    issues = []
    try:
        expected_contract = certificate.source_contract()
    except (RuntimeError, OSError, json.JSONDecodeError) as exc:
        print(f"order-ten compact curvature certificate: invalid sources: {exc}")
        return 1
    expected = certificate.deterministic_blocks()
    if not args.block_cache.exists():
        issues.append(f"missing block cache: {args.block_cache}")
        records = []
    else:
        try:
            records = load_records(args.block_cache)
        except RuntimeError as exc:
            issues.append(str(exc))
            records = []
    extrema = None
    if records:
        row_issues, extrema = validate_rows(
            records,
            expected,
            expected_contract["id"],
        )
        issues.extend(row_issues)
    if not args.artifact.exists():
        issues.append(f"missing artifact: {args.artifact}")
    elif records and extrema is not None:
        try:
            artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
            issues.extend(
                validate_artifact(
                    artifact,
                    records,
                    extrema,
                    expected_contract,
                    args.block_cache,
                )
            )
        except (json.JSONDecodeError, OSError, KeyError, RuntimeError) as exc:
            issues.append(f"artifact parse failure: {exc}")
    if args.rebuild_samples and records and not issues:
        issues.extend(
            validate_rebuild_samples(
                records,
                expected,
                expected_contract["id"],
            )
        )

    if issues:
        print(f"order-ten compact curvature certificate: {len(issues)} issues")
        for item in issues:
            print(f"- {item}")
        return 1
    suffix = " with seven exact block rebuilds" if args.rebuild_samples else ""
    print(
        "validated order-ten compact curvature certificate"
        f"{suffix}: {len(records)} contiguous blocks, 0 issues, "
        "1 compact first-summand theorem, 0 full-kernel claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
