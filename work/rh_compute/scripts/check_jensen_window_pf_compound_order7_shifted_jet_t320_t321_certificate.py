#!/usr/bin/env python3
"""Check the order-seven shifted-jet t=320..321 certificate."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.md"
)
EXPECTED_IDS = [
    "co7sj_01_exact_potential_quadrature",
    "co7sj_02_shifted_stable_jets",
    "co7sj_03_bootstrap_majorant",
    "co7sj_04_continuous_curvature",
    "co7sj_05_compact_handoff",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def issue(issues: list[str], message: str) -> None:
    issues.append(message)


def main() -> int:
    issues: list[str] = []
    if not ARTIFACT.exists():
        print(f"missing artifact: {ARTIFACT}")
        return 1
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate"
    ):
        issue(issues, "bad artifact kind")
    if artifact.get("date") != "2026-07-13":
        issue(issues, "bad artifact date")
    if artifact.get("status") != (
        "rigorous continuous fourth-nested first-summand curvature theorem on 320<=t<=321"
    ):
        issue(issues, "bad artifact status")

    parameters = artifact.get("parameters", {})
    expected_parameters = {
        "precision_bits": 384,
        "point_order": 10,
        "remainder_order": 11,
        "curvature_constant": 600,
        "start_t": "320",
        "end_t": "321",
        "block_count": 12,
        "block_width": "1/12",
    }
    for key, value in expected_parameters.items():
        if parameters.get(key) != value:
            issue(issues, f"bad parameter {key}")

    source = artifact.get("source_contract", {})
    cache_rel = source.get("h_cache", "")
    cache_path = REPO_ROOT / cache_rel
    if not cache_path.exists():
        issue(issues, "missing H source cache")
        cache_records = []
    else:
        if sha256(cache_path) != source.get("h_cache_sha256"):
            issue(issues, "H source cache hash mismatch")
        cache_records = [
            json.loads(line)
            for line in cache_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if len(cache_records) != source.get("h_cache_rows"):
            issue(issues, "H source cache row-count mismatch")
        if len(cache_records) != 920:
            issue(issues, "unexpected H source cache size")
        for index, record in enumerate(cache_records):
            if record.get("index") != index or record.get("passed") is not True:
                issue(issues, f"bad H cache row {index}")
                break
            if set(record.get("h_derivatives", {})) != {
                str(order) for order in range(2, 23)
            }:
                issue(issues, f"bad H derivative columns at row {index}")
                break
    if source.get("h_derivative_orders") != [2, 22]:
        issue(issues, "bad H derivative-order contract")
    if source.get("h_mode_range") != ["4619/5000", "933/1000"]:
        issue(issues, "bad H mode range")
    if source.get("tail_convexity_contract") != (
        "V''>0 for u>=1/100 and V'<0 for 0<u<=1/100"
    ):
        issue(issues, "bad tail convexity contract")
    if cache_records:
        flint.ctx.prec = 384
        first_t = compact.interval_from_text(cache_records[0]["t_left"])
        last_t = compact.interval_from_text(cache_records[-1]["t_right"])
        if not bool(first_t < 315 and last_t > 326):
            issue(issues, "H source cache does not cover the t+-5 endpoint collar")
        mode_left = Fraction(cache_records[0]["mode_left"])
        mode_right = Fraction(cache_records[-1]["mode_right"])
        if not bool(
            potential_jet_arb(arb_rational(mode_left), 1)[1] < 315
            and potential_jet_arb(arb_rational(mode_right), 1)[1] > 326
        ):
            issue(issues, "independent mode-to-t collar check failed")

    blocks = artifact.get("blocks", [])
    if len(blocks) != 12:
        issue(issues, "bad continuation block count")
    for index, block in enumerate(blocks):
        anchor = Fraction(320) + Fraction(index, 12)
        right = anchor + Fraction(1, 12)
        if block.get("index") != index:
            issue(issues, f"bad block index {index}")
        if block.get("anchor") != str(anchor) or block.get("right") != str(right):
            issue(issues, f"bad block endpoints {index}")
        if block.get("width") != "1/12" or block.get("passed") is not True:
            issue(issues, f"bad block status {index}")
        if Decimal(block.get("point_scaled_curvature", "600")) >= 600:
            issue(issues, f"point curvature did not clear at block {index}")
        scaled = Decimal(block.get("scaled_curvature_upper", "600"))
        margin = Decimal(block.get("curvature_margin_lower", "0"))
        if not scaled < 600 or not margin > 0:
            issue(issues, f"continued curvature did not clear at block {index}")
        if Decimal(block.get("curvature_upper", "0")) <= 0:
            issue(issues, f"bad curvature upper at block {index}")
        if Decimal(block.get("r_remainder_derivative_upper", "0")) <= 0:
            issue(issues, f"bad r remainder at block {index}")
        coordinates = block.get("point_coordinates", {})
        gaps = block.get("gaps", {})
        for name in ("B", "J", "R", "S", "T"):
            if Decimal(coordinates.get(name, "0")) <= 0:
                issue(issues, f"nonpositive point {name} at block {index}")
            gap = gaps.get(name, {})
            point_lower = Decimal(gap.get("point_lower", "0"))
            continued = Decimal(gap.get("continued_lower", "0"))
            floor = Decimal(gap.get("bootstrap_floor", "0"))
            remainder = Decimal(gap.get("remainder_derivative_upper", "0"))
            if not point_lower > continued > floor > 0 or not remainder > 0:
                issue(issues, f"bad {name} bootstrap at block {index}")
        quadrature = block.get("quadrature", {})
        if quadrature.get("shift_count") != 11:
            issue(issues, f"bad quadrature shift count at block {index}")
        if Decimal(quadrature.get("maximum_panel_error_upper", "1")) >= Decimal(
            "1e-15"
        ):
            issue(issues, f"panel error too large at block {index}")
        if Decimal(quadrature.get("maximum_tail_moment_upper", "1")) >= Decimal(
            "1e-15"
        ):
            issue(issues, f"tail error too large at block {index}")
        brackets = quadrature.get("mode_brackets", [])
        if len(brackets) != 11 or any(
            len(pair) != 2 or Fraction(pair[0]) >= Fraction(pair[1])
            for pair in brackets
        ):
            issue(issues, f"bad mode brackets at block {index}")
    for previous, current in zip(blocks, blocks[1:]):
        if previous.get("right") != current.get("anchor"):
            issue(issues, "continuation block gap")

    if blocks:
        largest = max(
            (row["scaled_curvature_upper"] for row in blocks), key=Decimal
        )
        weakest_margin = min(
            (row["curvature_margin_lower"] for row in blocks), key=Decimal
        )
        weakest_t = min(
            (row["gaps"]["T"]["continued_lower"] for row in blocks),
            key=Decimal,
        )
        summary = artifact.get("summary", {})
        if summary.get("largest_scaled_curvature_upper") != largest:
            issue(issues, "bad largest scaled summary")
        if summary.get("weakest_curvature_margin_lower") != weakest_margin:
            issue(issues, "bad weakest curvature summary")
        if summary.get("weakest_T_lower") != weakest_t:
            issue(issues, "bad weakest T summary")
        if summary.get("all_blocks_passed") is not True:
            issue(issues, "summary does not mark all blocks passed")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issue(issues, "bad certificate row ids")
    if len(rows) != 5:
        issue(issues, "bad certificate row count")
    else:
        if [row.get("readiness") for row in rows] != [
            "ready_to_apply",
            "ready_to_apply",
            "ready_to_apply",
            "ready_to_apply",
            "not_ready_to_apply",
        ]:
            issue(issues, "bad certificate readiness sequence")
        if rows[3].get("formula") != (
            "r_1''(t)<=600/t^2 for 320<=t<=321"
        ):
            issue(issues, "bad continuous theorem formula")
        if rows[4].get("formula") != "r_1''(t)<=600/t^2 for t>=321":
            issue(issues, "bad continuation handoff formula")

    if artifact.get("generator") != (
        "work/rh_compute/scripts/"
        "jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.py"
    ):
        issue(issues, "bad generator path")
    if artifact.get("checker") != (
        "work/rh_compute/scripts/"
        "check_jensen_window_pf_compound_order7_shifted_jet_t320_t321_certificate.py"
    ):
        issue(issues, "bad checker path")
    if not NOTE.exists():
        issue(issues, "missing markdown note")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for phrase in (
            "rigorous continuous fourth-nested first-summand curvature",
            "r_1''(t)<=600/t^2 for every real 320<=t<=321",
            "This is not a proof of the remaining",
        ):
            if phrase not in note:
                issue(issues, f"missing note phrase: {phrase}")

    if issues:
        for item in issues:
            print(f"ISSUE: {item}")
        print(
            "failed order-seven shifted-jet t=320..321 certificate: "
            f"{len(issues)} issues"
        )
        return 1
    print(
        "validated order-seven shifted-jet t=320..321 certificate: "
        "5 rows, 12 contiguous blocks, 11 shifts per anchor, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
