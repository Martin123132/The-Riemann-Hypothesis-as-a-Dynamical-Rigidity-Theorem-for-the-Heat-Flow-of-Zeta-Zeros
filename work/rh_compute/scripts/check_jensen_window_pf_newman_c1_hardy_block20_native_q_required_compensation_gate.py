#!/usr/bin/env python3
"""Validate the block-20 native-q required-compensation gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing native-q artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate",
        "native-q kind drift",
    )
    require(
        artifact["status"] == "rigorous_374_call_native_q_required_compensation_inversion_enclosed",
        "native-q status drift",
    )
    scope = artifact["scope"]
    require(scope["recursive_chain_count"] == 374, "native-q chain count drift")
    require(scope["precision_ladder_decimal_digits"] == [180, 260], "native-q precision drift")
    require(scope["precision_overlap_count"] == 374, "native-q precision-overlap drift")
    require(scope["prior_full_corrected_defect_overlap_count"] == 374, "native-q prior overlap drift")
    require(scope["parent_upper_index"] == 104, "native-q parent index drift")
    require(scope["child_upper_indices"] == [1, 2], "native-q child index drift")

    identity = artifact["identity"]
    require("q_need=P_raw-M*C_adapted" in identity["required_q"], "native-q required identity drift")
    require("subtract-one cancels" in identity["defect_relation"], "native-q affine cancellation drift")
    aggregate = artifact["aggregate"]
    require(sum(aggregate["orientation_histogram"].values()) == 374, "native-q orientation histogram drift")
    require(sum(aggregate["component_sign_histogram"].values()) == 374, "native-q sign histogram drift")
    require(aggregate["required_compensations_excluding_zero"] == 374, "native-q zero exclusion drift")
    minimum = Decimal(aggregate["minimum_required_compensation_magnitude_lower"])
    maximum = Decimal(aggregate["maximum_required_compensation_magnitude_upper"])
    require(Decimal("0.004") < minimum < Decimal("0.006"), "native-q minimum scale drift")
    require(Decimal("0.20") < maximum < Decimal("0.21"), "native-q maximum scale drift")
    require(Decimal(aggregate["maximum_relative_to_full_q_upper"]) > 0, "native-q relative target drift")
    require(Decimal(aggregate["maximum_required_q_cross_gap_abs_upper"]) >= 0, "negative q cross-gap bound")
    require(Decimal(aggregate["maximum_orientation_identity_gap_abs_upper"]) >= 0, "negative identity-gap bound")

    rows = artifact["rows"]
    require(len(rows) == 374, "native-q row count drift")
    require([row["chain"] for row in rows] == sorted(row["chain"] for row in rows), "native-q order drift")
    for row in rows:
        require(row["precision_overlap"], "native-q precision nonoverlap row")
        require(row["prior_full_corrected_defect_overlap"], "native-q prior nonoverlap row")
        require(row["required_compensation_excludes_zero"], "native-q compensation contains zero")
        require(row["orientation"] in {"direct", "conjugated", "direct_subtract_one", "conjugated_subtract_one"}, "native-q orientation drift")
        require(row["required_native_q_compensation_real_sign"] in {"positive", "negative", "contains_zero"}, "native-q real sign drift")
        require(row["required_native_q_compensation_imag_sign"] in {"positive", "negative", "contains_zero"}, "native-q imag sign drift")

    for key in (
        "telemetry",
        "special_residual",
        "t5_real_erfc_residual",
        "full_corrected_defect",
        "selector_atlas",
        "point_ball_builder",
    ):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"native-q {key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "native-q builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "native-q checker hash drift")
    require(file_hash(REPO_ROOT / artifact["source"]["path"]) == artifact["source"]["sha256"], "native-q source hash drift")

    for token in (
        "Status: rigorous finite native-q inversion",
        "q_need = P_raw - M C",
        "Subtract-one cancels",
        "maximum oriented-defect identity gap",
        "complex target",
        "particular omitted source term equals",
        "proof of an outer Hardy estimate",
    ):
        require(token in note, f"native-q note token missing: {token}")
    require("does not identify that correction" in artifact["proof_boundary"], "native-q boundary drift")
    print(
        "validated Hardy block-20 native-q compensation gate: "
        f"374 nonzero targets, range {aggregate['minimum_required_compensation_magnitude_lower']} to "
        f"{aggregate['maximum_required_compensation_magnitude_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
