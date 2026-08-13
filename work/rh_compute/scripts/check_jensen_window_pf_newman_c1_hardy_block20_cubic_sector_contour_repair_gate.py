#!/usr/bin/env python3
"""Validate the cubic sector-contour repair gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.py"
)
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing sector contour artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate",
        "sector contour kind drift",
    )
    require(
        artifact["status"] == "exact_sign_aware_cubic_nonsaddle_contour_existence_and_rotation_lemma",
        "sector contour status drift",
    )
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["recursive_call_count"] == len(rows) == 374, "sector contour row count drift")
    require(scope["precision_ladder_decimal_digits"] == [90, 150], "sector contour precision drift")
    require(scope["precision_overlap_count"] == 374, "sector contour overlap count drift")
    require(aggregate["strict_convexity_count"] == 374, "sector contour convexity count drift")
    require(aggregate["negative_phi3_count"] == aggregate["positive_phi3_count"] == 187, "sector contour sign count drift")
    require(aggregate["b_slanted_count"] == aggregate["b_vertical_count"] == 187, "sector contour 67b count drift")
    require(aggregate["c_slanted_count"] == aggregate["c_vertical_count"] == 187, "sector contour 67c count drift")
    for name in (
        "minimum_second_derivative",
        "minimum_b_linear_decay",
        "minimum_c_linear_decay",
        "minimum_cubic_decay",
        "minimum_slanted_quadratic_decay_lower",
    ):
        require(Decimal(aggregate[name]) > 0, f"sector contour nonpositive aggregate: {name}")

    contracts = artifact["ray_contracts"]
    require(contracts["67b_phi3_negative"]["angle"] == "5*pi/6", "sector contour 67b slanted angle drift")
    require(contracts["67b_phi3_positive"]["angle"] == "pi/2", "sector contour 67b vertical angle drift")
    require(contracts["67c_phi3_negative"]["angle"] == "pi/2", "sector contour 67c vertical angle drift")
    require(contracts["67c_phi3_positive"]["angle"] == "pi/6", "sector contour 67c slanted angle drift")
    chain_roster = [row["chain"] for row in rows]
    require(chain_roster == sorted(chain_roster) and len(set(chain_roster)) == 374, "sector contour row order drift")
    for row in rows:
        require(row["precision_overlap"], "sector contour row precision drift")
        require(Decimal(row["second_derivative_minimum"]["decimal"]) > 0, "sector contour row convexity drift")
        require(Decimal(row["b_contract"]["linear_decay"]["decimal"]) > 0, "sector contour row 67b linear drift")
        require(Decimal(row["b_contract"]["cubic_decay"]["decimal"]) > 0, "sector contour row 67b cubic drift")
        require(Decimal(row["c_contract"]["linear_decay"]["decimal"]) > 0, "sector contour row 67c linear drift")
        require(Decimal(row["c_contract"]["cubic_decay"]["decimal"]) > 0, "sector contour row 67c cubic drift")
        if row["phi3_sign"] == "negative":
            require(row["b_contract"]["angle"] == "5*pi/6", "sector contour negative 67b angle drift")
            require(row["c_contract"]["angle"] == "pi/2", "sector contour negative 67c angle drift")
        else:
            require(row["b_contract"]["angle"] == "pi/2", "sector contour positive 67b angle drift")
            require(row["c_contract"]["angle"] == "pi/6", "sector contour positive 67c angle drift")

    require("no residues" in artifact["contour_lemma"]["analyticity"], "sector contour analyticity drift")
    require("abs(log(r))" in artifact["contour_lemma"]["sum_kernel"], "sector contour summation drift")
    require(artifact["paper"]["equations"] == [65, 67, 82, 83, 86, 89, 90, 91, 93, 94], "sector contour paper equation drift")
    for group in (artifact["paper"], artifact["source"], *artifact["sources"].values()):
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"sector contour dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"sector contour dependency hash drift: {group['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "sector contour builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "sector contour checker hash drift")
    for token in (
        "Status: exact sign-aware contour lemma",
        "5*pi/6",
        "pi/6",
        "quadratic endpoint-ray value",
        "1+|log(r)|",
        "not a proof",
    ):
        require(token in note, f"sector contour note token missing: {token}")
    require("does not bound the full-cubic versus quadratic endpoint remainder" in artifact["proof_boundary"], "sector contour boundary drift")
    print("validated cubic sector-contour repair gate: 374 strict-convexity calls, 187+187 sign-aware slanted rays")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
