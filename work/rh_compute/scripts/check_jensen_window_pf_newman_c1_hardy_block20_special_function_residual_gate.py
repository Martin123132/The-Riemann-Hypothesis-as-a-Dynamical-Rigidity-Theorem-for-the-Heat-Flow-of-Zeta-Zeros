#!/usr/bin/env python3
"""Validate the rigorous block-20 special-function residual gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing special residual artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate",
        "special residual kind drift",
    )
    require(
        artifact["status"] == "rigorous_finite_point_special_residuals_and_correlated_q_replacement_enclosed",
        "special residual status drift",
    )
    scope = artifact["scope"]
    require(scope["recursive_chain_count"] == 374, "special residual chain count drift")
    require(scope["psi_evaluation_count"] == 2244, "special residual PSI count drift")
    require(scope["erf_evaluation_count"] == 748, "special residual ERF count drift")
    require(scope["saved_t1_t2_t4_bit_replay_count"] == 1122, "special residual replay drift")
    require(scope["precision_decimal_digits"] >= 160, "special residual precision too low")

    identity = artifact["psi_correlation_identity"]
    require(identity["exactly_cancelled_calls"] == ["P2", "P4"], "PSI cancellation drift")
    require(identity["surviving_calls"] == ["P1", "P3", "P5", "P6"], "PSI survivor drift")
    require("endpoint*(P6-P1)" in identity["identity"], "PSI identity drift")
    require("(1-i)*x/sqrt(2)" in artifact["pi_provenance"]["rigorous_erf_target"], "ERF pi provenance drift")

    aggregate = artifact["aggregate"]
    require(sum(aggregate["psi_source_branch_histogram"].values()) == 2244, "PSI branch histogram drift")
    require(sum(aggregate["erf_source_branch_histogram"].values()) == 748, "ERF branch histogram drift")
    for key in (
        "maximum_psi_residual_abs_upper",
        "maximum_psi_pair_P6_P1_abs_upper",
        "maximum_psi_pair_P3_P5_abs_upper",
        "maximum_erf_residual_abs_upper",
        "maximum_psi_q_delta_abs_upper",
        "maximum_erf_q_delta_abs_upper",
        "maximum_total_q_delta_abs_upper",
        "local_parent_state_shift_abs_upper",
    ):
        require(Decimal(aggregate[key]) > 0, f"nonpositive special residual bound: {key}")
    require(Decimal(aggregate["maximum_psi_pair_P6_P1_abs_upper"]) < Decimal("1e-25"), "P6-P1 residual symmetry lost")
    require(Decimal(aggregate["maximum_psi_pair_P3_P5_abs_upper"]) < Decimal("1e-25"), "P3-P5 residual symmetry lost")
    require(aggregate["q_delta_threshold_counts"]["greater_than_1e-4"] == 24, "material q-shift count drift")
    require(aggregate["q_delta_threshold_counts"]["greater_than_1e-3"] == 0, "unexpected milliscale q shift")
    require(
        aggregate["maximum_total_q_delta_abs_upper"] == aggregate["local_parent_state_shift_abs_upper"],
        "local affine propagation drift",
    )

    rows = artifact["rows"]
    require(len(rows) == 374, "special residual row count drift")
    require([row["chain"] for row in rows] == sorted(row["chain"] for row in rows), "special residual order drift")
    for row in rows:
        require(len(row["psi"]) == 6 and len(row["erf"]) == 2, "special residual call roster drift")
        require(Decimal(row["correlated_replacement"]["total_q_delta_abs_upper"]) > 0, "zero q replacement drift")
        require(Decimal(row["correlated_replacement"]["psi_residual_pair_P6_P1_abs_upper"]) >= 0, "P6-P1 pair bound drift")
        require(Decimal(row["correlated_replacement"]["psi_residual_pair_P3_P5_abs_upper"]) >= 0, "P3-P5 pair bound drift")
        require(len(row["correlated_replacement"]["source_q_hex"]) == 2, "source q payload drift")

    for key in ("probe_result", "probe_output", "telemetry", "selector_atlas", "local_cell_builder"):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"special residual {key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "special residual builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "special residual checker hash drift")
    require(file_hash(REPO_ROOT / artifact["source"]["path"]) == artifact["source"]["sha256"], "special residual source drift")

    for token in (
        "Status: rigorous finite-point PSI/ERF residuals",
        "1122 / 1122",
        "no unexplained pi is inserted",
        "P2=PSI(1-fracL)",
        "cancel exactly",
        "maximum total special q shift",
        "calls with total shift > 1e-4",
        "not yet the full q error",
        "prize-level conclusion",
    ):
        require(token in note, f"special residual note token missing: {token}")
    require("374 finite block-20 calls only" in artifact["proof_boundary"], "special residual boundary drift")
    print(
        "validated Hardy block-20 special-function residual gate: "
        f"2244 PSI, 748 ERF, max correlated q shift {aggregate['maximum_total_q_delta_abs_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
