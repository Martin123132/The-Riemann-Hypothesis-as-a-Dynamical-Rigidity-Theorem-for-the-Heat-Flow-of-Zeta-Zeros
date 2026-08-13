#!/usr/bin/env python3
"""Validate the rigorous block-20 direct MIT=1 kernel gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing direct-kernel artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate", "direct kind drift")
    require(artifact["status"] == "rigorous_50_call_true_two_pi_direct_kernel_audit", "direct status drift")
    scope = artifact["scope"]
    require(scope["direct_chain_count"] == 50, "direct count drift")
    require(scope["upper_index"] == 104 and scope["terms_per_kernel"] == 105, "direct length drift")
    require(scope["precision_ladder_decimal_digits"] == [180, 260], "direct precision drift")
    require(scope["precision_overlap_count"] == 50, "direct overlap drift")
    require("last-decimal unit" in scope["logged_state_contract"], "direct print enclosure drift")
    require("p=4*atan(1)" in artifact["phase_definition"]["source"], "direct pi provenance drift")
    require(artifact["phase_definition"]["mathematical"] == "tpm=-2*pi", "direct true-pi target drift")

    aggregate = artifact["aggregate"]
    for key in (
        "maximum_source_roundoff_gap_abs_upper",
        "maximum_phase_normalization_gap_abs_upper",
        "maximum_total_true_gap_abs_upper",
        "maximum_phase_lipschitz_upper",
        "minimum_true_kernel_magnitude_lower",
    ):
        require(Decimal(aggregate[key]) > 0, f"nonpositive direct bound: {key}")
    require(Decimal(aggregate["maximum_source_roundoff_gap_abs_upper"]) < Decimal("1e-20"), "direct source roundoff too large")
    require(Decimal(aggregate["maximum_phase_normalization_gap_abs_upper"]) < Decimal("1e-20"), "direct phase shift too large")
    require(Decimal(aggregate["maximum_total_true_gap_abs_upper"]) < Decimal("1e-20"), "direct total gap too large")
    require(aggregate["all_phase_gaps_within_lipschitz_budget"] is True, "direct Lipschitz guard drift")

    rows = artifact["rows"]
    require(len(rows) == 50, "direct row count drift")
    require([row["chain"] for row in rows] == sorted(row["chain"] for row in rows), "direct row order drift")
    for row in rows:
        require(row["length_upper_index"] == 104 and row["term_count"] == 105, "direct row length drift")
        require(row["precision_overlap"] is True, "direct precision overlap lost")
        require(Decimal(row["phase_normalization_gap_abs_upper"]) <= Decimal(row["phase_lipschitz_upper"]), "direct row Lipschitz failure")

    for key in ("telemetry", "selector_atlas", "initial_adapter", "point_ball_builder"):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"direct {key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "direct builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "direct checker hash drift")
    require(file_hash(REPO_ROOT / artifact["source"]["path"]) == artifact["source"]["sha256"], "direct source drift")

    for token in (
        "Status: rigorous finite exact-point audit",
        "no `q` correction",
        "p=4*atan(1)",
        "one full last-decimal unit",
        "No arbitrary circle or polygon is",
        "phase-Lipschitz budget",
        "all 50 saved points",
        "remaining recursive `t5` saddle error",
        "prize-level conclusion",
    ):
        require(token in note, f"direct-kernel note token missing: {token}")
    require("50 direct block-20 calls only" in artifact["proof_boundary"], "direct boundary drift")
    print(
        "validated Hardy block-20 direct-kernel gate: "
        f"50 kernels, max true/logged gap {aggregate['maximum_total_true_gap_abs_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
