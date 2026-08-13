#!/usr/bin/env python3
"""Validate the block-20 W2--W5 source/paper correspondence gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing W2--W5 artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate",
        "W2--W5 kind drift",
    )
    require(
        artifact["status"] == "published_w2_w5_orientation_and_finite_component_residual_rigorously_enclosed",
        "W2--W5 status drift",
    )
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["recursive_call_count"] == len(rows) == 374, "W2--W5 row count drift")
    require(scope["precision_ladder_decimal_digits"] == [90, 150], "W2--W5 precision ladder drift")
    require(scope["precision_overlap_count"] == 374, "W2--W5 precision overlap drift")
    require(aggregate["w3_call_count"] + aggregate["w4_call_count"] == 374, "W2--W5 selector count drift")
    require(aggregate["w3_call_count"] > 0 and aggregate["w4_call_count"] > 0, "W2--W5 selector branch lost")
    require(
        aggregate["improved_call_count"] + aggregate["worsened_call_count"] + aggregate["indeterminate_call_count"] == 374,
        "W2--W5 classification count drift",
    )
    require(0 <= aggregate["paper_residual_excluding_zero_count"] <= 374, "W2--W5 exclusion count drift")
    require(Decimal(aggregate["minimum_paper_residual_magnitude_lower"]) >= 0, "W2--W5 minimum residual drift")
    require(Decimal(aggregate["maximum_paper_residual_magnitude_upper"]) > 0, "W2--W5 maximum residual drift")
    require(Decimal(aggregate["maximum_source_component_assembly_gap_upper"]) < Decimal("1e-30"), "W2--W5 source assembly drift")
    require(Decimal(aggregate["maximum_component_transport_identity_gap_upper"]) < Decimal("1e-80"), "W2--W5 component identity drift")
    require(Decimal(aggregate["maximum_residual_two_construction_gap_upper"]) < Decimal("1e-80"), "W2--W5 residual identity drift")
    require(Decimal(aggregate["maximum_endpoint_source_to_pi_transport_upper"]) < Decimal("1e-28"), "W2--W5 endpoint pi transport drift")
    require(Decimal(aggregate["maximum_xi_source_rounding_transport_upper"]) < Decimal("1e-28"), "W2--W5 xi transport drift")
    require(Decimal(aggregate["maximum_parent_source_to_pi_transport_upper"]) < Decimal("1e-26"), "W2--W5 parent pi transport drift")
    require(Decimal(aggregate["maximum_eq69_source_to_pi_transport_upper"]) < Decimal("1e-26"), "W2--W5 integral pi transport drift")

    chain_roster = [row["chain"] for row in rows]
    require(chain_roster == sorted(chain_roster) and len(set(chain_roster)) == 374, "W2--W5 row order drift")
    for row in rows:
        require(row["selector"] in ("W3", "W4"), "W2--W5 row selector drift")
        require(row["precision_overlap"], "W2--W5 row precision drift")
        require(row["magnitude_classification"] in ("improved", "worsened", "indeterminate"), "W2--W5 row class drift")
        require(Decimal(row["source_component_assembly_gap_upper"]) < Decimal("1e-30"), "W2--W5 row assembly drift")
        require(Decimal(row["component_transport_identity_gap_upper"]) < Decimal("1e-80"), "W2--W5 row component identity drift")
        require(Decimal(row["residual_two_construction_gap_upper"]) < Decimal("1e-80"), "W2--W5 row residual identity drift")
        require(Decimal(row["xi_source_rounding_transport_upper"]) < Decimal("1e-28"), "W2--W5 row xi transport drift")

    correspondence = artifact["correspondence"]
    require(correspondence["W2"]["paper"] == "equation (89)", "W2 correspondence drift")
    require("conj(t4)" in correspondence["W3"]["source"] and "conj(t4)" in correspondence["W4"]["source"], "W3/W4 correspondence drift")
    require(correspondence["W5"]["source"] == "conj(t1)", "W5 correspondence drift")
    require(artifact["paper"]["pages"] == [29, 30], "W2--W5 paper page drift")
    require(artifact["paper"]["equations"] == [89, 93, 94, 96, 97], "W2--W5 paper equation drift")

    groups = [artifact["paper"], artifact["source"], *artifact["sources"].values()]
    for group in groups:
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"W2--W5 dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"W2--W5 dependency hash drift: {group['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "W2--W5 builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "W2--W5 checker hash drift")

    for token in (
        "Status: exact orientation ledger",
        "source conj(t2)",
        "source conj(t4)",
        "source conj(t1)",
        "Both constructions overlap on all",
        "Pi enters only through the Fourier character",
        "This is a rigorous finite exact-input audit",
        "PF-infinity, RH",
    ):
        require(token in note, f"W2--W5 note token missing: {token}")
    require("does not integrate the exact infinite nonsaddle families" in artifact["proof_boundary"], "W2--W5 boundary drift")
    print(
        "validated Hardy block-20 W2--W5 correspondence gate: "
        f"{aggregate['w3_call_count']} W3, {aggregate['w4_call_count']} W4, "
        f"{aggregate['paper_residual_excluding_zero_count']}/374 residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
