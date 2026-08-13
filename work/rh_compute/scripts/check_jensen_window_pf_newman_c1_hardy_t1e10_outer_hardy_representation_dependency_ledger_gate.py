#!/usr/bin/env python3
"""Validate the fail-closed outer-Hardy representation dependency ledger."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate as gate


FORBIDDEN_PROMOTION_KEYS = {
    "complete_hardy_error_upper",
    "hardy_error_upper",
    "xi_error_upper",
    "rh_proved",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def forbidden_keys(value: Any) -> set[str]:
    hits: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            if str(key) in FORBIDDEN_PROMOTION_KEYS:
                hits.add(str(key))
            hits.update(forbidden_keys(child))
    elif isinstance(value, list):
        for child in value:
            hits.update(forbidden_keys(child))
    return hits


def main() -> int:
    require(gate.RESULT.is_file(), "missing outer-Hardy dependency result")
    require(gate.NOTE.is_file(), "missing outer-Hardy dependency note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate",
        "kind drift",
    )
    require(
        artifact["status"] == "outer_hardy_dependency_audit_validated_with_explicit_open_remainder_columns",
        "status drift",
    )
    require(bool(artifact["passed"]), "dependency audit failed")
    require(not forbidden_keys(artifact), f"forbidden promotion keys present: {sorted(forbidden_keys(artifact))}")

    first = artifact["first_gate"]
    require(first["name"] == "complete_outer_hardy_remainder_present", "first-gate name drift")
    require(not bool(first["satisfied"]), "outer remainder was silently marked complete")
    require(not bool(first["hardy_certification_claim_permitted"]), "Hardy certification was silently permitted")

    paper = artifact["paper_audit"]
    require(paper["pdf_pages"] == [44, 45] and paper["equations"] == [126, 127], "paper location drift")
    require(bool(paper["paper_declares_hybrid_not_exact"]), "paper non-exactness boundary lost")
    require(not bool(paper["paper_supplies_explicit_uniform_constant_for_all_big_o_terms_here"]), "paper constants overpromoted")
    require(int(paper["missing_token_count"]) == 0, "paper token audit failed")

    source_rows = artifact["source_audit"]
    require(len(source_rows) == len(gate.SOURCE_ASSERTIONS), "source assertion roster drift")
    require(all(int(row["missing_token_count"]) == 0 for row in source_rows), "source assertion failed")
    require([row["id"] for row in source_rows] == [row["id"] for row in gate.SOURCE_ASSERTIONS], "source assertion order drift")

    dependencies = artifact["dependency_ledger"]
    require(len(dependencies) == 18, "dependency roster drift")
    require(len({row["id"] for row in dependencies}) == len(dependencies), "duplicate dependency id")
    require({row["state"] for row in dependencies} == {"exact", "enclosed", "open"}, "dependency state roster drift")
    require(
        {row["kind"] for row in dependencies} >= {"normalization", "conjugation", "truncation", "omitted_block", "source_arithmetic"},
        "dependency kind roster drift",
    )
    open_ids = [row["id"] for row in dependencies if row["state"] == "open"]
    require("hybrid_equations_126_127" in open_ids, "hybrid representation was silently closed")
    require("omitted_rs_corrections_and_remainder" in open_ids, "Riemann-Siegel remainder was silently closed")
    require("assembled_outer_hardy_remainder" in open_ids, "assembled outer remainder was silently closed")
    require(all(row.get("required_upgrade") for row in dependencies if row["state"] == "open"), "open row lacks required upgrade")

    outputs = artifact["output_rows"]
    require(len(outputs) == gate.EXPECTED_OUTPUTS, "output roster drift")
    for index, row in enumerate(outputs, start=1):
        require(int(row["output_index"]) == index, "output index drift")
        offset = Decimal(index - 8) / Decimal(100)
        require(row["output_label"] == f"t{offset:+.2f}", f"output label drift at {index}")
        require(Decimal(row["target_t"]) == Decimal("1e10") + offset, f"target height drift at {index}")
        require(row["open_dependency_ids"] == open_ids, f"open dependency projection drift at {index}")
        require(not bool(row["hardy_z_certified"]), f"Hardy output {index} was overpromoted")
        require(not bool(row["xi_certified"]), f"Xi output {index} was overpromoted")
        require(not bool(row["requested_tolerance_certified"]), f"tolerance output {index} was overpromoted")

    findings = artifact["critical_findings"]
    require(not bool(findings["paper_representation_exact"]), "paper representation was overpromoted")
    require(not bool(findings["source_parameter_et_is_proved_hardy_error_bound"]), "et was overpromoted")
    require(not bool(findings["source_shifted_theta_exactly_recomputed"]), "shifted theta was overpromoted")
    require(not bool(findings["source_shifted_a_exactly_recomputed"]), "shifted a was overpromoted")
    require(bool(findings["source_visibly_retains_only_leading_rs_endpoint_correction"]), "RS correction audit drift")
    require(not bool(findings["current_internal_bound_is_complete_hardy_error"]), "internal bound was overpromoted")

    aggregate = artifact["aggregate"]
    counts = aggregate["dependency_state_counts"]
    require(sum(int(value) for value in counts.values()) == len(dependencies), "dependency aggregate drift")
    require(int(aggregate["open_dependency_count"]) == len(open_ids), "open dependency aggregate drift")
    require(int(aggregate["outputs_hardy_z_certified"]) == 0, "Hardy aggregate drift")
    require(int(aggregate["outputs_xi_certified"]) == 0, "Xi aggregate drift")
    require(int(aggregate["outputs_with_explicit_open_dependencies"]) == gate.EXPECTED_OUTPUTS, "output open-boundary aggregate drift")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not exact", "not a proved error theorem", "hardy_z_certified=false", "Arb Hardy-Z calibration", "not evidence for RH"):
        require(token in note, f"note boundary token missing: {token}")

    print(
        "validated outer-Hardy dependency audit: "
        f"15 outputs, exact={counts['exact']}, enclosed={counts['enclosed']}, open={counts['open']}, Hardy-certified=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
