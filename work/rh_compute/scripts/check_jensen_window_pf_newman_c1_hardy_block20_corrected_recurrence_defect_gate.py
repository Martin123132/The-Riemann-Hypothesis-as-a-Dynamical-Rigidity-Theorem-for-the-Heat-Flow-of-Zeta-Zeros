#!/usr/bin/env python3
"""Validate the all-374 special-corrected recurrence-defect gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing corrected-defect artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate",
        "corrected-defect kind drift",
    )
    require(
        artifact["status"] == "rigorous_374_call_special_corrected_recurrence_defects_enclosed",
        "corrected-defect status drift",
    )
    scope = artifact["scope"]
    require(scope["recursive_chain_count"] == 374, "corrected-defect count drift")
    require(scope["precision_ladder_decimal_digits"] == [180, 260], "corrected-defect precision drift")
    require(scope["precision_overlap_count"] == 374, "corrected-defect overlap drift")
    require(scope["parent_upper_index"] == 104 and scope["child_upper_indices"] == [1, 2], "corrected shape drift")

    aggregate = artifact["aggregate"]
    for key in (
        "minimum_source_defect_magnitude_lower",
        "minimum_corrected_defect_magnitude_lower",
        "maximum_source_defect_magnitude_upper",
        "maximum_corrected_defect_magnitude_upper",
        "maximum_source_roundoff_gap_abs_upper",
        "maximum_defect_change_abs_upper",
        "maximum_relative_change_upper",
    ):
        require(Decimal(aggregate[key]) > 0, f"nonpositive corrected-defect bound: {key}")
    require(Decimal(aggregate["maximum_source_roundoff_gap_abs_upper"]) < Decimal("1e-20"), "corrected source replay gap too large")
    require(Decimal(aggregate["maximum_defect_change_abs_upper"]) < Decimal("1e-3"), "corrected special shift too large")
    classes = aggregate["rigorous_change_classification"]
    require(sum(classes.values()) == 374, "corrected classification count drift")
    require(aggregate["corrected_defects_excluding_zero"] == 374, "corrected zero exclusion lost")

    rows = artifact["rows"]
    require(len(rows) == 374, "corrected-defect row count drift")
    require([row["chain"] for row in rows] == sorted(row["chain"] for row in rows), "corrected row order drift")
    for row in rows:
        require(row["child_length"] in (1, 2), "corrected child length drift")
        require(row["precision_overlap"] is True, "corrected precision overlap lost")
        require(row["corrected_defect_excludes_zero"] is True, "corrected row admits zero")
        require(row["change_classification"] in classes, "corrected classification label drift")

    for key in ("telemetry", "special_residual", "selector_atlas", "point_ball_builder"):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"corrected {key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "corrected builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "corrected checker hash drift")

    for token in (
        "Status: rigorous finite exact-point recurrence defects",
        "dyadic Arb",
        "maximum special-function defect change",
        "rigorously improved calls",
        "corrected defects excluding zero",
        "does not remove the finite discrepancy",
        "remaining nonzero defect",
        "prize-level conclusion",
    ):
        require(token in note, f"corrected-defect note token missing: {token}")
    require("Nonzero defects diagnose a remaining approximation obligation" in artifact["proof_boundary"], "corrected boundary drift")
    print(
        "validated Hardy block-20 corrected recurrence-defect gate: "
        f"374 defects, 374 corrected nonzero, {classes['improved']} improved/{classes['worsened']} worsened"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
