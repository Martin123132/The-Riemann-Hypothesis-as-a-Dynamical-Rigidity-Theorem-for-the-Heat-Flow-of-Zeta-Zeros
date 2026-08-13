#!/usr/bin/env python3
"""Validate the t5-erfc-corrected block-20 recurrence-defect gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing t5-erfc corrected-defect artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate",
        "t5-erfc corrected-defect kind drift",
    )
    require(
        artifact["status"]
        == "rigorous_374_call_special_and_intrinsic_erfc_corrected_recurrence_defects_enclosed",
        "t5-erfc corrected-defect status drift",
    )
    scope = artifact["scope"]
    require(scope["recursive_chain_count"] == 374, "t5-erfc corrected chain count drift")
    require(scope["intrinsic_real_erfc_evaluation_count"] == 913, "t5-erfc corrected call count drift")
    require(scope["precision_ladder_decimal_digits"] == [180, 260], "t5-erfc corrected precision drift")
    require(scope["precision_overlap_count"] == 374, "t5-erfc precision-overlap drift")
    require(scope["prior_special_corrected_defect_overlap_count"] == 374, "prior corrected-defect overlap drift")
    require(scope["parent_upper_index"] == 104, "t5-erfc parent index drift")
    require(scope["child_upper_indices"] == [1, 2], "t5-erfc child indices drift")

    aggregate = artifact["aggregate"]
    require(aggregate["full_corrected_defects_excluding_zero"] == 374, "full corrected zero exclusion drift")
    require(Decimal(aggregate["minimum_full_corrected_defect_magnitude_lower"]) > 0, "nonpositive full defect minimum")
    require(Decimal(aggregate["maximum_full_corrected_defect_magnitude_upper"]) > 0, "nonpositive full defect maximum")
    require(Decimal(aggregate["maximum_intrinsic_erfc_defect_increment_abs_upper"]) < Decimal("1e-30"), "intrinsic erfc increment unexpectedly large")
    require(Decimal(aggregate["maximum_intrinsic_erfc_relative_change_upper"]) < Decimal("1e-28"), "intrinsic erfc relative change unexpectedly large")
    require(Decimal(aggregate["minimum_global_scale_separation_lower"]) > Decimal("1e30"), "t5-erfc scale separation lost")
    classes = aggregate["intrinsic_erfc_change_classification"]
    require(sum(classes.values()) == 374, "t5-erfc change classification drift")

    rows = artifact["rows"]
    require(len(rows) == 374, "t5-erfc corrected row count drift")
    require([row["chain"] for row in rows] == sorted(row["chain"] for row in rows), "t5-erfc corrected order drift")
    for row in rows:
        require(row["precision_overlap"], "t5-erfc precision nonoverlap row")
        require(row["prior_special_corrected_defect_overlap"], "prior corrected-defect nonoverlap row")
        require(row["full_corrected_defect_excludes_zero"], "full corrected defect contains zero")
        require(Decimal(row["intrinsic_real_erfc_defect_increment_abs_upper"]) >= 0, "negative erfc increment bound")
        require(row["change_classification"] in {"improved", "worsened", "indeterminate"}, "change class drift")

    for key in (
        "telemetry",
        "special_residual",
        "t5_real_erfc_residual",
        "prior_corrected_defect",
        "selector_atlas",
        "point_ball_builder",
    ):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"t5-erfc corrected {key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "t5-erfc corrected builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "t5-erfc corrected checker hash drift")

    for token in (
        "Status: rigorous finite recurrence propagation",
        "All imported correction balls use their stored dyadic endpoints",
        "maximum intrinsic-real-erfc defect increment",
        "over thirty orders of magnitude",
        "endpoint-saddle, finite-`ip`",
        "prize-level theorem",
    ):
        require(token in note, f"t5-erfc corrected note token missing: {token}")
    require("identify but do not derive or bound" in artifact["proof_boundary"], "t5-erfc corrected boundary drift")
    print(
        "validated Hardy block-20 t5-erfc-corrected recurrence-defect gate: "
        f"374 nonzero defects, scale separation {aggregate['minimum_global_scale_separation_lower']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
