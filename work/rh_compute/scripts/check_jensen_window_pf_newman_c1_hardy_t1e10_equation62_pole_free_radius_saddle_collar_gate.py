#!/usr/bin/env python3
"""Independently check the pole-free five-saddle collar artifact."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb

import jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing saddle-collar result")
    require(gate.NOTE.is_file(), "missing saddle-collar note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate", "kind drift")
    require(artifact["status"] == "pole_free_radius_deformation_isolates_exact_five_saddle_collar", "status drift")
    require(bool(artifact["passed"]), "saddle-collar gate failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")

    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    for output_index, row in enumerate(rows, start=1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        fresh = gate.radius_row(target_t, CHECK_PRECISION)
        require(row["output_index"] == output_index, "output index drift")
        require(row["target_t"] == fresh["target_t"], "target height drift")
        require(row["crossed_saddles"] == fresh["crossed_saddles"] == list(gate.COLLAR), "crossed saddle roster drift")
        require(row["far_pole_gap"] == fresh["far_pole_gap"], "far pole gap drift")
        require(row["branch_shell"]["certificate"] == fresh["branch_shell"]["certificate"], "branch shell certificate drift")
        for key in ("a_ball", "plus_a_radial_distance_ball", "minus_a_radial_distance_ball"):
            require(arb(row["branch_shell"][key]).overlaps(arb(fresh["branch_shell"][key])), f"branch shell ball drift: {key}")
        for selector in ("source_selector", "cell_selector"):
            require(row[selector]["classical_tail_start"] == fresh[selector]["classical_tail_start"], "tail start drift")
            for key in ("left_intercept_ball", "radius_ball", "D_lower_limit_u0_ball", "continuous_lower_root_ball", "far_intercept_ball"):
                require(arb(row[selector][key]).overlaps(arb(fresh[selector][key])), f"selector ball drift: {selector}/{key}")
        for saved, current in zip(row["saddle_rows"], fresh["saddle_rows"]):
            for key in ("n", "inside_source_selector_D_range", "inside_cell_selector_D_range", "classification"):
                require(saved[key] == current[key], f"saddle discrete drift: {key}")
            for key in ("alpha_of_n_ball", "u_saddle_ball", "b15a_equation124_threshold_identity_ball"):
                require(arb(saved[key]).overlaps(arb(current[key])), f"saddle ball drift: {key}")
            require(arb(current["b15a_equation124_threshold_identity_ball"]).contains(0), "fresh threshold identity excludes zero")

    aggregate = artifact["aggregate"]
    require(aggregate["output_count"] == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(aggregate["crossed_saddles_per_output"] == len(gate.COLLAR), "collar length aggregate drift")
    require(aggregate["total_crossed_saddle_certificates"] == gate.EXPECTED_OUTPUTS * len(gate.COLLAR), "collar certificate aggregate drift")
    require(aggregate["crossed_saddle_roster"] == list(gate.COLLAR), "collar aggregate roster drift")
    require(aggregate["common_left_pole_gap"] == [gate.TERMINAL_ALPHA, gate.NEXT_ODD_ALPHA], "left pole gap aggregate drift")
    require(aggregate["common_far_pole_gap"] == [gate.FAR_ODD_BELOW, gate.FAR_ODD_ABOVE], "far pole gap aggregate drift")

    identities = artifact["exact_identities"]
    for key in ("D_finite_interval", "saddle_threshold", "full_contour_compensation", "S_definition"):
        require(isinstance(identities[key], str) and identities[key], f"missing exact identity: {key}")
    decision = artifact["decision"]
    require(bool(decision["pole_free_radius_deformation_exists_at_all_outputs"]), "pole-free deformation lost")
    require(not bool(decision["denominator_branch_points_enter_deformation_annulus"]), "branch shell contamination")
    require(not bool(decision["alpha_pole_roster_changes_under_deformation"]), "pole roster falsely changes")
    require(bool(decision["exactly_five_D_saddles_cross_lower_limit"]), "five-saddle result lost")
    require(bool(decision["diagnostic_midpoint_defines_a_contour_resplitting_candidate"]), "resplitting candidate lost")
    require(not bool(decision["diagnostic_midpoint_is_a_published_cutoff_rule"]), "diagnostic midpoint promoted to published cutoff")
    require(not bool(decision["five_classical_terms_can_be_added_without_remainder_transport"]), "uncompensated promotion")
    require(not bool(decision["height_uniform_transported_remainder_proved"]), "height-uniform overpromotion")
    require(not bool(decision["finite_gate_has_rh_implication"]), "RH overpromotion")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file: {path}")
            require(file_hash(path) == record["sha256"], f"hash drift: {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "diagnostic midpoint", "159777,159779", "37941,37942,37943,37944,37945", "B15a", "exact compensation law", "no RH implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated pole-free radius saddle collar independently: "
        f"outputs={gate.EXPECTED_OUTPUTS}, crossed={list(gate.COLLAR)}, pole-roster-unchanged=15"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
