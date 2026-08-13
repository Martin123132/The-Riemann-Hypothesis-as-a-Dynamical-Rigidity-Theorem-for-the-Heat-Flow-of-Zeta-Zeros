#!/usr/bin/env python3
"""Independently check the certified five-saddle D transport."""

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

from flint import acb, arb

import jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing five-saddle D result")
    require(gate.NOTE.is_file(), "missing five-saddle D note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate", "kind drift")
    require(artifact["status"] == "five_saddle_D_transport_certified_and_residual_compensation_target_isolated", "status drift")
    require(bool(artifact["passed"]), "five-saddle D gate failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")
    require(artifact["scope"]["collar"] == list(gate.COLLAR), "collar roster drift")
    algebra = gate.phase_identity_audit()
    require(bool(algebra["identity_exact"]), "fresh stable phase identity failed")
    require(algebra == artifact["stable_phase_identity"], "saved phase identity drift")

    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    gaps: list[arb] = []
    ratios: list[arb] = []
    phase_differences: list[arb] = []
    for output_index, row in enumerate(rows, start=1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        fresh = gate.output_transport(target_t, CHECK_PRECISION)
        require(row["output_index"] == output_index, "output index drift")
        require(row["target_t"] == fresh["target_t"], "target height drift")
        for key in (
            "source_D_lower_limit_ball",
            "cell_D_lower_limit_ball",
            "D_interval_width_ball",
            "alternating_D_sum_real_ball",
            "alternating_D_sum_imag_ball",
            "transformed_D_collar_real_ball",
            "transformed_D_collar_imag_ball",
            "paper_phase_classical_collar_ball",
            "exact_theta_classical_collar_ball",
            "D_real_minus_paper_classical_ball",
            "D_to_paper_classical_absolute_ratio_ball",
            "exact_theta_minus_paper_phase_collar_ball",
        ):
            require(arb(row[key]).overlaps(arb(fresh[key])), f"fresh output ball drift: {key}")
        require(len(row["mode_rows"]) == len(gate.COLLAR), "mode roster length drift")
        for saved, current in zip(row["mode_rows"], fresh["mode_rows"]):
            require(saved["n"] == current["n"], "mode n drift")
            for key in (
                "D_source_minus_cell_real_ball",
                "D_source_minus_cell_imag_ball",
                "B18_Y_main_real_ball",
                "B18_Y_main_imag_ball",
                "D_difference_minus_Y_real_ball",
                "D_difference_minus_Y_imag_ball",
                "D_to_Y_ratio_real_ball",
                "D_to_Y_ratio_imag_ball",
                "q_at_source_limit_real_ball",
                "q_at_source_limit_imag_ball",
                "q_at_cell_limit_real_ball",
                "q_at_cell_limit_imag_ball",
            ):
                require(arb(saved[key]).overlaps(arb(current[key])), f"fresh mode ball drift at n={saved['n']}: {key}")
        gap = arb(fresh["D_real_minus_paper_classical_ball"])
        ratio = arb(fresh["D_to_paper_classical_absolute_ratio_ball"])
        phase_difference = abs(arb(fresh["exact_theta_minus_paper_phase_collar_ball"]))
        require(gap > gate.REQUESTED_TOLERANCE, f"fresh transport gap failed at output={output_index}")
        require(arb("0.6") < ratio < arb("0.8"), f"fresh transport ratio failed at output={output_index}")
        gaps.append(gap)
        ratios.append(ratio)
        phase_differences.append(phase_difference)

    aggregate = artifact["aggregate"]
    require(aggregate["output_count"] == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(aggregate["D_integrals_per_output"] == len(gate.COLLAR), "integrals/output aggregate drift")
    require(aggregate["total_certified_D_integrals"] == gate.EXPECTED_OUTPUTS * len(gate.COLLAR), "integral aggregate drift")
    require(aggregate["D_transport_gap_above_0p005_count"] == gate.EXPECTED_OUTPUTS, "gap count aggregate drift")
    require(min(gaps, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_D_real_minus_classical_ball"])), "minimum gap aggregate drift")
    require(max(gaps, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_D_real_minus_classical_ball"])), "maximum gap aggregate drift")
    require(min(ratios, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_D_to_classical_absolute_ratio_ball"])), "minimum ratio aggregate drift")
    require(max(ratios, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_D_to_classical_absolute_ratio_ball"])), "maximum ratio aggregate drift")
    require(max(phase_differences, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_exact_theta_paper_phase_difference_absolute_ball"])), "phase difference aggregate drift")

    decision = artifact["decision"]
    require(bool(decision["five_D_integrals_certified_at_all_outputs"]), "D integration result lost")
    require(not bool(decision["five_D_transport_equals_five_classical_terms"]), "false D/classical equality")
    require(bool(decision["D_transport_gap_exceeds_0p005_at_all_outputs"]), "D gap result lost")
    require(not bool(decision["exact_theta_phase_refinement_closes_gap"]), "phase refinement falsely closes gap")
    require(bool(decision["remaining_channels_must_be_transported_jointly"]), "joint transport obligation lost")
    require(bool(decision["transport_belongs_to_diagnostic_alternative_cutoff"]), "diagnostic transport classification lost")
    require(not bool(decision["diagnostic_cutoff_is_published_or_validated"]), "diagnostic cutoff overpromoted")
    require(not bool(decision["height_uniform_transported_remainder_proved"]), "height-uniform overpromotion")
    require(not bool(decision["finite_gate_has_rh_implication"]), "RH overpromotion")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file: {path}")
            require(file_hash(path) == record["sha256"], f"hash drift: {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "diagnostic", "Lewis", "40-term", "37941,...,37945", "strictly above `0.005`", "cannot simply be", "no RH implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated five-saddle D transport independently: "
        f"integrals={aggregate['total_certified_D_integrals']}, gap>0.005={len(gaps)}/{gate.EXPECTED_OUTPUTS}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
