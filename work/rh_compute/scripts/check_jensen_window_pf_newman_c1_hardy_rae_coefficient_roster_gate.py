#!/usr/bin/env python3
"""Validate the coefficient-to-physical-RAE roster inversion gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.py"
)
CHECKER = Path(__file__).resolve()
SOURCE = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
)
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER, SOURCE):
        require(path.is_file(), f"missing RAE-map artifact: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")

    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate",
        "RAE-map kind drift",
    )
    require(
        artifact["status"]
        == "rigorous_fixture_physical_roster_inversion_with_continuous_coverage_open",
        "RAE-map status drift",
    )
    require(artifact["source"]["sha256"] == SOURCE_SHA256, "RAE-map source hash drift")
    require(len(artifact["source"]["locations"]) == 19, "source location inventory drift")
    provenance = artifact["source"]["pi_provenance"]
    require(provenance["ordinary_real_constant"] == "p=4*atan(1)", "pi origin drift")
    require("Pi2n(1)=2*pi" in provenance["physical_scale_constant"], "PARI pi origin drift")

    theorem = artifact["exact_inverse_theorem"]
    require(theorem["branch_cancellation"] == "xx=2*a2-3*a3 on both source branches", "branch inverse drift")
    require(theorem["r_from_pc"] == "r=(sqrt(pc)+1/sqrt(pc))/2", "physical inverse drift")
    require(theorem["exact_rational_fixture_count"] == 8, "exact fixture count drift")
    require(len(theorem["exact_rational_fixtures"]) == 8, "exact fixture roster drift")

    fixture = artifact["fixture"]
    require(fixture["block"] == 20, "fixture block drift")
    require(fixture["rn1_start"] == 656853, "fixture RN1 drift")
    require(fixture["mt"] == 209 and fixture["aenums"] == 212, "fixture controls drift")
    require(fixture["first_rae"] == 657064, "first physical pivot drift")
    require(fixture["spacing"] == 420, "physical roster spacing drift")
    require(fixture["full_last_rae"] == 745684, "full block endpoint drift")
    require(fixture["last_telemetry_rae"] == 670084, "telemetry endpoint drift")
    require(
        fixture["rn1_end"] == fixture["rn1_start"] + fixture["spacing"] * fixture["aenums"],
        "checkpoint block increment drift",
    )

    ladder = artifact["precision_ladder"]
    require(ladder["bits"] == [192, 320], "precision ladder drift")
    require(ladder["stable_linear_modulo_lifts"] is True, "modulo lift instability")
    require(ladder["stable_recovered_roster"] is True, "physical roster instability")

    aggregate = artifact["aggregate"]
    expected_roster = [fixture["first_rae"] + fixture["spacing"] * index for index in range(32)]
    require(aggregate["pair_count"] == 32 and aggregate["branch_count"] == 64, "aggregate count drift")
    require(aggregate["recovered_roster"] == expected_roster, "recovered roster drift")
    require(
        Decimal(aggregate["maximum_coefficient_replay_gap"]["ball"]["upper_decimal"])
        < Decimal(2) ** Decimal(-88),
        "coefficient replay tolerance failed",
    )
    require(
        Decimal(aggregate["maximum_recovered_rae_error"]["ball"]["upper_decimal"])
        < Decimal(2) ** Decimal(-55),
        "RAE recovery tolerance failed",
    )
    require(
        Decimal(aggregate["minimum_linear_modulo_margin"]["ball"]["lower_decimal"]) > 0,
        "linear modulo margin failed",
    )

    rows = artifact["rows"]
    require(len(rows) == 32, "RAE-map row count drift")
    for index, row in enumerate(rows, 1):
        target = expected_roster[index - 1]
        require(row["sum_index"] == index and row["target_rae"] == target, "physical row sequence drift")
        require(len(row["branches"]) == 2, f"sum index {index} branch count drift")
        require([branch["branch"] for branch in row["branches"]] == [1, 2], "branch ordering drift")
        for branch in row["branches"]:
            require(branch["nearest_physical_pivot"] == target, "nearest physical pivot drift")
            require(len(branch["initial_coefficients_hex"]) == 3, "saved coefficient degree drift")
            require(
                Decimal(branch["recovered_rae_error"]["upper_decimal"]) < Decimal("0.5"),
                "physical pivot uniqueness drift",
            )
            require(
                all(Decimal(gap["upper_decimal"]) < Decimal(2) ** Decimal(-88) for gap in branch["ideal_formula_replay_gaps"].values()),
                "row replay gap drift",
            )

    handoff = artifact["next_handoff"]
    require("Interval-evaluate" in handoff["first"], "coverage handoff drift")
    require("all 212" in handoff["second"], "full roster handoff drift")
    require("failed coverage claim" in handoff["falsification_rule"], "falsification rule drift")

    for key in ("telemetry", "checkpoint", "input"):
        dependency = REPO_ROOT / artifact["sources"][key]["path"]
        require(file_hash(dependency) == artifact["sources"][key]["sha256"], f"{key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: exact inversion and rigorous 32-pivot fixture recovery",
        "xx = 2*a2 - 3*a3",
        "r  = (sqrt(pc) + 1/sqrt(pc))/2",
        "p=4*atan(1)",
        "a = sqrt(8*t/pi)",
        "first `32` values",
        "does **not** yet prove that an interval",
        "remaining 180 pivots",
        "prize-level conclusion",
    ):
        require(token in note, f"RAE-map note token missing: {token}")

    require("continuous q-cell coverage" in artifact["proof_boundary"], "coverage boundary missing")
    print(
        "validated Hardy coefficient-to-RAE roster gate: "
        f"{aggregate['pair_count']} pivots, {aggregate['branch_count']} branches, "
        f"spacing {fixture['spacing']}, continuous coverage open"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
