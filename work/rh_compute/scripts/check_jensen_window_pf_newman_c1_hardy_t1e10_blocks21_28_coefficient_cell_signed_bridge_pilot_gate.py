#!/usr/bin/env python3
"""Validate the later-block signed coefficient-cell bridge pilot."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate.md"
EXPECTED_CHAINS = (780, 1456)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file(), "missing later coefficient-cell bridge pilot result")
    require(NOTE.is_file(), "missing later coefficient-cell bridge pilot note")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate", "kind drift")
    require(artifact["status"] == "two_worst_ratio_later_recursive_signed_coefficient_cells_closed", "status drift")
    require(bool(artifact["passed"]), "pilot is not passed")
    require(artifact["precisions_decimal_digits"] == [70, 110], "precision ladder drift")
    require(artifact["signed_bridge_identity"] == "Dcorr=Ppi-I69pi-Qpaper=Qexact-Qpaper", "signed identity drift")
    rows = artifact["rows"]
    require(tuple(int(row["chain"]) for row in rows) == EXPECTED_CHAINS, "witness roster drift")
    require([(int(row["block"]), int(row["parent_length"])) for row in rows] == [(21, 119), (23, 154)], "witness geometry drift")
    for row in rows:
        require(row["selector"] == "W4" and row["phi3_sign"] == "negative", "witness branch drift")
        require(bool(row["radii_nonzero"]), "zero-radius cell")
        require(int(row["transport_extra_halvings"]) > 0, "cell was not transport-refined")
        require(bool(row["precision_overlap"]), "precision nonoverlap")
        require(bool(row["point_correction_overlap"]), "saved point escaped cell enclosure")
        require(bool(row["cell_correction_within_retained_majorant"]), "cell correction exceeds retained majorant")
        require(bool(row["signed_bridge_closed_on_analytic_cell"]), "signed cell bridge not closed")
        require(Fraction(row["cell_to_point_inflation_upper"]) <= Fraction(5, 4), "cell inflation exceeds cap")
        require(Fraction(row["maximum_ray_tail_upper"]) < Fraction(1, 10**30), "ray tail target drift")
        cell = row["refined_cell"]
        require(int(cell["halvings"]) == int(cell["inherited_selector_cell_halvings"]) + int(cell["transport_extra_halvings"]), "halving accounting drift")
        for radius in cell["radii"].values():
            exact = radius["exact"]
            require(Fraction(int(exact["numerator"]), int(exact["denominator"])) > 0, "nonpositive exact radius")
    aggregate = artifact["aggregate"]
    for key in ("closed_cell_count", "precision_overlap_count", "point_overlap_count", "nonzero_radius_cell_count", "within_retained_majorant_count"):
        require(int(aggregate[key]) == 2, f"aggregate {key} drift")
    require(int(aggregate["exact_algebra_checks"]) == 9, "algebra audit count drift")
    require(Fraction(aggregate["maximum_cell_to_point_inflation_upper"]) <= Fraction(5, 4), "aggregate inflation drift")
    require(bool(aggregate["transported_inflation_within_finite_margin"]), "finite margin not preserved")
    require(Fraction(aggregate["total_transported_inflation_upper"]) < Fraction(aggregate["finite_exact_point_margin"]), "transported margin arithmetic drift")
    require(int(aggregate["remaining_later_recursive_calls"]) == 1038, "remaining-call count drift")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed dependency {path}")
            require(file_hash(path) == record["sha256"], f"dependency hash drift {path}")
    note = NOTE.read_text(encoding="utf-8")
    for phrase in (
        "not a full recurrence, height-uniform theorem, or RH proof",
        "source binary128 rounding",
        "1,040 later recursive calls",
        "outer Hardy representation/remainder",
    ):
        require(phrase in note, f"proof-boundary note phrase missing: {phrase}")
    print("validated later coefficient-cell signed bridge pilot: 2/2 nonzero cells, 9 algebra checks, finite margin preserved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
