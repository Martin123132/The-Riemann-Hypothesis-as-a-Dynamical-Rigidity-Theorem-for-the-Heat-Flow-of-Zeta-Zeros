#!/usr/bin/env python3
"""Independently check the single-carrier adiabatic benchmark."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_single_carrier_adiabatic_benchmark"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER_PATH = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"
SUCCESSOR_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
WRONSKIAN_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json"
)
HISTORICAL_SUCCESSOR_SHA256 = (
    "c9451ba387d3f0e3f17a914719fc1ec899cfd41bf1667beb19073d1b89c29e90"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> int:
    issues: list[str] = []
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")

    if artifact.get("kind") != STEM:
        issues.append("kind drifted")
    if artifact.get("builder_sha256") != file_hash(BUILDER_PATH):
        issues.append("builder hash mismatch")
    stored_sources = artifact.get("source_sha256", {})
    current_successor_hash = file_hash(SUCCESSOR_RESULT)
    if stored_sources.get("successor_lemma") not in {
        current_successor_hash,
        HISTORICAL_SUCCESSOR_SHA256,
    }:
        issues.append("successor source snapshot is unrecognized")
    if stored_sources.get("critical_component_wronskian_gate") != (
        file_hash(WRONSKIAN_RESULT)
    ):
        issues.append("Wronskian source hash contract mismatch")
    successor = json.loads(SUCCESSOR_RESULT.read_text(encoding="utf-8"))
    if successor.get("kind") != (
        "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"
    ):
        issues.append("current successor semantic source kind drifted")
    if successor.get("exact", {}).get("shell", {}).get("parameters") != (
        "t_j=1/(5*j), R_j=j+38"
    ):
        issues.append("current successor semantic shell drifted")
    if "antecedent open" not in successor.get("status", ""):
        issues.append("current successor semantic guard drifted")

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    expected_ids = [
        "single_carrier_first_jet",
        "phase_matrix_determinant",
        "carrier_clearance_floor",
        "carrier_heat_jet",
        "uniform_log_jet_hypotheses",
        "adiabatic_condition_number",
        "cofinal_shell_cost",
        "gronwall_transport",
        "four_carrier_obstruction",
        "xi_crossing_small_ball_handoff",
    ]
    if ids != expected_ids:
        issues.append("row ids or order drifted")

    u, v, ell = sp.symbols("u v ell", real=True, nonzero=True)
    matrix = sp.Matrix([[1, 0], [u / ell, -v / ell]])
    if sp.simplify(matrix.det() + v / ell) != 0:
        issues.append("phase-matrix determinant identity failed")
    frobenius_sq = sum(entry**2 for entry in matrix)
    expected_frobenius = 1 + (u**2 + v**2) / ell**2
    if sp.simplify(frobenius_sq - expected_frobenius) != 0:
        issues.append("phase-matrix Frobenius identity failed")

    r, r_x, r_xx = sp.symbols("r r_x r_xx")
    second_log_jet = r**2 + r_x
    third_from_recurrence = sp.expand(
        r * second_log_jet + 2 * r * r_x + r_xx
    )
    expected_third = r**3 + 3 * r * r_x + r_xx
    if sp.simplify(third_from_recurrence - expected_third) != 0:
        issues.append("third logarithmic jet identity failed")

    C1, C2, C3, c = sp.symbols(
        "C1 C2 C3 c",
        positive=True,
    )
    B2 = C1**2 + C2
    B3 = C1**3 + 3 * C1 * C2 + C3
    K = sp.sqrt(1 + C1**2) * sp.sqrt(B2**2 + B3**2) / c
    if K.has(sp.zoo, sp.nan):
        issues.append("condition-number constant is nonfinite")

    j = sp.symbols("j", positive=True)
    shell_cost = (
        sp.log((j + 38) / (4 * sp.pi)) ** 2
        / (5 * j * (j + 1))
    )
    if sp.limit(shell_cost, j, sp.oo) != 0:
        issues.append("cofinal delta*L^2 limit is not zero")

    x = sp.symbols("x", real=True)
    I = sp.I
    toy = (
        sp.exp(-3 * I * x)
        - sp.exp(-4 * I * x)
        + I * sp.exp(-I * x)
        - I * sp.exp(-2 * I * x) / 2
    )
    toy_0 = sp.simplify(toy.subs(x, 0))
    toy_1 = sp.simplify(sp.diff(toy, x).subs(x, 0))
    toy_2_real = sp.simplify(
        sp.re(sp.diff(toy, x, 2).subs(x, 0))
    )
    toy_rate_imag = sp.simplify(sp.im(toy_1 / toy_0))
    if toy_0 != I / 2:
        issues.append(f"toy E(0)={toy_0}, expected i/2")
    if toy_1 != I:
        issues.append(f"toy E'(0)={toy_1}, expected i")
    if toy_2_real != 7:
        issues.append(
            f"toy Re(E)''(0)={toy_2_real}, expected 7"
        )
    if toy_rate_imag != 0:
        issues.append("toy effective phase speed is not zero")

    summary = artifact.get("exact_summary", {})
    expected_summary = {
        "row_count": 10,
        "exact_benchmark_reductions": 8,
        "four_carrier_obstructions": 1,
        "open_xi_small_ball_handoffs": 1,
        "cofinal_cost": "delta_j*ell_j^2->0",
        "generic_multi_carrier_promotion": False,
        "all_j_xi_theorem": False,
    }
    if summary != expected_summary:
        issues.append("exact summary drifted")
    open_row = rows[-1] if rows else {}
    if open_row.get("status") != "open and not proved":
        issues.append("Xi small-ball handoff is not explicitly open")

    required_note = (
        "delta_j*ell_j^2 -> 0",
        "Why It Is Not Xi Yet",
        "Im(E_toy'(0)/E_toy(0))=0",
        "crossing-restricted",
        "remain open",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note missing marker {marker!r}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "not proved",
        "does not prove an all-j",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof boundary missing {marker!r}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman single-carrier adiabatic benchmark: "
        "10 rows, 8 exact benchmark reductions, "
        "1 four-carrier obstruction, 1 open Xi small-ball handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
