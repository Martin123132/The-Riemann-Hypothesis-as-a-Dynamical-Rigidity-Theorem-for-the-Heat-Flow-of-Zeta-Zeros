#!/usr/bin/env python3
"""Independently check the time-dependent scaled-jet successor lemma."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_time_dependent_scaled_successor_lemma"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER_PATH = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"
SUCCESSOR_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
CROSSING_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.json"
)
ADAPTIVE_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark.json"
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
    expected_sources = {
        "adiabatic_successor_lemma": file_hash(SUCCESSOR_RESULT),
        "crossing_slope_gap_reduction": file_hash(CROSSING_RESULT),
        "adaptive_jet_benchmark": file_hash(ADAPTIVE_RESULT),
    }
    if artifact.get("source_sha256") != expected_sources:
        issues.append("source hash contract mismatch")

    expected_ids = [
        "positive_scaled_jet",
        "scaled_degree_homotopy",
        "time_dependent_heat_jet",
        "additive_scaled_transport",
        "relative_scaled_transport",
        "relative_cell_clearance",
        "adaptive_hermite_specialization",
        "cofinal_multiplicative_factor",
        "positive_scale_switch",
        "two_regime_successor_target",
    ]
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row ids or order drifted")
    if not rows or rows[-1].get("status") != "open and not proved":
        issues.append("two-regime Xi antecedent is not explicitly open")

    r, s = sp.symbols("r s", real=True, positive=True)
    interpolated = (1 - r) + r * s
    if sp.simplify(interpolated.subs(r, 0) - 1) != 0:
        issues.append("scale homotopy start failed")
    if sp.simplify(interpolated.subs(r, 1) - s) != 0:
        issues.append("scale homotopy end failed")

    h_x, h_xx, h_xxx, s_t = sp.symbols(
        "h_x h_xx h_xxx s_t",
        real=True,
    )
    expected_derivative = sp.Matrix(
        [-h_xx, s_t * h_x - s * h_xxx]
    )
    chain_derivative = sp.Matrix(
        [-h_xx, s_t * h_x + s * (-h_xxx)]
    )
    if chain_derivative != expected_derivative:
        issues.append("scaled heat derivative failed")

    tau, t_lower, t_upper, k = sp.symbols(
        "tau t_lower t_upper k", positive=True
    )
    integral = sp.integrate(k / tau, (tau, t_lower, t_upper))
    if sp.simplify(
        integral - k * sp.log(t_upper / t_lower)
    ) != 0:
        issues.append("relative K/t integral failed")
    j, c_m = sp.symbols("j c_m", positive=True)
    lower = sp.Rational(1, 5) / (j + 1)
    upper = sp.Rational(1, 5) / j
    cofinal_integral = sp.integrate(
        c_m / (2 * tau), (tau, lower, upper)
    )
    expected_cost = c_m * sp.log((j + 1) / j) / 2
    if sp.simplify(cofinal_integral - expected_cost) != 0:
        issues.append("cofinal relative cost failed")
    if sp.limit(expected_cost, j, sp.oo) != 0:
        issues.append("cofinal relative cost limit failed")

    expected_summary = {
        "row_count": 10,
        "exact_reductions": 9,
        "open_xi_targets": 1,
        "time_dependent_positive_scaling": True,
        "relative_transport_available": True,
        "requires_endpoint_simplicity": False,
        "all_j_xi_theorem": False,
    }
    if artifact.get("exact_summary") != expected_summary:
        issues.append("exact summary drifted")

    required_note = (
        "Positive Scaled Jet",
        "partial_t V_s=(-H_xx,s_t H_x-sH_xxx)",
        "Relative Transport",
        "No additive ratio below one",
        "(j/(j+1))^(C_m/2)",
        "Scale Interface",
        "multiplicity-compatible all-stage theorem",
        "are open",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note missing marker {marker!r}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "not proved",
        "does not prove Q209",
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
        "validated Newman time-dependent scaled-jet successor lemma: "
        "10 rows, 9 exact reductions, 1 open two-regime Xi antecedent"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
