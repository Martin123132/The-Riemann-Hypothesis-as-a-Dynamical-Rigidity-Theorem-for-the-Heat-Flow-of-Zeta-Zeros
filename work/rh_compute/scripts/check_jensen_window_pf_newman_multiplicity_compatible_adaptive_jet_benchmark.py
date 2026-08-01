#!/usr/bin/env python3
"""Independently check the multiplicity-compatible adaptive-jet benchmark."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER_PATH = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"
SUCCESSOR_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
SINGLE_CARRIER_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_single_carrier_adiabatic_benchmark.json"
)
CROSSING_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.json"
)
ATTAINMENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
DELTA_LOCALIZATION_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_positive_boundary_delta_localization_gate.json"
)
HISTORICAL_SUCCESSOR_SHA256 = (
    "c9451ba387d3f0e3f17a914719fc1ec899cfd41bf1667beb19073d1b89c29e90"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def check_hermite_identities(issues: list[str]) -> None:
    t, y = sp.symbols("t y", positive=True, real=True)
    z = sp.symbols("z", real=True)
    s = sp.sqrt(2 * t)
    for m in range(1, 13):
        heat_polynomial = sp.factorial(m) * sum(
            (-t) ** k
            * y ** (m - 2 * k)
            / (
                sp.factorial(k)
                * sp.factorial(m - 2 * k)
            )
            for k in range(m // 2 + 1)
        )
        hermite_form = sp.expand(
            s**m * sp.hermite_prob(m, y / s)
        )
        if sp.simplify(heat_polynomial - hermite_form) != 0:
            issues.append(f"m={m}: Hermite rescaling failed")
        if sp.simplify(
            sp.diff(heat_polynomial, t)
            + sp.diff(heat_polynomial, y, 2)
        ) != 0:
            issues.append(f"m={m}: backward heat equation failed")

        first = sp.diff(heat_polynomial, y)
        if sp.simplify(
            first
            - m * s ** (m - 1) * sp.hermite_prob(m - 1, y / s)
        ) != 0:
            issues.append(f"m={m}: first Hermite jet failed")

        adaptive_second = s * first
        derivative_vector = (
            sp.diff(heat_polynomial, t),
            sp.diff(adaptive_second, t),
        )
        he_m2 = (
            sp.hermite_prob(m - 2, z) if m >= 2 else sp.Integer(0)
        )
        he_m3 = (
            sp.hermite_prob(m - 3, z) if m >= 3 else sp.Integer(0)
        )
        expected_a = -m * (m - 1) * he_m2
        expected_b = (
            m * sp.hermite_prob(m - 1, z)
            - m * (m - 1) * (m - 2) * he_m3
        )
        substitutions = {y: s * z}
        scaled_actual = tuple(
            sp.simplify(value.subs(substitutions) / s ** (m - 2))
            for value in derivative_vector
        )
        if scaled_actual != (
            sp.simplify(expected_a),
            sp.simplify(expected_b),
        ):
            issues.append(f"m={m}: adaptive heat jet failed")

        current = sp.Poly(sp.hermite_prob(m, z), z)
        previous = sp.Poly(sp.hermite_prob(m - 1, z), z)
        if sp.degree(sp.gcd(current, previous)) != 0:
            issues.append(f"m={m}: consecutive Hermite gcd is nontrivial")


def main() -> int:
    issues: list[str] = []
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")

    if artifact.get("kind") != STEM:
        issues.append("kind drifted")
    if artifact.get("builder_sha256") != file_hash(BUILDER_PATH):
        issues.append("builder hash mismatch")
    expected_sources = {
        "single_carrier_adiabatic_benchmark": file_hash(
            SINGLE_CARRIER_RESULT
        ),
        "crossing_slope_gap_reduction": file_hash(CROSSING_RESULT),
        "positive_boundary_attainment_lemma": file_hash(
            ATTAINMENT_RESULT
        ),
        "positive_boundary_delta_localization_gate": file_hash(
            DELTA_LOCALIZATION_RESULT
        ),
    }
    stored_sources = artifact.get("source_sha256", {})
    for key, expected in expected_sources.items():
        if stored_sources.get(key) != expected:
            issues.append(f"source hash contract mismatch: {key}")
    current_successor_hash = file_hash(SUCCESSOR_RESULT)
    if stored_sources.get("adiabatic_successor_lemma") not in {
        current_successor_hash,
        HISTORICAL_SUCCESSOR_SHA256,
    }:
        issues.append("successor source snapshot is unrecognized")
    successor = json.loads(SUCCESSOR_RESULT.read_text(encoding="utf-8"))
    if successor.get("kind") != (
        "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"
    ):
        issues.append("current successor semantic source kind drifted")
    if "antecedent open" not in successor.get("status", ""):
        issues.append("current successor semantic guard drifted")
    if successor.get("exact", {}).get("heat_jet", {}).get(
        "time_derivative_norm_squared"
    ) != "h_xx**2 + h_xxx**2/ell**2":
        issues.append("current successor semantic heat-jet drifted")

    expected_ids = [
        "fixed_x_uniform_floor_guard",
        "backward_heat_monomial",
        "hermite_rescaling",
        "real_simple_splitting",
        "adaptive_first_jet",
        "adaptive_heat_jet",
        "multiplicity_condition_number",
        "cofinal_relative_transport",
        "time_dependent_scaling_transfer",
        "xi_domain_split_target",
    ]
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row ids or order drifted")
    if not rows or rows[-1].get("status") != "open and not proved":
        issues.append("two-regime Xi target is not explicitly open")

    check_hermite_identities(issues)

    t = sp.symbols("t", positive=True)
    s = sp.sqrt(2 * t)
    if sp.simplify(sp.diff(s, t) - 1 / s) != 0:
        issues.append("adaptive scale derivative failed")
    j = sp.symbols("j", positive=True)
    if sp.limit(sp.log(1 + 1 / j), j, sp.oo) != 0:
        issues.append("cofinal logarithmic step cost failed")

    expected_summary = {
        "row_count": 10,
        "exact_benchmark_reductions": 9,
        "open_xi_targets": 1,
        "adaptive_condition_number": "C_m/(2t)",
        "cofinal_step_cost": "K_m*log(1+1/j)->0",
        "requires_endpoint_simplicity": False,
        "all_j_xi_theorem": False,
    }
    if artifact.get("exact_summary") != expected_summary:
        issues.append("exact summary drifted")

    required_note = (
        "Why The Uniform Floor Is Too Strong",
        "existing positive-boundary attainment",
        "new contribution",
        "Multiple-Zero Thought Experiment",
        "P_m(t,y)=s^m He_m(z)",
        "C_m/(2t)",
        "(j/(j+1))^K_m",
        "two-regime theorem",
        "compatible with finite-multiplicity zeros",
        "not a local factorization theorem for Xi",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note missing marker {marker!r}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "No Xi",
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
        "validated Newman multiplicity-compatible adaptive-jet benchmark: "
        "10 rows, 9 exact reductions, 1 open two-regime Xi handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
