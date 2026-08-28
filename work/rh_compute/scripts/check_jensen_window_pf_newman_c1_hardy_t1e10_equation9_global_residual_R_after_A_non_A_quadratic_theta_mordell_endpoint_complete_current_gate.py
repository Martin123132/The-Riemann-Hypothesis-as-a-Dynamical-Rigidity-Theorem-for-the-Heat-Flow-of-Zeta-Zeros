#!/usr/bin/env python3
"""Independently check the endpoint-complete Mordell current gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
L = 2_481_422
K = L - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "Mordell current artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    require(artifact["scope"]["L"] == L and artifact["scope"]["K"] == K, "source roster drift")

    a, z, tau, r, k = sp.symbols("a z tau r k", real=True)
    g0, gw = sp.symbols("G G_w")
    c_factor = sp.exp(sp.I * sp.pi / 4 - sp.I * sp.pi * z**2 / (2 * tau)) / sp.sqrt(2 * tau)
    main = a * c_factor * g0 / tau + (
        sp.diff(c_factor, z) * g0 + c_factor * gw / (2 * tau)
    ) / (sp.pi * sp.I)
    expected = c_factor * ((a - z) * g0 + gw / (2 * sp.pi * sp.I)) / tau
    require(sp.simplify(main - expected) == 0, "independent main-current cancellation failed")
    phase_left = (a**2 - z**2) / (2 * tau) + z * k / tau - k**2 / (2 * tau)
    phase_right = a * (r + k) / tau - (r + k) ** 2 / (2 * tau)
    require(sp.simplify((phase_left - phase_right).subs(a, z + r)) == 0, "independent dual phase failed")

    decisions = artifact["decision"]
    for key in (
        "Kuznetsov_Theorem_1_specialized_to_non_A_current",
        "normalized_shift_preserves_absolute_dual_first_moment",
        "transformed_zeroth_and_first_moments_joined_exactly",
        "two_Mordell_endpoint_terms_joined_exactly_with_derivatives",
        "endpoint_complete_one_step_identity_built",
        "full_length_one_step_cross_checked_against_reference_oracle",
        "Mordell_quadrature_rows_are_diagnostic_only",
    ):
        require(decisions.get(key) is True, f"missing exact/diagnostic decision: {key}")
    for key in (
        "certified_Mordell_h_and_h_prime_evaluator_built",
        "small_tau_Euler_Maclaurin_branch_built",
        "recursive_interval_theta_evaluator_built",
        "physical_quadrature_completed",
        "non_A_bound_proved",
        "R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    short_rows = artifact["short_high_precision_rows"]
    require(len(short_rows) == 2 and all(row["passed"] is True for row in short_rows), "short roster drift")
    for row in short_rows:
        require(float(row["absolute_error"]) <= float(row["allowance"]), "short row exceeds allowance")
        z_value = Fraction(row["normalized_z"])
        require(Fraction(-1, 2) <= z_value < Fraction(1, 2), "short normalized z drift")
        require(row["transformed_m"] == floor_fraction(2 * row["n"] * Fraction(row["tau"])), "short m drift")

    full_rows = artifact["full_one_step_rows"]
    require(len(full_rows) == 2 and all(row["diagnostic_only"] is True for row in full_rows), "full roster drift")
    for row in full_rows:
        x = Fraction(row["x"])
        s = Fraction(row["s"])
        tau_value = x / 2
        a_value = x * (A + 2 * s) / 2
        shift = floor_fraction(a_value + Fraction(1, 2))
        require(row["integer_shift_r"] == shift, "full phase shift drift")
        require(Fraction(row["normalized_z"]) == a_value - shift, "full normalized phase drift")
        require(row["transformed_m"] == floor_fraction(2 * K * tau_value), "full m drift")
        require(row["transformed_term_count"] == row["transformed_m"] + 1, "full transformed count drift")
        require(row["transformed_mass_normalized_error"] <= 2.0e-10, "transformed diagnostic drift")
        require(row["complete_mass_normalized_error"] <= 2.0e-10, "complete diagnostic drift")
        period = row["minimal_transformed_period"]
        w = Fraction(row["transformed_w"])
        sigma = Fraction(row["transformed_sigma"])
        require((2 * sigma * period).denominator == 1, "period slope condition failed")
        require((w * period + sigma * period * period).denominator == 1, "period constant condition failed")

    source = artifact["primary_source"]
    require(source["arxiv"] == "1306.4081v2", "primary source drift")
    require(source["imported_numerical_constants"] is False, "source constant boundary drift")
    require(source["practical_quadrature_rigorous_error_bound_claimed_by_source"] is False, "quadrature boundary drift")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("(ME1)", "(ME3)", "(ME4)", "(ME6)", "r+k", "h_z"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("not interval enclosures" in normalized and "rh" in normalized, "proof boundary missing")
    print("independently checked endpoint-complete Mordell one-step current identity", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
