#!/usr/bin/env python3
"""Independently check the quadratic-theta Poisson current self-duality."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate"
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


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "self-duality artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    require(artifact["scope"]["L"] == L and artifact["scope"]["K"] == K, "theta roster drift")
    require(artifact["scope"]["Hiary_moments"] == [0, 1], "moment roster drift")

    n, v, x, c = sp.symbols("n v x c", positive=True, real=True)
    phase = x * (c + 2 * n) ** 2 / 4 - 2 * v * n
    n_star = v / x - c / 2
    require(sp.simplify(sp.diff(phase, n).subs(n, n_star)) == 0, "independent stationary point failed")
    require(sp.simplify((c + 2 * n_star) - 2 * v / x) == 0, "independent current transport failed")
    require(sp.simplify(phase.subs(n, n_star) - c * v + v**2 / x) == 0, "independent dual phase failed")

    mapping = artifact["theta_mapping"]
    require(mapping["parameters"] == "K=L-1, c=A+2s, a=x*c/2, b=x/2", "parameter map drift")
    require(mapping["S0"].startswith("S_0=F(K,0"), "S0 map drift")
    require(mapping["S1"].startswith("S_1=K*F(K,1"), "S1 map drift")

    decisions = artifact["decision"]
    for key in (
        "non_A_theta_current_is_Hiary_j0_j1_instance",
        "physical_quadratic_parameter_already_in_contracting_range",
        "every_main_step_contracts_length_by_at_least_one_half",
        "completed_saddle_current_is_pure_absolute_dual_first_moment",
        "large_constant_phase_cancels_in_completed_saddle_current",
        "zeroth_and_first_reindexed_dual_moments_must_remain_tied",
    ):
        require(decisions.get(key) is True, f"missing exact self-duality decision: {key}")
    for key in (
        "endpoint_complete_R0_R1_evaluator_built",
        "small_b_Euler_Maclaurin_branch_built",
        "complete_fast_theta_evaluator_built",
        "uniform_error_theorem_proved",
        "physical_quadrature_completed",
        "non_A_bound_proved",
        "R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    rows = artifact["normalization_identity_rows"]
    require(len(rows) == 4 and sum(row["conjugated"] for row in rows) == 2, "normalization branch roster drift")
    require({row["branch"] for row in rows} == {"unit_period", "half_shift_then_conjugate", "half_shift", "unit_shift_then_conjugate"}, "normalization branches incomplete")

    traces = artifact["recursion_traces"]
    require(len(traces) == 8, "recursion trace count drift")
    for trace in traces:
        for row in trace["rows"]:
            a = Fraction(row["a"])
            b = Fraction(row["b"])
            require(0 <= a < 1 and 0 <= b <= Fraction(1, 4), "trace normalization range drift")
            require(row["p"] == ceil_fraction(a), "trace p drift")
            require(row["q"] == floor_fraction(a + 2 * b * row["K"]), "trace q drift")
            if "contracted_K" in row:
                require(row["contracted_K"] == row["q"], "contracted length drift")
                require(row["q"] <= (row["K"] + 1) // 2, "trace did not halve")
    require(artifact["trace_summary"]["maximum_main_steps"] <= 22, "too many exact recursion steps")

    pilot = artifact["floating_saddle_only_rows"]
    require(len(pilot) == 6 and all(row["diagnostic_only"] is True for row in pilot), "pilot classification drift")
    for row in pilot:
        x_value = Fraction(row["x"])
        s_value = Fraction(row["s"])
        c_value = A + 2 * s_value
        p = ceil_fraction(x_value * c_value / 2)
        q = floor_fraction(x_value * c_value / 2 + x_value * K)
        require(row["stationary_integer_range"] == [p, q], "pilot stationary range drift")
        require(row["stationary_term_count"] == max(0, q - p + 1), "pilot stationary count drift")

    source = artifact["primary_source"]
    require(source["arxiv"] == "0711.5002v4" and source["imported_numerical_constants"] is False, "primary-source boundary drift")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("(PS1)", "(PS3)", "(PS4)", "(PS5)", "(PS6)", "c+2n_v=2v/x"):
        require(token in note, f"note token missing: {token}")
    normalized_note = " ".join(note.split()).lower()
    require("no endpoint/nonstationary remainder evaluator" in normalized_note and "rh" in normalized_note, "proof boundary missing")
    print("independently checked quadratic-theta Poisson current self-duality and contraction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
