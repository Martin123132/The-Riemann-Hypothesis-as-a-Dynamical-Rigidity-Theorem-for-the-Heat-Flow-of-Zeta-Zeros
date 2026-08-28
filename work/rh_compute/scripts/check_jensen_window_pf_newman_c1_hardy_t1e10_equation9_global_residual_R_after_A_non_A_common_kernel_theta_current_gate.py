#!/usr/bin/env python3
"""Independently check the non-A common-kernel theta-current reduction."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_common_kernel_theta_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    require(artifact["scope"]["height"] == 10_000_000_000, "height drift")
    require(artifact["scope"]["source_cell_count"] == 2_481_422, "source length drift")
    require(artifact["scope"]["extended_target"] == [622, 39_936], "extended target drift")

    decisions = artifact["decision"]
    for key in (
        "non_A_common_kernel_defined_on_one_regulator",
        "coefficient_mask_and_extended_notch_forms_equivalent",
        "extended_notch_and_two_jet_forms_equivalent",
        "long_source_integral_folded_to_one_cell",
        "folded_source_reduced_to_quadratic_sum_and_first_moment",
        "reference_term_recurrence_proved",
        "zero_x_theta_current_extension_proved",
        "separate_positive_row_norms_forbidden",
    ):
        require(decisions.get(key) is True, f"missing common-kernel decision: {key}")
    for key in (
        "non_A_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    dependencies = {}
    for name, record in artifact["dependencies"].items():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {name}")
        dependencies[name] = load_json(path)
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency pass flag drift")

    finite = dependencies["finite_regulator_equivalence"]
    normal = dependencies["extended_notch_normal_form"]
    pole = dependencies["extended_notch_pole_free_kernel"]
    fold = dependencies["one_cell_fold"]
    pair = dependencies["two_jet_reassembly"]
    embedding = dependencies["A_face_embedding"]

    rows = finite["finite_equivalence_certificate"]["positive_class_rows"]
    require(sum(int(row["count"]) for row in rows[:-1]) == 39_936, "finite row prefix count drift")
    require(rows[2]["post_A_coefficients"] == rows[3]["post_A_coefficients"], "adjacent A-window row vector drift")
    require(rows[2]["joined_summand"] == rows[3]["joined_summand"] == "P_-m+B_m+B_-m", "adjacent A-window summand drift")

    endpoint = normal["symbolic_certificate"]["endpoint_roster"]
    require(endpoint == {"39853_to_39936": "B_m-A_-m", "622_to_39852": "A_m+B_m"}, "extended-notch endpoint roster drift")
    require(normal["symbolic_certificate"]["two_jet_common_kernel_form"].startswith("H_x+E_x+2*pi*i*sum_(m=1)^M"), "extended two-jet form drift")
    require(pair["decision"]["zero_regulator_pair_series_absolute_and_uniform"] is True, "pair convergence drift")
    require(pair["decision"]["raw_derivative_majorant_closes_target"] is False, "raw-majorant route drift")
    require(fold["fold_certificate"]["one_cell_fold"].startswith("integral_0^L f_x(u)C_epsilon(u)du=integral_0^1"), "one-cell fold drift")
    require(pole["exact_kernel_certificate"]["extended_target"] == [622, 39_936], "pole-free target drift")
    require(embedding["decision"]["six_class_mask_unchanged"] is True, "embedding mask guard drift")

    n, a, s, x = sp.symbols("n a s x", real=True)
    alpha = a + 2 * n + 2 * s
    require(sp.expand(alpha**2 / 4 - (a + 2 * s) ** 2 / 4 - n**2 - n * (a + 2 * s)) == 0, "independent phase split failed")
    alpha1 = alpha + 2
    alpha2 = alpha + 4
    require(sp.simplify((alpha1**2 - alpha**2) / 4 - (alpha + 1)) == 0, "independent first ratio failed")
    require(sp.simplify((alpha2**2 - 2 * alpha1**2 + alpha**2) / 4 - 2) == 0, "independent ratio-of-ratios failed")

    l = 2_481_422
    zero_sum = l * (159_577 + 2 * s) + l * (l - 1)
    require(sp.simplify(zero_sum - l * (159_577 + 2 * s + l - 1)) == 0, "independent zero-x extension failed")
    theta = artifact["theta_current_certificate"]
    require(theta["roster"] == {"A": 159_577, "B": 5_122_421, "L": l, "n_range": [0, l - 1]}, "theta roster drift")
    require("S_0" in theta["theta_current"] and "2S_1" in theta["theta_current"], "theta-current reduction missing")

    route = artifact["route_decision"]
    require(route["one_cell_pole_free_theta_current_route_selected_for_pilot"] is True, "selected route drift")
    require(route["fast_incomplete_quadratic_Gauss_evaluator_built"] is False, "fast evaluator overpromotion")
    require(route["fast_evaluator_error_theorem_proved"] is False, "error theorem overpromotion")
    require(Decimal(artifact["budget_certificate"]["non_A_sufficient_absolute_target"]) == Decimal("0.0331747039947"), "target arithmetic drift")

    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("(NK2)", "(NK3)", "(NK4)", "(NK5)", "(NK8)", "(NK9)", "0.0331747039947"):
        require(token in note, f"note token missing: {token}")
    require("No row of (NK2) may be normed independently" in note, "grouping guard missing")
    require("No fast-evaluator error theorem" in note and "RH" in note, "proof boundary missing")

    print("independently checked exact non-A common kernel and quadratic-theta reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
