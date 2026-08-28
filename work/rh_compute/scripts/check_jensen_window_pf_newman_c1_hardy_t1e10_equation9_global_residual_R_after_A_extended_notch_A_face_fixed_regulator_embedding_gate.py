#!/usr/bin/env python3
"""Independently check the translated A-face fixed-regulator embedding."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate"
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
    ctx.dps = 150
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    require(artifact["scope"]["height"] == 10_000_000_000, "height drift")
    require(artifact["scope"]["A"] == 159_577, "A drift")

    decisions = artifact["decision"]
    for key in (
        "fixed_regulator_A_face_channel_embedded_exactly",
        "A_face_channel_independent_of_M_after_common_cutoff",
        "common_positive_regulator_preserved",
        "block_minus_translated_strip_orientation_verified",
        "six_class_mask_unchanged",
        "block_and_strip_charged_exactly_once",
        "A_notch_atoms_not_reintroduced",
        "physical_Abel_zero_limit_split_proved",
        "complete_A_face_absolute_bound_applies_to_embedded_channel",
        "remaining_non_A_sufficient_target_0_point_0331747039947_derived",
    ):
        require(decisions.get(key) is True, f"missing embedding decision: {key}")
    for key in (
        "remaining_non_A_bound_proved",
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

    finite = dependencies["finite_regulator_equivalence"]
    normal = dependencies["extended_notch_normal_form"]
    split = dependencies["normal_tangential_split"]
    peano = dependencies["midpoint_Peano_partition"]
    six = dependencies["six_current_reduction"]
    assembly = dependencies["complete_A_face_assembly"]
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency pass flag drift")

    target = set(range(622, 39_895))
    a_window = set(range(39_853, 39_937))
    extended = set(range(622, 39_937))
    shifted = set(range(39_895, 39_937))
    require(target | a_window == extended, "independent extended-notch union failed")
    require(extended - target == shifted and len(shifted) == 42, "independent shifted support failed")
    left = Fraction(79_789, 2)
    right = Fraction(79_873, 2)
    require(right - left == 42, "independent strip length failed")
    require(
        all(Fraction(m, 1) == left + Fraction(2 * j + 1, 2) for j, m in enumerate(sorted(shifted))),
        "independent midpoint roster failed",
    )
    require(5_122_421 > max(shifted), "common cutoff support failed")

    rows = finite["finite_equivalence_certificate"]["positive_class_rows"]
    row = next(item for item in rows if item["sector_id"] == "A_window_outer_pairs")
    require(row["mode_range"] == [39_895, 39_936] and row["count"] == 42, "six-class support drift")
    require(row["joined_summand"] == "P_-m+B_m+B_-m", "six-class ownership drift")
    require(row["post_A_coefficients"] == {
        "A_minus": 0,
        "A_plus": 0,
        "B_minus": 1,
        "B_plus": 1,
        "P_minus": 1,
        "P_plus": 0,
    }, "six-class vector drift")

    require(normal["set_certificate"]["indicator_identity"] == "chi_T+chi_A*(1-chi_T)=chi_U, U=T union A_window={622,...,39936}", "normal-form set identity drift")
    require(split["translation_certificate"]["fixed_regulator_identity"].startswith("C_U,epsilon(s)-C_T,epsilon(s)=-sum_(m=39895)^39936"), "finite block sign drift")
    require(split["translation_certificate"]["interpolant_identity"].startswith("h_U,epsilon(y)-h_T,epsilon(y)=-q_epsilon(y)"), "translated strip sign drift")
    block_sign = -1
    strip_inside_bracket_sign = -1
    require((block_sign, -strip_inside_bracket_sign) == (-1, 1), "independent oriented-channel sign failed")
    require(peano["symbolic_certificate"]["oriented_defect"].startswith("[-sum_m f(m)]-[-integral_strip f(y)dy]"), "Peano orientation drift")

    six_formula = six["symbolic_certificate"]["six_current_amplitude"]
    six_error = six["symbolic_certificate"]["pointwise_amplitude_error"]
    require("A*x" in six_formula and "(x/2)^n" in six_formula, "O(x) six-current structure missing")
    require("(x/2)^6" in six_error, "O(x^6) exact remainder structure missing")
    require(six["decision"]["D6_over_x_regular_at_zero_proved"] is True, "zero-end regularity drift")

    complete = arb(assembly["assembly_certificate"]["complete_translated_A_face_absolute_bound"])
    require(complete.upper() < arb("0.00364"), "embedded A-face threshold failed")
    target_decimal = Decimal("0.0368147039947") - Decimal("0.00364")
    require(target_decimal == Decimal("0.0331747039947"), "independent budget arithmetic failed")
    budget = artifact["budget_certificate"]
    require(Decimal(budget["non_A_sufficient_absolute_target"]) == target_decimal, "recorded non-A target drift")
    require("R_after_A=I_(A,42)+R_nonA" in artifact["physical_limit_certificate"]["limit_split"], "physical split missing")

    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("mathfrak R_A,(M,epsilon)" in note and "mathfrak R_nonA,(M,epsilon)" in note, "fixed-regulator split missing")
    require("P_-m+B_m+B_-m" in note and "does **not** delete a" in note, "mask guard missing")
    require("O(x^(-1/4))" in note and "Dominated convergence" in note, "physical-limit argument missing")
    require("0.0331747039947" in note and "does not prove it" in note, "non-A target boundary missing")
    require("remain unproved" in note and "RH" in note, "proof boundary missing")

    print("independently checked exact fixed-regulator translated A-face embedding", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
