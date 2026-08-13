#!/usr/bin/env python3
"""Validate the cubic level-one-to-level-two Legendre adapter gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.py"
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


def fraction(record: dict[str, int]) -> Fraction:
    value = Fraction(int(record["numerator"]), int(record["denominator"]))
    require(value.denominator > 0, "invalid rational record")
    return value


def dyadic(value: list[int]) -> Fraction:
    require(isinstance(value, list) and len(value) == 2, "invalid dyadic record")
    mantissa, exponent = int(value[0]), int(value[1])
    return Fraction(mantissa * 2**exponent, 1) if exponent >= 0 else Fraction(mantissa, 2 ** (-exponent))


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER, SOURCE):
        require(path.is_file(), f"missing cubic-child artifact: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate",
        "cubic-child kind drift",
    )
    require(
        artifact["status"] == "exact_cubic_child_adapter_with_finite_bit_and_analytic_q_tail_open",
        "cubic-child status drift",
    )
    require(artifact["source"]["sha256"] == SOURCE_SHA256, "result source hash drift")
    require(len(artifact["source"]["locations"]) == 16, "source location inventory drift")

    lemma = artifact["exact_cubic_lemma"]
    require("u=y*X+3*a3*X^2" in lemma["legendre_dual"], "Legendre equation drift")
    require("H3(u)=u^2/(2*y)-a3*u^3/y^3" in lemma["source_truncation"], "H3 formula drift")
    require("9*a3^2/(2*y^5)" in lemma["first_omitted_term"], "h4 formula drift")
    require("binomial expansion" in lemma["translation"], "translation identity drift")
    require("integer" in lemma["integer_normalization"], "integer normalization drift")
    require("same positive normalized phase" in lemma["orientation"], "orientation conclusion drift")
    require("exp(i*tpp*r0)" in lemma["constant_phase"], "constant phase drift")

    denominator = artifact["kernel_denominator_lemma"]
    require("1+h1*u" in denominator["centered_form"], "centered denominator drift")
    require("denn+k*(rminorcor-k*rminorcor2)" in denominator["source_form"], "source denominator drift")
    require(denominator["saved_checks"] >= 128, "saved denominator coverage drift")
    require(denominator["synthetic_raised_order_checks"] == 168, "synthetic denominator count drift")
    require(fraction(denominator["minimum_saved_absolute_denominator"]) > Fraction(9, 10), "denominator margin too small")

    saved = artifact["saved_fixture"]
    require(saved["chain_count"] == 64, "saved chain count drift")
    require(saved["branch_counts"] == {"1": 32, "2": 32}, "saved branch count drift")
    require(saved["mit_histogram"] == {"2": 64}, "MIT roster drift")
    require(saved["degree_pair_histogram"] == {"3_to_3": 64}, "degree pair drift")
    require(saved["parent_length_range"] == [104, 104], "parent length drift")
    require(saved["child_length_range"] == [1, 2], "child length drift")
    require(saved["integer_phase_checks"] == denominator["saved_checks"], "phase/denominator roster mismatch")
    require(len(saved["rows"]) == 64, "saved row count drift")
    require([row["chain"] for row in saved["rows"]] == list(range(1, 65)), "saved sequence gap")
    require(all(row["parent_degree"] == row["child_degree"] == 3 for row in saved["rows"]), "saved degree row drift")
    require(all(fraction(row["source_raise_indicator"]) < Fraction(1, 1000) for row in saved["rows"]), "raise threshold failure")
    require(fraction(saved["maximum_source_raise_indicator"]) < Fraction(1, 100000), "raise indicator scale drift")
    require(fraction(saved["maximum_exact_raw_coefficient_gap"]) < Fraction(1, 10**28), "coefficient gap too large")
    for margin in saved["minimum_selector_margins"].values():
        require(fraction(margin["exact"]) > 0, "saved selector tie")

    rigorous = artifact["rigorous_point_reconstruction"]
    require(rigorous["precision_ladder_dps"] == [96, 160], "precision ladder drift")
    kernel_gap = dyadic(rigorous["maximum_weighted_kernel_source_roundoff_gap_upper"]["upper_dyadic"])
    tail = dyadic(rigorous["maximum_exact_legendre_tail_at_saved_child_indices"]["upper_dyadic"])
    z_value = dyadic(rigorous["maximum_dimensionless_branch_parameter_abs"]["upper_dyadic"])
    discriminant = dyadic(rigorous["minimum_stationary_discriminant"]["lower_dyadic"])
    require(kernel_gap < Fraction(1, 10**28), "weighted-kernel gap too large")
    require(Fraction(0) < tail < Fraction(1, 10000), "Legendre tail scale drift")
    require(z_value < Fraction(1, 100), "stationary branch parameter too large")
    require(discriminant > Fraction(99, 100), "stationary discriminant margin drift")
    require("not a uniform" in rigorous["scope"], "point-scope guard drift")

    synthetic = artifact["synthetic_normalization"]
    require(synthetic["fixture_count"] >= 60, "synthetic normalization coverage drift")
    require(synthetic["integer_phase_checks"] == 17 * synthetic["fixture_count"], "synthetic parity count drift")
    require(synthetic["ix_parities"] == [0, 1], "synthetic parity branches missing")
    require(synthetic["xr_signs"] == [-1, 1], "synthetic orientation branches missing")
    require(synthetic["conjugation_values"] == [False, True], "synthetic conjugation branches missing")

    formal = artifact["formal_raised_order"]
    require(formal["case_count"] == 20, "formal series case count drift")
    require(formal["coefficient_checks"] == 76, "formal coefficient count drift")
    require(formal["zero_extension_checks"] > 100, "zero-extension coverage drift")
    require(formal["degrees"] == list(range(3, 8)), "formal degree coverage drift")
    require("source formulas above quartic" in formal["conclusion"], "raised-order scope guard drift")

    conclusion = artifact["route_conclusion"]
    require("exact cubic formal Legendre transform" in conclusion["closed"], "closed scope drift")
    require("All 64 MIT=2 chains" in conclusion["saved_source_case"], "saved case drift")
    require("L(1)=1 or 2" in conclusion["new_observation"], "length distinction missing")
    require("heuristic, not a tail theorem" in conclusion["open"], "tail nonpromotion guard missing")
    require("outward-rounded interval cell" in conclusion["next_target"], "next target drift")

    sources = artifact["sources"]
    for key in ("telemetry", "initial_adapter"):
        path = REPO_ROOT / sources[key]["path"]
        require(file_hash(path) == sources[key]["sha256"], f"{key} dependency hash drift")
    require(file_hash(BUILDER) == sources["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == sources["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: exact cubic transformed-child adapter; not a proof and analytic q/tail enclosure open",
        "`H(u)=u^2/(2*y)-a3*u^3/y^3+[9*a3^2/(2*y^5)]u^4+...`",
        "The raise-order test therefore uses `L(1)`, not `L(0)`.",
        "These are finite exact-input point statements, not interval theorems.",
        "does not turn the source's first-omitted-term criterion into a rigorous tail bound",
        "prize-level conclusion",
    ):
        require(token in note, f"cubic-child note token missing: {token}")

    print(
        "validated Hardy cubic child Legendre adapter gate: "
        f"{saved['chain_count']} chains, {saved['integer_phase_checks']} parity checks, "
        f"{denominator['saved_checks']} denominator checks, {formal['case_count']} formal series cases"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
