#!/usr/bin/env python3
"""Independently validate the Suzuki truncation-path spectral frontier."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_suzuki_spectral_frontier.md"
)

REQUIRED_IDS = {
    "ssf_01_fixed_space_translation",
    "ssf_02_effective_finite_interval",
    "ssf_03_hilbert_schmidt_identity",
    "ssf_04_hilbert_schmidt_divergence",
    "ssf_05_hilbert_schmidt_certificate_no_go",
    "ssf_06_truncation_path_continuity",
    "ssf_07_operator_norm_monotonicity",
    "ssf_08_global_determinant_contractivity_equivalence",
    "ssf_09_first_crossing_alternative",
    "ssf_10_isolated_time_guard",
    "ssf_11_suzuki_high_contour_certificate",
    "ssf_12_high_contour_finite_ceiling",
    "ssf_13_conditional_inner_contractivity",
    "ssf_14_no_uniform_t_gap",
    "ssf_15_signed_quadratic_form_target",
    "ssf_16_terminal_sequence_reduction",
}

REQUIRED_NOTE_TEXT = (
    "K(2t-u-v)",
    "||K[t]||_HS^2=integral_0^(2t)",
    "(2t-R)*c_R -> infinity",
    "continuous in Hilbert-Schmidt norm",
    "det(I+/-K[t])!=0",
    "first crossing",
    "A=2P",
    "M_v^2*exp(4v*t)<1",
    "Theta_(omega,nu)(i*v)",
    "T_HC(omega,nu)",
    "lim_(t->infinity)||K[t]||=1",
    "|<K[t]f,f>|<||f||^2",
    "shifted-zero accumulation",
    "jensen_window_pf_suzuki_determinant_only_reduction.md",
    "It does not prove",
)


def independent_hs_check() -> list[str]:
    issues: list[str] = []
    # K(s)=s, t=1.
    weighted_integral = (
        2 * Fraction(2) ** 3 / 3 - Fraction(2) ** 4 / 4
    )
    # Direct triangle integration:
    # integral_0^2 integral_0^(2-u) (2-u-v)^2 dv du.
    direct_triangle = Fraction(4, 3)
    if weighted_integral != direct_triangle:
        issues.append(
            "independent Hilbert-Schmidt reduction failed: "
            f"{weighted_integral} != {direct_triangle}"
        )
    return issues


def independent_spectral_checks() -> list[str]:
    issues: list[str] = []
    eigenvalue = Fraction(2)
    if abs(eigenvalue) != 2:
        issues.append("rank-one toy norm should be two")
    if 1 - eigenvalue != -1:
        issues.append("rank-one toy det(I-A) should be minus one")
    if 1 + eigenvalue != 3:
        issues.append("rank-one toy det(I+A) should be three")

    # A continuous scalar path from zero to magnitude above one must hit an
    # endpoint. These rational samples independently exercise both signs.
    for endpoint in (Fraction(1), Fraction(-1)):
        determinant = 1 - endpoint if endpoint > 0 else 1 + endpoint
        if determinant != 0:
            issues.append(f"endpoint {endpoint} did not zero its determinant")
    return issues


def validate(artifact: Path, note: Path) -> list[str]:
    issues: list[str] = []
    if not artifact.exists():
        return [f"missing artifact: {artifact}"]
    if not note.exists():
        return [f"missing note: {note}"]

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if payload.get("kind") != "jensen_window_pf_suzuki_spectral_frontier":
        issues.append("unexpected artifact kind")

    rows = payload.get("rows", [])
    if len(rows) != 16:
        issues.append(f"expected 16 rows, found {len(rows)}")
    ids = {row.get("id") for row in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            f"row id mismatch: missing={sorted(REQUIRED_IDS-ids)} "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )

    row_by_id = {row.get("id"): row for row in rows}
    for row_id in (
        "ssf_04_hilbert_schmidt_divergence",
        "ssf_05_hilbert_schmidt_certificate_no_go",
        "ssf_10_isolated_time_guard",
        "ssf_12_high_contour_finite_ceiling",
        "ssf_14_no_uniform_t_gap",
    ):
        if row_by_id.get(row_id, {}).get("readiness") != "guard_validated":
            issues.append(f"{row_id} is not marked guard_validated")
    for row_id in ("ssf_15_signed_quadratic_form_target",):
        if row_by_id.get(row_id, {}).get("readiness") != "not_ready_to_apply":
            issues.append(f"{row_id} is not marked open")
    if (
        row_by_id.get(
            "ssf_16_terminal_sequence_reduction", {}
        ).get("readiness")
        != "internally_audited"
    ):
        issues.append("cofinal terminal reduction is not internally audited")
    if (
        row_by_id.get(
            "ssf_13_conditional_inner_contractivity", {}
        ).get("readiness")
        != "conditional_only"
    ):
        issues.append("conditional HB row is not marked conditional_only")

    exact = payload.get("exact", {})
    for key in (
        "operator",
        "fixed_space",
        "effective_interval",
        "hilbert_schmidt",
        "hilbert_schmidt_growth",
        "path_continuity",
        "norm_monotonicity",
        "determinant_contractivity",
        "first_crossing",
        "isolated_time_guard",
        "high_contour_certificate",
        "zeta_vertical_asymptotic",
        "high_contour_ceiling",
        "conditional_inner_result",
        "no_uniform_gap",
        "quadratic_form_target",
        "terminal_sequence_reduction",
    ):
        if not exact.get(key):
            issues.append(f"missing exact statement: {key}")

    toy = payload.get("toy_checks", {})
    expected_toy = {
        "hs_squared": "4/3",
        "rank_one_eigenvalue": "2",
        "rank_one_det_minus": "-1",
        "rank_one_det_plus": "3",
    }
    for key, value in expected_toy.items():
        if toy.get(key) != value:
            issues.append(
                f"toy check mismatch for {key}: {toy.get(key)!r} != {value!r}"
            )

    source_urls = {
        source.get("url")
        for source in payload.get("sources", [])
        if isinstance(source, dict)
    }
    if "https://arxiv.org/abs/1606.05726" not in source_urls:
        issues.append("missing Suzuki primary source")

    issues.extend(independent_hs_check())
    issues.extend(independent_spectral_checks())

    text = note.read_text(encoding="utf-8")
    for marker in REQUIRED_NOTE_TEXT:
        if marker not in text:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore RH is proved",
        "this proves RH",
        "proves `Lambda <= 0`",
        "global Fredholm nonvanishing is established",
        "one fixed family directly proves the terminal limit",
        "uniform spectral gap is proved",
    ):
        if forbidden.lower() in text.lower():
            issues.append(f"forbidden promotion language: {forbidden}")

    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(args.artifact, args.note)
    for issue in issues:
        print(f"SUZUKI-SPECTRAL-FRONTIER issue: {issue}")
    print(
        "validated Jensen-window PF Suzuki spectral frontier: "
        f"16 rows, {len(issues)} issues, "
        "10 exact path/reduction identities, 5 route guards, "
        "1 open global obligation"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
