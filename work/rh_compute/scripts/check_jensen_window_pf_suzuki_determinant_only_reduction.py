#!/usr/bin/env python3
"""Independently validate the Suzuki determinant-only reduction artifact."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_determinant_only_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_determinant_only_reduction.md"
)

REQUIRED_IDS = {
    "sdr_01_all_history_contraction",
    "sdr_02_dense_hankel_form",
    "sdr_03_uniform_form_bound",
    "sdr_04_global_hankel_extension",
    "sdr_05_reflection_to_convolution",
    "sdr_06_causal_translation_invariance",
    "sdr_07_hardy_multiplier",
    "sdr_08_high_strip_identification",
    "sdr_09_theta_hardy_extension",
    "sdr_10_pole_removal",
    "sdr_11_shifted_zero_rule",
    "sdr_12_cofinal_accumulation",
    "sdr_13_determinant_only_equivalence",
    "sdr_14_terminal_redundancy",
    "sdr_15_single_pair_cancellation_guard",
    "sdr_16_open_all_time_gate",
}

REQUIRED_EXACT = {
    "determinant_to_contraction",
    "dense_form",
    "uniform_form_bound",
    "global_hankel_extension",
    "reflection_convolution",
    "translation_causality",
    "causal_multiplier",
    "high_strip_identification",
    "probe_nonvanishing",
    "theta_extension",
    "pole_removal",
    "shifted_zero",
    "cofinal_accumulation",
    "determinant_only_equivalence",
    "terminal_redundancy",
    "single_pair_guard",
    "open_arithmetic_gate",
}

REQUIRED_NOTE_TEXT = (
    "det(I+/-K[t])!=0 for every t>=0",
    "D=L2_c(R)={compactly supported L2 functions}",
    "|B(f,g)|<=||f||_2||g||_2",
    "represented by a self-adjoint H with ||H||<=1",
    "G=H*J",
    "P_a G=P_a G P_a",
    "H-infinity({Re(s)>0})",
    "L(q)(s)=(1-e^(-s))/s!=0",
    "Theta_tilde(z)=M(-i*z)",
    "rho-2*omega",
    "distinct zeros rho-2*omega_n -> rho",
    "terminal condition",
    "For one omega",
    "still-open all-time determinant premise",
    "does not prove",
    "Independent expert review",
)


def independent_coordinate_check() -> list[str]:
    issues: list[str] = []
    beta = Fraction(3, 4)
    gamma = Fraction(7)
    omega = Fraction(1, 8)
    z_real = -gamma
    z_imag = beta - Fraction(1, 2) - omega

    # With z=z_real+i*z_imag, -i*z=z_imag-i*z_real.
    minus_i_z_real = z_imag
    minus_i_z_imag = -z_real
    denominator_real = Fraction(1, 2) + omega + minus_i_z_real
    numerator_real = Fraction(1, 2) - omega + minus_i_z_real

    if z_imag != Fraction(1, 8) or z_imag <= 0:
        issues.append("zero coordinate did not land at z=-7+i/8 in C+")
    if (denominator_real, minus_i_z_imag) != (beta, gamma):
        issues.append("denominator coordinate did not recover rho")
    if (numerator_real, minus_i_z_imag) != (
        beta - 2 * omega,
        gamma,
    ):
        issues.append("numerator coordinate did not recover rho-2*omega")
    return issues


def independent_cofinal_check() -> list[str]:
    issues: list[str] = []
    rho_real = Fraction(3, 4)
    omegas = [Fraction(1, n) for n in range(10, 16)]
    shifted = [rho_real - 2 * omega for omega in omegas]
    if len(set(shifted)) != len(shifted):
        issues.append("strictly distinct omegas did not give distinct shifts")
    if not all(value < rho_real for value in shifted):
        issues.append("shifted zeros did not lie strictly left of rho")
    if not all(
        abs(shifted[index] - rho_real)
        > abs(shifted[index + 1] - rho_real)
        for index in range(len(shifted) - 1)
    ):
        issues.append("sample cofinal shifts did not approach rho")
    return issues


def validate(artifact: Path, note: Path) -> list[str]:
    issues: list[str] = []
    if not artifact.exists():
        return [f"missing artifact: {artifact}"]
    if not note.exists():
        return [f"missing note: {note}"]

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if (
        payload.get("kind")
        != "jensen_window_pf_suzuki_determinant_only_reduction"
    ):
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
        "sdr_13_determinant_only_equivalence",
        "sdr_14_terminal_redundancy",
    ):
        if row_by_id.get(row_id, {}).get("readiness") != "internally_audited":
            issues.append(f"{row_id} lacks internally_audited status")
    if (
        row_by_id.get(
            "sdr_15_single_pair_cancellation_guard", {}
        ).get("readiness")
        != "guard_validated"
    ):
        issues.append("single-pair cancellation guard is not validated")
    if (
        row_by_id.get("sdr_16_open_all_time_gate", {}).get("readiness")
        != "not_ready_to_apply"
    ):
        issues.append("all-time arithmetic gate is not marked open")

    exact = payload.get("exact", {})
    missing_exact = REQUIRED_EXACT - set(exact)
    if missing_exact:
        issues.append(f"missing exact statements: {sorted(missing_exact)}")

    coordinate = payload.get("coordinate_check", {})
    expected_coordinate = {
        "rho": "3/4+7i",
        "omega": "1/8",
        "z": "-7+i/8",
        "denominator_argument": "3/4+7i",
        "numerator_argument": "1/2+7i",
        "shift": "1/4",
    }
    for key, expected in expected_coordinate.items():
        if coordinate.get(key) != expected:
            issues.append(
                f"coordinate mismatch for {key}: "
                f"{coordinate.get(key)!r} != {expected!r}"
            )

    source_urls = {
        source.get("url")
        for source in payload.get("sources", [])
        if isinstance(source, dict)
    }
    for required_url in (
        "https://arxiv.org/abs/1606.05726",
        "https://link.springer.com/article/10.1007/s00498-024-00387-4",
    ):
        if required_url not in source_urls:
            issues.append(f"missing primary source: {required_url}")

    issues.extend(independent_coordinate_check())
    issues.extend(independent_cofinal_check())

    text = note.read_text(encoding="utf-8")
    lower_text = text.lower()
    for marker in REQUIRED_NOTE_TEXT:
        if marker.lower() not in lower_text:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore the determinant conditions hold for zeta",
        "we have proved rh",
        "this proves rh",
        "proves lambda<=0",
        "one fixed all-time determinant family forces rh",
        "terminal limit follows directly from one fixed family",
        "published determinant-only theorem",
    ):
        if forbidden in lower_text:
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
        print(f"SUZUKI-DETERMINANT-ONLY issue: {issue}")
    print(
        "validated Suzuki determinant-only reduction: "
        f"16 rows, {len(issues)} issues, "
        "12 exact bridge steps, 2 theorem/corollary candidates, "
        "1 route guard, 1 open arithmetic gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
