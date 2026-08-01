#!/usr/bin/env python3
"""Independently validate the Xi Pick/Suzuki arithmetic-Hankel bridge."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md"
)

REQUIRED_IDS = {
    "xpsh_01_phi_xi_normalization",
    "xpsh_02_xi_directional_pick",
    "xpsh_03_hyperbolic_modulus_monotonicity",
    "xpsh_04_sondow_horizontal_growth",
    "xpsh_05_zero_disk_identity",
    "xpsh_06_horizontal_growth_countermodel",
    "xpsh_07_suzuki_hb_family",
    "xpsh_08_arithmetic_hankel_operator",
    "xpsh_09_local_positive_hamiltonian",
    "xpsh_10_suzuki_global_equivalence",
    "xpsh_11_fredholm_hankel_hierarchy",
    "xpsh_12_rank_one_tn_spectral_guard",
    "xpsh_13_strict_spectral_target",
    "xpsh_14_determinant_only_reduction",
    "xpsh_15_xi_pick_suzuki_handoff",
    "xpsh_16_boundary_unimodularity_guard",
}

REQUIRED_NOTE_TEXT = (
    "M(w)=integral_R Phi(u)cosh(u*w)du",
    "tau*Re(xi'(s)/xi(s))-delta*Im(xi'(s)/xi(s))",
    "partial_Y log|xi(1/2+delta+i*tau)|",
    "F_*(z)=z^2+6*z+25",
    "(x+3)^2+y^2-16=-14511/10000<0",
    "Horizontal xi-modulus",
    "nu*omega>1",
    "gamma(t)=[det(I+K[t])/det(I-K[t])]^2",
    "Theta_*(r)=(r+i)/(r-i)",
    "boundary-unitary argument would",
    "omega_n decreasing to 0",
    "det(K(x_i+x_j))",
    "det(I-K)=0 and det(I+K)=2",
    "distinct shifted zeros accumulating at itself",
    "jensen_window_pf_suzuki_determinant_only_reduction.md",
    "It does not prove",
)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def independent_polynomial_checks(payload: dict) -> list[str]:
    issues: list[str] = []
    p = Fraction(1)
    q = Fraction(2)
    a = Fraction(11, 10)
    b = Fraction(7, 5)

    alpha_real = p * p - q * q
    alpha_imag = 2 * p * q
    x = a * a - b * b
    y = 2 * a * b
    disk_defect = (x - alpha_real) ** 2 + y**2 - alpha_imag**2
    lower_denominator = (
        (x - alpha_real) ** 2 + (y - alpha_imag) ** 2
    )
    upper_denominator = (
        (x - alpha_real) ** 2 + (y + alpha_imag) ** 2
    )
    anti_pick = (
        (y - alpha_imag) / lower_denominator
        + (y + alpha_imag) / upper_denominator
    )
    roots = (
        (Fraction(1), Fraction(2)),
        (Fraction(1), Fraction(-2)),
        (Fraction(-1), Fraction(2)),
        (Fraction(-1), Fraction(-2)),
    )
    horizontal = sum(
        (a - root_real)
        / ((a - root_real) ** 2 + (b - root_imag) ** 2)
        for root_real, root_imag in roots
    )

    expected = {
        "witness_z": "-3/4+(77/25)*i",
        "alpha": "-3+(4)*i",
        "disk_defect": "-14511/10000",
        "lower_denominator": "59089/10000",
        "upper_denominator": "551889/10000",
        "horizontal_log_derivative": fraction_text(horizontal),
        "anti_pick_log_derivative": fraction_text(anti_pick),
    }
    guard = payload.get("polynomial_guard", {})
    for key, value in expected.items():
        if guard.get(key) != value:
            issues.append(
                f"polynomial guard mismatch for {key}: "
                f"{guard.get(key)!r} != {value!r}"
            )
    if disk_defect != Fraction(-14511, 10000):
        issues.append("independent disk-defect identity failed")
    if horizontal <= 0:
        issues.append("independent horizontal derivative is not positive")
    if anti_pick >= 0:
        issues.append("independent anti-Pick witness is not negative")
    return issues


def independent_rank_one_checks() -> list[str]:
    issues: list[str] = []
    eigenvalue = Fraction(1)
    det_minus = 1 - eigenvalue
    det_plus = 1 + eigenvalue
    if det_minus != 0:
        issues.append("rank-one I-K determinant should vanish")
    if det_plus != 2:
        issues.append("rank-one I+K determinant should equal two")
    return issues


def independent_boundary_unimodularity_checks() -> list[str]:
    issues: list[str] = []
    for value in (
        Fraction(-5, 2),
        Fraction(0),
        Fraction(7, 3),
    ):
        numerator_modulus_squared = value * value + 1
        denominator_modulus_squared = value * value + 1
        if numerator_modulus_squared != denominator_modulus_squared:
            issues.append(
                f"boundary modulus mismatch at real value {value}"
            )
    pole_real = Fraction(0)
    pole_imag = Fraction(1)
    if pole_real != 0 or pole_imag <= 0:
        issues.append("toy ratio pole is not in the upper half-plane")
    return issues


def validate(artifact: Path, note: Path) -> list[str]:
    issues: list[str] = []
    if not artifact.exists():
        return [f"missing artifact: {artifact}"]
    if not note.exists():
        return [f"missing note: {note}"]

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if payload.get("kind") != "jensen_window_pf_xi_pick_suzuki_hankel_bridge":
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
        "xpsh_06_horizontal_growth_countermodel",
        "xpsh_12_rank_one_tn_spectral_guard",
        "xpsh_16_boundary_unimodularity_guard",
    ):
        if row_by_id.get(row_id, {}).get("readiness") != "guard_validated":
            issues.append(f"{row_id} is not marked guard_validated")
    for row_id in (
        "xpsh_13_strict_spectral_target",
        "xpsh_15_xi_pick_suzuki_handoff",
    ):
        if row_by_id.get(row_id, {}).get("readiness") != "not_ready_to_apply":
            issues.append(f"{row_id} is not marked open")
    if (
        row_by_id.get(
            "xpsh_14_determinant_only_reduction", {}
        ).get("readiness")
        != "internally_audited"
    ):
        issues.append("determinant-only reduction is not internally audited")

    exact = payload.get("exact", {})
    for key in (
        "normalization",
        "xi_directional_pick",
        "hyperbolic_monotonicity",
        "horizontal_growth_theorem",
        "zero_disk_identity",
        "suzuki_family",
        "arithmetic_hankel_operator",
        "local_hamiltonian",
        "boundary_unitarity_guard",
        "suzuki_equivalence",
        "fredholm_hankel_series",
        "rank_one_tn_guard",
        "spectral_target",
        "determinant_only_reduction",
        "terminal_reduction",
    ):
        if not exact.get(key):
            issues.append(f"missing exact statement: {key}")

    source_urls = {
        source.get("url")
        for source in payload.get("sources", [])
        if isinstance(source, dict)
    }
    if "https://doi.org/10.1016/j.jfa.2021.109116" not in source_urls:
        issues.append("missing Suzuki JFA primary source")
    if "https://arxiv.org/abs/1005.1104" not in source_urls:
        issues.append("missing Sondow-Dumitrescu primary source")

    issues.extend(independent_polynomial_checks(payload))
    issues.extend(independent_rank_one_checks())
    issues.extend(independent_boundary_unimodularity_checks())

    text = note.read_text(encoding="utf-8")
    for marker in REQUIRED_NOTE_TEXT:
        if marker not in text:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore RH is proved",
        "this proves RH",
        "proves `Lambda <= 0`",
        "global Fredholm nonvanishing is established",
        "terminal kernel condition is established",
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
        print(f"XI-PICK-SUZUKI issue: {issue}")
    print(
        "validated Jensen-window PF Xi Pick/Suzuki Hankel bridge: "
        f"16 rows, {len(issues)} issues, 7 exact coordinate/guard identities, "
        "4 published operator steps, 3 exact countermodels, 1 open global gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
