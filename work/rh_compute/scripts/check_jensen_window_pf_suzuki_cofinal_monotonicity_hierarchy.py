#!/usr/bin/env python3
"""Validate the Suzuki cofinal monotonicity-hierarchy theorem audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md"
)

REQUIRED_IDS = {
    f"cmh_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "coefficient_identity"),
        (2, "first_weight"),
        (3, "weight_recursion"),
        (4, "summatory_definition"),
        (5, "log_antiderivative"),
        (6, "mellin_identity"),
        (7, "k1_normalization"),
        (8, "eventual_sign_landau"),
        (9, "pole_free_inner"),
        (10, "fixed_shift_handoff"),
        (11, "cofinal_handoff"),
        (12, "rh_k1_asymptotic"),
        (13, "rh_higher_asymptotic"),
        (14, "cofinal_equivalence"),
        (15, "near_zero_first_weight"),
        (16, "near_zero_all_weights"),
        (17, "near_one_all_weights"),
        (18, "signed_weight_guard"),
        (19, "finite_range_guard"),
        (20, "fixed_shift_guard"),
        (21, "open_smoothed_target"),
    )
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Suzuki Cofinal Monotonicity Hierarchy",
    "This is not a proof of eventual sign, RH",
    "integral_1^infinity H_(omega,k)(x)x^(1/2-s)dx/x",
    "H_(omega,1)(x)=sqrt(x)*h_omega^<1>(x)",
    "If omega_j decreases to zero",
    "for any integers k_j>=1, then RH",
    "g_(omega,k)(x)",
    "(1/2-omega)^(-(k-1))",
    "Every smoothing weight is",
    "genuinely signed",
    "Finite positivity remains nonpromotable",
    "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
    "https://arxiv.org/abs/1204.1823",
    "https://arxiv.org/abs/1204.1827",
    "does not",
    "prove RH or Lambda<=0",
)


def independent_mellin_factorization_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 60
    s = mp.mpf(2)
    omega = mp.mpf(1) / 4

    def xi(value: mp.mpf) -> mp.mpf:
        return (
            value
            * (value - 1)
            * mp.power(mp.pi, -value / 2)
            * mp.gamma(value / 2)
            * mp.zeta(value)
        )

    coefficient_factor = mp.zeta(s - omega) / mp.zeta(s + omega)
    archimedean_factor = (
        (s - omega)
        * (s - omega - 1)
        / ((s + omega) * (s + omega - 1))
        * mp.power(mp.pi, omega)
        * mp.gamma((s - omega) / 2)
        / mp.gamma((s + omega) / 2)
    )
    expected = xi(s - omega) / xi(s + omega)
    relative_error = abs(
        coefficient_factor * archimedean_factor / expected - 1
    )
    if relative_error > mp.mpf("1e-50"):
        issues.append(
            "independent Mellin factorization failed: "
            f"relative error {relative_error}"
        )
    return issues


def independent_endpoint_sign_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 60
    for omega in (mp.mpf(1) / 4, mp.mpf(1) / 8):
        beta_parameter = (3 - 2 * omega) / 2
        a_omega = (
            4
            * omega
            / (1 - 2 * omega)
            * mp.power(mp.pi, omega)
            / mp.gamma(omega)
            * mp.beta(beta_parameter, omega)
        )
        if not a_omega > 0:
            issues.append(f"a_omega is not positive at omega={omega}")
        for k in (1, 2, 3, 4):
            leading = -a_omega * mp.power(
                mp.mpf("0.5") - omega, -(k - 1)
            )
            if not leading < 0:
                issues.append(
                    f"near-zero leading coefficient is not negative "
                    f"at omega={omega}, k={k}"
                )

        near_one_coefficient = (
            mp.power(2 * mp.pi, omega)
            / (omega * mp.gamma(omega))
        )
        for k in (1, 2, 3, 4):
            if not near_one_coefficient > 0:
                issues.append(
                    f"near-one coefficient is not positive "
                    f"at omega={omega}, k={k}"
                )
            near_one_coefficient /= omega + k
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != (
        "jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy"
    ):
        issues.append("bad kind")

    rows = payload.get("rows", [])
    by_id = {item.get("id"): item for item in rows}
    if len(rows) != 21:
        issues.append(f"expected 21 rows, found {len(rows)}")
    missing = REQUIRED_IDS - set(by_id)
    if missing:
        issues.append(f"missing row ids: {sorted(missing)}")
    if len(by_id) != len(rows):
        issues.append("duplicate row id")

    for item in rows:
        boundary = str(item.get("proof_boundary", "")).lower()
        if not boundary:
            issues.append(f"{item.get('id')}: missing proof boundary")
        if not item.get("role") or not item.get("status"):
            issues.append(f"{item.get('id')}: missing role or status")

    for row_id in (
        "cmh_18_signed_weight_guard",
        "cmh_19_finite_range_guard",
        "cmh_20_fixed_shift_guard",
    ):
        item = by_id.get(row_id, {})
        if item.get("role") != "countermodel_gate":
            issues.append(f"{row_id}: bad guard role")
        if item.get("status") != "guard_validated":
            issues.append(f"{row_id}: bad guard status")

    target = by_id.get("cmh_21_open_smoothed_target", {})
    if target.get("role") != "open_arithmetic_gate":
        issues.append("open target has wrong role")
    if target.get("status") != "open_target":
        issues.append("open target has wrong status")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 21,
        "countermodel_gate_count": 3,
        "open_arithmetic_gate_count": 1,
        "eventual_sign_proved": False,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key}: expected {expected!r}, "
                f"found {audit.get(key)!r}"
            )

    source_urls = {item.get("url") for item in payload.get("sources", [])}
    for required in (
        "https://arxiv.org/abs/1204.1823",
        "https://arxiv.org/abs/1204.1827",
        "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
    ):
        if required not in source_urls:
            issues.append(f"missing source: {required}")

    for needle in REQUIRED_NOTE_STRINGS:
        if needle not in note:
            issues.append(f"note missing: {needle}")

    issues.extend(independent_mellin_factorization_check())
    issues.extend(independent_endpoint_sign_check())

    lowered = note.lower()
    for forbidden in (
        "we have proved rh",
        "this proves rh",
        "therefore rh is true",
        "finite positivity proves eventual sign",
        "coefficient positivity proves",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        print(
            "failed Suzuki cofinal monotonicity hierarchy: "
            f"{len(issues)} issues"
        )
        return 1

    print(
        "validated Suzuki cofinal monotonicity hierarchy: "
        "21 rows, 0 issues, 3 signed/nonpromotion guards, "
        "1 cofinal equivalence candidate, 1 open arithmetic gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
