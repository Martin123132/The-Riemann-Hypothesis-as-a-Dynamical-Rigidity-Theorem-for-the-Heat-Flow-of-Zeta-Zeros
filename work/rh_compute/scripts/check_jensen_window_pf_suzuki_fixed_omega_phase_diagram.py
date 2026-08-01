#!/usr/bin/env python3
"""Independently validate the fixed-omega Suzuki phase diagram."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_fixed_omega_phase_diagram.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md"
)

REQUIRED_IDS = {
    f"sfp_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "determinant_property"),
        (2, "all_history_to_multiplier"),
        (3, "boundary_symmetry"),
        (4, "power_independence"),
        (5, "pole_coordinate"),
        (6, "multiplicity_cancellation"),
        (7, "no_poles_to_inner"),
        (8, "inner_to_strict_compressions"),
        (9, "fixed_pair_equivalence"),
        (10, "cancellation_set"),
        (11, "generic_shift"),
        (12, "zero_width"),
        (13, "phase_diagram"),
        (14, "local_failure_if_not_rh"),
        (15, "closure_criterion"),
        (16, "chain_countermodel"),
        (17, "multiplicity_countermodel"),
        (18, "growth_guard"),
        (19, "summatory_l2"),
        (20, "summatory_sign"),
        (21, "open_scalar_target"),
    )
}

REQUIRED_EXACT = {
    "determinant_property",
    "fixed_pair_equivalence",
    "power_independence",
    "boundary_unimodularity",
    "pole_coordinate",
    "multiplicity_rule",
    "no_poles_to_inner",
    "inner_to_determinants",
    "determinants_to_inner",
    "cancellation_set",
    "zero_width",
    "phase_diagram",
    "local_exclusion",
    "closure_criterion",
    "summatory_coefficients",
    "summatory_function",
    "summatory_l2",
    "summatory_sign",
    "open_scalar_target",
    "growth_guard",
}

REQUIRED_NOTE_TEXT = (
    "D(omega)",
    "Theta_omega is meromorphic inner in C+",
    "m_xi(rho-2*omega)>=m_xi(rho)",
    "Phragmen-Lindelof",
    "C_xi is countable",
    "delta_xi",
    "0 belongs to closure(D)",
    "F(z-i)/F(z+i)=(z-4i)/(z+4i)",
    "pole at z=2i",
    "exp(-i*z)",
    "c_omega(n)",
    "h_omega^<1>(x)",
    "belongs to L2(1,infinity)",
    "one sign for all sufficiently large x",
    "finite samples do not certify it",
    "independent expert review",
    "not a proof",
)


def independent_coordinate_check() -> list[str]:
    issues: list[str] = []
    beta = sp.Rational(3, 4)
    gamma = sp.Integer(7)
    omega = sp.Rational(1, 8)
    z = -gamma + sp.I * (beta - sp.Rational(1, 2) - omega)
    denominator = sp.Rational(1, 2) + omega - sp.I * z
    numerator = sp.Rational(1, 2) - omega - sp.I * z
    rho = beta + sp.I * gamma
    if sp.simplify(denominator - rho) != 0:
        issues.append("denominator coordinate does not recover rho")
    if sp.simplify(numerator - (rho - 2 * omega)) != 0:
        issues.append("numerator coordinate does not recover rho-2*omega")
    if sp.im(z) <= 0:
        issues.append("test pole coordinate is not in C+")
    return issues


def independent_polynomial_check() -> list[str]:
    issues: list[str] = []
    z = sp.symbols("z")
    i = sp.I
    f = lambda w: (w**2 + 1) * (w**2 + 9)
    ratio = sp.cancel(f(z - i) / f(z + i))
    if sp.simplify(ratio - (z - 4 * i) / (z + 4 * i)) != 0:
        issues.append("chain countermodel factorization failed")

    f_def = lambda w: (w**2 + 1) * (w**2 + 9) ** 2
    ratio_def = sp.cancel(f_def(z - i) / f_def(z + i))
    expected = (
        (z - 4 * i) ** 2
        * (z + 2 * i)
        / ((z - 2 * i) * (z + 4 * i) ** 2)
    )
    if sp.simplify(ratio_def - expected) != 0:
        issues.append("multiplicity countermodel factorization failed")
    denominator = sp.denom(sp.cancel(ratio_def))
    if sp.simplify(denominator.subs(z, 2 * i)) != 0:
        issues.append("multiplicity defect does not leave a pole at 2i")
    return issues


def independent_local_exclusion_check() -> list[str]:
    issues: list[str] = []
    displacement = sp.Rational(1, 4)
    isolation_radius = sp.Rational(1, 10)
    epsilon = min(displacement, isolation_radius / 2)
    probes = [epsilon / 4, epsilon / 2, 3 * epsilon / 4]
    if not all(0 < 2 * omega < isolation_radius for omega in probes):
        issues.append("local-exclusion probes leave the isolating disk")
    if not all(omega < displacement for omega in probes):
        issues.append("local-exclusion probes do not map rho into C+")
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
        != "jensen_window_pf_suzuki_fixed_omega_phase_diagram"
    ):
        issues.append("unexpected artifact kind")

    rows = payload.get("rows", [])
    if len(rows) != 21:
        issues.append(f"expected 21 rows, found {len(rows)}")
    ids = {row.get("id") for row in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            f"row id mismatch: missing={sorted(REQUIRED_IDS-ids)} "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )

    row_by_id = {row.get("id"): row for row in rows}
    for row_id in (
        "sfp_09_fixed_pair_equivalence",
        "sfp_11_generic_shift",
        "sfp_13_phase_diagram",
        "sfp_14_local_failure_if_not_rh",
        "sfp_15_closure_criterion",
    ):
        if row_by_id.get(row_id, {}).get("readiness") != "internally_audited":
            issues.append(f"{row_id} is not internally audited")
    for row_id in (
        "sfp_16_chain_countermodel",
        "sfp_17_multiplicity_countermodel",
        "sfp_18_growth_guard",
    ):
        if row_by_id.get(row_id, {}).get("readiness") != "guard_validated":
            issues.append(f"{row_id} is not guard validated")
    if (
        row_by_id.get("sfp_21_open_scalar_target", {}).get("readiness")
        != "not_ready_to_apply"
    ):
        issues.append("scalar arithmetic target is not marked open")

    exact = payload.get("exact", {})
    missing_exact = REQUIRED_EXACT - set(exact)
    if missing_exact:
        issues.append(f"missing exact statements: {sorted(missing_exact)}")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 21,
        "exact_countermodel_count": 3,
        "published_scalar_target_count": 2,
        "open_arithmetic_gate_count": 1,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit mismatch for {key}: {audit.get(key)!r} != {expected!r}"
            )

    urls = {
        source.get("url")
        for source in payload.get("sources", [])
        if isinstance(source, dict)
    }
    for required in (
        "https://arxiv.org/abs/1606.05726",
        "https://arxiv.org/abs/1204.1827",
        "https://arxiv.org/abs/1204.1823",
        "https://link.springer.com/article/10.1007/s00498-024-00387-4",
    ):
        if required not in urls:
            issues.append(f"missing primary source: {required}")

    issues.extend(independent_coordinate_check())
    issues.extend(independent_polynomial_check())
    issues.extend(independent_local_exclusion_check())

    text = note.read_text(encoding="utf-8")
    lower = text.lower()
    for marker in REQUIRED_NOTE_TEXT:
        if marker.lower() not in lower:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "we have proved rh",
        "this proves rh",
        "therefore rh is true",
        "proves lambda<=0",
        "finite positivity establishes eventual sign",
        "published fixed-omega determinant equivalence",
    ):
        if forbidden in lower:
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
        print(f"SUZUKI-FIXED-OMEGA issue: {issue}")
    print(
        "validated Suzuki fixed-omega phase diagram: "
        f"21 rows, {len(issues)} issues, "
        "3 exact countermodel guards, 2 published scalar targets, "
        "1 open arithmetic gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
