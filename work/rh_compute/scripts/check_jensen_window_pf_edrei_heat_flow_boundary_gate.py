#!/usr/bin/env python3
"""Independently validate the Edrei heat-flow boundary gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md"

REQUIRED_IDS = {
    "ehfb_01_radial_heat_operator",
    "ehfb_02_log_derivative_pde",
    "ehfb_03_edrei_moment_ode",
    "ehfb_04_stieltjes_generating_pde",
    "ehfb_05_rank_one_boundary",
    "ehfb_06_rank_one_orientation",
    "ehfb_07_double_zero_heat_polynomial",
    "ehfb_08_double_zero_discriminant",
    "ehfb_09_all_shift_minor_crossing",
    "ehfb_10_backward_invariance_guard",
    "ehfb_11_xi_specific_handoff",
}

REQUIRED_NOTE_TEXT = (
    "partial_lambda F_lambda=(4*z*partial_z^2+2*partial_z)F_lambda",
    "partial_lambda R=4*z*R_zz+8*z*R*R_z+6*R_z+4*R^2",
    "a_n'= -2*(n+1)*(2*n+3)*a_(n+1)",
    "partial_lambda Delta_s=8*m^2*(m-1)*beta^(2*s+5)",
    "Disc_z=32*lambda*beta^3*(1+3*lambda*beta)",
    "generic backward Stieltjes-cone invariance is false",
    "Xi/Phi-specific all-order rigidity",
    "This artifact proves",
    "It does not prove",
)


def flow_rhs(values: list[Fraction], n: int) -> Fraction:
    convolution = sum(values[k] * values[n - k] for k in range(n + 1))
    return (
        -2 * (n + 1) * (2 * n + 3) * values[n + 1]
        + 4 * (n + 1) * convolution
    )


def independent_rank_one_checks() -> list[str]:
    issues: list[str] = []
    beta = Fraction(3, 5)
    for multiplicity in (1, 2, 3):
        values = [
            Fraction(multiplicity) * beta ** (n + 1) for n in range(7)
        ]
        derivatives = [flow_rhs(values, n) for n in range(6)]
        for shift in (0, 1, 2):
            observed = (
                derivatives[shift] * values[shift + 2]
                + values[shift] * derivatives[shift + 2]
                - 2 * values[shift + 1] * derivatives[shift + 1]
            )
            expected = (
                8
                * multiplicity**2
                * (multiplicity - 1)
                * beta ** (2 * shift + 5)
            )
            if observed != expected:
                issues.append(
                    f"rank-one orientation mismatch m={multiplicity} s={shift}: "
                    f"{observed} != {expected}"
                )
    return issues


def independent_quadratic_checks() -> list[str]:
    issues: list[str] = []
    beta = Fraction(1)
    for lam, expected_sign in ((Fraction(1, 10), 1), (Fraction(-1, 10), -1)):
        leading = beta**2
        linear = 2 * beta + 12 * lam * beta**2
        constant = 1 + 4 * lam * beta + 12 * lam**2 * beta**2
        discriminant = linear**2 - 4 * leading * constant
        expected = 32 * lam * beta**3 * (1 + 3 * lam * beta)
        if discriminant != expected:
            issues.append(f"quadratic discriminant mismatch at lambda={lam}")
        if constant <= 0:
            issues.append(f"quadratic constant is not positive at lambda={lam}")
        if expected_sign * discriminant <= 0:
            issues.append(f"quadratic discriminant has wrong sign at lambda={lam}")
        for shift in (0, 1, 2):
            observed_minor = (
                beta ** (2 * shift + 2)
                * discriminant
                / constant ** (shift + 3)
            )
            if expected_sign * observed_minor <= 0:
                issues.append(
                    f"shifted minor has wrong sign at lambda={lam}, s={shift}"
                )
    return issues


def validate(artifact: Path, note: Path) -> list[str]:
    issues: list[str] = []
    if not artifact.exists():
        return [f"missing artifact: {artifact}"]
    if not note.exists():
        return [f"missing note: {note}"]

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if payload.get("kind") != "jensen_window_pf_edrei_heat_flow_boundary_gate":
        issues.append("unexpected artifact kind")
    rows = payload.get("rows", [])
    if len(rows) != 11:
        issues.append(f"expected 11 rows, found {len(rows)}")
    ids = {row.get("id") for row in rows}
    if ids != REQUIRED_IDS:
        issues.append(f"row id mismatch: missing={sorted(REQUIRED_IDS-ids)} extra={sorted(ids-REQUIRED_IDS)}")

    row_by_id = {row.get("id"): row for row in rows}
    if row_by_id.get("ehfb_10_backward_invariance_guard", {}).get("readiness") != "rejected_by_countermodel":
        issues.append("generic backward-invariance shortcut is not rejected")
    if row_by_id.get("ehfb_11_xi_specific_handoff", {}).get("readiness") != "not_ready_to_apply":
        issues.append("Xi/Phi handoff is not marked open")

    exact = payload.get("exact", {})
    for key in (
        "radial_heat",
        "log_derivative_pde",
        "moment_ode",
        "generating_pde",
        "rank_one_orientation",
        "quadratic_discriminant",
        "quadratic_shifted_minor",
        "backward_failure",
        "open_handoff",
    ):
        if not exact.get(key):
            issues.append(f"missing exact statement: {key}")

    issues.extend(independent_rank_one_checks())
    issues.extend(independent_quadratic_checks())

    text = note.read_text(encoding="utf-8")
    for marker in REQUIRED_NOTE_TEXT:
        if marker not in text:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore RH",
        "proves RH",
        "proves `Lambda <= 0`",
        "the Xi/Phi invariant is proved",
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
        print(f"EDREI-HEAT-BOUNDARY issue: {issue}")
    print(
        "validated Jensen-window PF Edrei heat-flow boundary gate: "
        f"11 rows, {len(issues)} issues, 4 exact flow identities, "
        "9 rank-one orientation checks, 2 exact heat witnesses, "
        "1 rejected generic backward-invariance shortcut, 1 open Xi/Phi handoff"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
