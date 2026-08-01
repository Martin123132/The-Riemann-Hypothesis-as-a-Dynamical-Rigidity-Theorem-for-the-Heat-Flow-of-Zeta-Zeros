#!/usr/bin/env python3
"""Independently validate the Phi Pick-kernel endpoint target."""

from __future__ import annotations

import argparse
import cmath
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_phi_pick_kernel_target.md"

REQUIRED_IDS = {
    "ppkt_01_phi_transform",
    "ppkt_02_kernel_atoms",
    "ppkt_03_log_derivative_ratio",
    "ppkt_04_pick_numerator",
    "ppkt_05_polarized_double_kernel",
    "ppkt_06_pick_endpoint_equivalence",
    "ppkt_07_first_quadrant_form",
    "ppkt_08_tilted_variance_wall",
    "ppkt_09_origin_kurtosis_wall",
    "ppkt_10_two_scale_pick_failure",
    "ppkt_11_positive_resolvent_route",
    "ppkt_12_xi_phi_pick_handoff",
}

REQUIRED_NOTE_TEXT = (
    "P_Phi(z):=-Im(F'(z)*conj(F(z)))",
    "P_Phi(z)=double_integral",
    "F in LP+ <=> P_Phi(z)>0 for every Im(z)>0",
    "2*|w|^2*P_Phi(w^2)",
    "R'(x)=E_Pi[r_x']+Var_Pi(r_x)",
    "R'(0)=(E[y^2]-3*E[y]^2)/12",
    "R'(0)=7/19200>0",
    "R(z)=gamma+<v,(I+z*T)^(-1)v>",
    "endpoint-equivalent, not a completed estimate",
    "It does not prove",
)


def independent_mixture_check() -> list[str]:
    issues: list[str] = []
    weights = (Fraction(9, 10), Fraction(1, 10))
    scales = (Fraction(1, 4), Fraction(5, 2))
    moment1 = sum(weight * scale for weight, scale in zip(weights, scales))
    moment2 = sum(weight * scale**2 for weight, scale in zip(weights, scales))
    derivative = (moment2 - 3 * moment1**2) / 12
    if moment1 != Fraction(19, 40):
        issues.append(f"first scale moment mismatch: {moment1}")
    if moment2 != Fraction(109, 160):
        issues.append(f"second scale moment mismatch: {moment2}")
    if derivative != Fraction(7, 19200) or derivative <= 0:
        issues.append(f"Pick-wall derivative witness mismatch: {derivative}")
    return issues


def kernel_atom(u: float, z: complex) -> tuple[complex, complex]:
    root = cmath.sqrt(z)
    c_value = cmath.cosh(u * root)
    d_value = u * cmath.sinh(u * root) / (2 * root)
    return c_value, d_value


def independent_polarization_checks() -> list[str]:
    issues: list[str] = []
    atoms = (0.4, 1.1, 1.7)
    weights = (0.2, 0.5, 0.3)
    for z in (0.3 + 0.4j, 1.2 + 0.7j, -0.8 + 1.1j):
        values = [kernel_atom(u, z) for u in atoms]
        f_value = sum(weight * pair[0] for weight, pair in zip(weights, values))
        derivative = sum(weight * pair[1] for weight, pair in zip(weights, values))
        direct = -(derivative * f_value.conjugate()).imag
        polarized = 0.0
        for i, weight_i in enumerate(weights):
            c_i, d_i = values[i]
            for j, weight_j in enumerate(weights):
                c_j, d_j = values[j]
                kernel = -0.5 * (
                    d_i * c_j.conjugate() + d_j * c_i.conjugate()
                ).imag
                polarized += weight_i * weight_j * kernel
        if abs(direct - polarized) > 2e-13 * max(1.0, abs(direct)):
            issues.append(f"polarization mismatch at z={z}: {direct} vs {polarized}")
    return issues


def validate(artifact: Path, note: Path) -> list[str]:
    issues: list[str] = []
    if not artifact.exists():
        return [f"missing artifact: {artifact}"]
    if not note.exists():
        return [f"missing note: {note}"]

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if payload.get("kind") != "jensen_window_pf_phi_pick_kernel_target":
        issues.append("unexpected artifact kind")
    rows = payload.get("rows", [])
    if len(rows) != 12:
        issues.append(f"expected 12 rows, found {len(rows)}")
    ids = {row.get("id") for row in rows}
    if ids != REQUIRED_IDS:
        issues.append(f"row id mismatch: missing={sorted(REQUIRED_IDS-ids)} extra={sorted(ids-REQUIRED_IDS)}")

    row_by_id = {row.get("id"): row for row in rows}
    if row_by_id.get("ppkt_06_pick_endpoint_equivalence", {}).get("readiness") != "ready_to_apply":
        issues.append("Pick endpoint equivalence is not marked ready")
    if row_by_id.get("ppkt_10_two_scale_pick_failure", {}).get("readiness") != "guard_validated":
        issues.append("two-scale mixture guard is not validated")
    for row_id in ("ppkt_11_positive_resolvent_route", "ppkt_12_xi_phi_pick_handoff"):
        if row_by_id.get(row_id, {}).get("readiness") != "not_ready_to_apply":
            issues.append(f"open route is not marked not-ready: {row_id}")

    exact = payload.get("exact", {})
    for key in (
        "phi_transform",
        "log_derivative",
        "pick_numerator",
        "polarization",
        "pick_equivalence",
        "first_quadrant",
        "tilted_identity",
        "first_wall",
        "operator_route",
        "open_target",
    ):
        if not exact.get(key):
            issues.append(f"missing exact statement: {key}")

    sources = payload.get("sources", [])
    source_blob = json.dumps(sources)
    for marker in (
        "web.math.ku.dk/~berg/manus/castellon.pdf",
        "10.1016/j.jmaa.2022.126432",
        "jensen_window_pf_edrei_stieltjes_equivalence_gate.md",
    ):
        if marker not in source_blob:
            issues.append(f"missing source: {marker}")

    issues.extend(independent_mixture_check())
    issues.extend(independent_polarization_checks())

    text = note.read_text(encoding="utf-8")
    for marker in REQUIRED_NOTE_TEXT:
        if marker not in text:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore RH",
        "proves RH",
        "proves `Lambda <= 0`",
        "P_Phi(z)>0 has been proved",
        "the positive resolvent has been constructed",
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
        print(f"PHI-PICK-TARGET issue: {issue}")
    print(
        "validated Jensen-window PF Phi Pick-kernel target: "
        f"12 rows, {len(issues)} issues, 7 exact identities, "
        "3 independent polarization checks, 1 exact mixture guard, "
        "2 open structural routes"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
