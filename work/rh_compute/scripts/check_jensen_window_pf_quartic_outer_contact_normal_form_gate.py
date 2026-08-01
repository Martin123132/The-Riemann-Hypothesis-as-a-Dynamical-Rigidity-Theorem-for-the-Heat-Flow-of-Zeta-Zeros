#!/usr/bin/env python3
"""Check the exact quartic outer-contact normal-form gate."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_contact_normal_form_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_outer_contact_normal_form_gate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_quartic_outer_contact_normal_form_gate.py"
)


def load_builder():
    spec = importlib.util.spec_from_file_location("outer_contact_normal", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load outer-contact normal-form builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    builder = load_builder()
    rebuilt = builder.build_payload()
    if payload != rebuilt:
        issues.append("stored payload differs from exact reconstruction")
    if payload.get("kind") != (
        "jensen_window_pf_quartic_outer_contact_normal_form_gate"
    ):
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "exact outer-contact normal form with open contact-cell sign"
    ):
        issues.append("artifact status changed")
    if len(payload.get("rows", [])) != 13:
        issues.append("expected thirteen gate rows")

    delta, q, s = sp.symbols("delta q s")
    a = 1 - delta
    curvature = 4 * q * delta**2
    p = sp.factor(4 * a - 3 * a**2 + curvature)
    A = sp.factor(-3 * a**2 + 8 * a + p)
    B = sp.factor(-a**2 + 2 * a + p)
    x2 = sp.factor(A / 6)
    x3 = sp.factor(18 * a * B / A**2)
    x4 = sp.factor(2 * p * A / (3 * B**2))
    x5 = 1 - delta**2 * s
    U = sp.factor(B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2))
    T = sp.factor((1 - U) / delta**2)
    e2 = sp.factor((1 - x2) / delta**2)
    e3 = sp.factor((1 - x3) / delta**2)
    e4 = sp.factor((1 - x4) / delta**2)
    gap1 = sp.factor(
        (delta**2 * e3) ** 2
        - x3**2 * (delta**2 * e2) * (delta**2 * e4)
    )
    gap2 = sp.factor(
        (delta**2 * e4) ** 2
        - x4**2 * (delta**2 * e3) * (delta**2 * s)
    )
    K = sp.factor(A / 2)
    L = sp.factor(B / 2)
    if sp.factor(gap1 - delta**6 * 8 * q**3 / K**3) != 0:
        issues.append("independent G_1 scaling failed")
    h2 = sp.factor(gap2 / delta**4)
    remainder = sp.factor(8 * delta**2 * q**3 / (27 * L**3))
    if sp.factor(h2 - x4**2 * e3 * (T - s) - remainder) != 0:
        issues.append("independent G_2 decomposition failed")
    if sp.factor(x5 - U - delta**2 * (T - s)) != 0:
        issues.append("independent outer-threshold scaling failed")

    r, qr = sp.symbols("r qr", positive=True)
    right_poly = (
        4 * r**2 * qr**2
        - 12 * r**2 * qr
        + 9 * r**2
        + 6 * r * qr
        - 9 * r
        + 18 * qr
        - 18
    )
    if sp.factor(sp.diff(right_poly, qr, 2) - 8 * r**2) != 0:
        issues.append("right-branch polynomial is not convex")
    if sp.factor(right_poly.subs(qr, 0) - 9 * (r - 2) * (r + 1)) != 0:
        issues.append("right-branch q=0 endpoint changed")
    if sp.factor(right_poly.subs(qr, 1) - r * (r - 3)) != 0:
        issues.append("right-branch q=1 endpoint changed")

    adjacent = payload.get("exact", {}).get(
        "adjacent_quintic_discriminant", {}
    )
    q_factor = 15 * q**2 - 40 * q + 24
    expected_R0 = sp.factor(
        256 * delta**6 * q**4 * (16 * q - 15) * L**2
    )
    expected_R1 = sp.factor(
        96 * delta**3 * a**3 * q_factor * L * p**2
    )
    expected_R2 = sp.factor(9 * a**6 * p**4)
    expected_slack_discriminant = sp.factor(
        147456
        * delta**6
        * a**6
        * (6 - q) ** 2
        * (1 - q) ** 3
        * L**2
        * p**4
    )
    for key, expected_value in (
        ("R_0", expected_R0),
        ("R_1", expected_R1),
        ("R_2", expected_R2),
        ("quadratic_discriminant", expected_slack_discriminant),
    ):
        try:
            stored_value = sp.sympify(adjacent.get(key, ""))
        except (sp.SympifyError, TypeError):
            issues.append(f"could not parse adjacent-quintic {key}")
            continue
        if sp.factor(stored_value - expected_value) != 0:
            issues.append(f"adjacent-quintic {key} changed")
    witness_substitution = {
        delta: sp.Rational(13, 200),
        q: sp.Rational(1, 2),
    }
    witness_epsilon = sp.Rational(1, 100000)
    witness_R = sp.factor(
        (
            expected_R0
            + expected_R1 * witness_epsilon
            + expected_R2 * witness_epsilon**2
        ).subs(witness_substitution)
    )
    if witness_R >= 0:
        issues.append("rational survivor left the low-q discriminant collar")
    if adjacent.get("contact_phase_boundary") != "q=15/16":
        issues.append("adjacent-quintic contact phase boundary changed")

    summary = payload.get("summary", {})
    if summary != {
        "rows": 13,
        "contact_variables": 3,
        "excluded_outer_branches": 1,
        "normalized_defects": 4,
        "initial_gap_scalings": 2,
        "first_extension_regimes": 2,
        "normalized_terminal_compatibilities": 1,
        "adjacent_quintic_discriminant_factorizations": 1,
        "slack_quadratic_coefficients": 3,
        "low_q_negative_discriminant_collars": 1,
        "uniform_all_contact_theorems": 0,
    }:
        issues.append("summary counts changed")

    note = args.note.read_text(encoding="utf-8")
    for marker in (
        "exact contact-level reduction",
        "0<delta<1, 0<q<1",
        "x_3<x_2",
        "1-x_j=delta^2 e_j",
        "x_5-U=delta^2(T-s)",
        "G_1=delta^6 h_1",
        "G_2=delta^4 h_2",
        "Disc(P_5)=positive_prefactor*epsilon^2*R(epsilon)",
        "Disc_epsilon(R)=",
        "0<q<15/16",
        "0<epsilon<epsilon_+ => Disc(P_5)<0",
        "0<G_3<min(W_6,C_6)",
        "Y_6<=1",
        "Delta_14=delta^4 delta14_hat",
        "complete outer-contact obstruction",
        "`Lambda<=0`",
        "RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated quartic outer-contact normal-form gate: "
        "13 rows, 0 issues, 3 contact variables, "
        "1 excluded outer branch, 4 normalized defects, "
        "2 initial gap scalings, 2 first-extension regimes, "
        "1 normalized Delta_14, "
        "1 adjacent-quintic discriminant factorization, "
        "1 low-q negative-discriminant collar, "
        "0 uniform all-contact theorems"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
