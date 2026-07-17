#!/usr/bin/env python3
"""Independently validate the weighted strong-log-concavity countermodel gate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate as target
import flint


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.md"
)

EXPECTED_IDS = [
    "nwslc_01_countermodel_family",
    "nwslc_02_strong_curvature",
    "nwslc_03_root_variable_concavity",
    "nwslc_04_explicit_admissible_kernel",
    "nwslc_05_base_fourier_transform",
    "nwslc_06_endpoint_laguerre_witness",
    "nwslc_07_explicit_arb_witness",
    "nwslc_08_weighted_correlation_failure",
    "nwslc_09_xi_specific_handoff",
]


def parse_ball(text: str) -> flint.arb:
    return flint.arb(text.replace("E", "e"))


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from reconstruction")
    if stored.get("kind") != (
        "jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate"
    ):
        issues.append("kind drifted")

    rows = stored.get("rows", [])
    if len(rows) != 9:
        issues.append("expected nine rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")

    u, x = sp.symbols("u x", real=True)
    polynomial = 1 + u**4 / 10
    curvature = sp.factor(sp.diff(sp.log(polynomial), u, 2))
    expected_curvature = sp.factor(
        4 * sp.Rational(1, 10) * u**2 * (3 - u**4 / 10)
        / (1 + u**4 / 10) ** 2
    )
    if sp.simplify(curvature - expected_curvature) != 0:
        issues.append("log-polynomial curvature identity failed")

    strong_constant = sp.simplify(2 - 6 / sp.sqrt(10))
    if not strong_constant.is_positive:
        issues.append("strong-curvature constant is not positive")
    z = sp.symbols("z", nonnegative=True)
    elementary_difference = sp.factor(
        3 * z / (1 + z**2)
        - z * (3 - z**2) / (1 + z**2) ** 2
    )
    if sp.simplify(elementary_difference - 4 * z**3 / (1 + z**2) ** 2) != 0:
        issues.append("elementary curvature majorant failed")
    if sp.simplify(
        sp.Rational(1, 2)
        - z / (1 + z**2)
        - (z - 1) ** 2 / (2 * (z**2 + 1))
    ) != 0:
        issues.append("z/(1+z^2) upper bound failed")

    gaussian = sp.sqrt(sp.pi) * sp.exp(-x**2 / 4)
    transform = sp.factor(gaussian + sp.diff(gaussian, x, 4) / 10)
    expected_transform = (
        sp.sqrt(sp.pi)
        * sp.exp(-x**2 / 4)
        * (x**4 - 12 * x**2 + 172)
        / 160
    )
    if sp.simplify(transform - expected_transform) != 0:
        issues.append("Gaussian transform identity failed")
    laguerre = sp.factor(
        sp.diff(transform, x) ** 2 - transform * sp.diff(transform, x, 2)
    )
    witness = sp.simplify(laguerre.subs(x, 3))
    expected_witness = -sp.Rational(743, 51200) * sp.pi * sp.exp(
        -sp.Rational(9, 2)
    )
    if sp.simplify(witness - expected_witness) != 0 or not witness.is_negative:
        issues.append("exact negative Laguerre witness failed")

    exact = stored.get("exact", {})
    r = sp.symbols("r", positive=True)
    root_log = -r - r**2 / 10 + sp.log(1 + r**2 / 10)
    root_second = sp.factor(sp.diff(root_log, r, 2))
    expected_root_second = -r**2 * (r**2 + 30) / (5 * (r**2 + 10) ** 2)
    if sp.simplify(root_second - expected_root_second) != 0:
        issues.append("explicit root-variable concavity identity failed")
    theta_root = sp.cosh(4 * sp.sqrt(r))
    theta_root_second = sp.factor(sp.diff(theta_root, r, 2))
    expected_theta_root_second = (
        4 * sp.sqrt(r) * sp.cosh(4 * sp.sqrt(r))
        - sp.sinh(4 * sp.sqrt(r))
    ) / r ** sp.Rational(3, 2)
    if sp.simplify(theta_root_second - expected_theta_root_second) != 0:
        issues.append("theta-factor root-variable convexity identity failed")
    y = sp.symbols("y", positive=True)
    theta_numerator = 4 * y * sp.cosh(4 * y) - sp.sinh(4 * y)
    theta_numerator_derivative = sp.factor(sp.diff(theta_numerator, y))
    if sp.simplify(theta_numerator_derivative - 16 * y * sp.sinh(4 * y)) != 0:
        issues.append("theta-factor positive numerator derivative failed")
    root_data = exact.get("explicit_root_concavity", {})
    if root_data.get("base_second_derivative") != str(root_second):
        issues.append("stored base root-variable concavity formula drifted")
    stored_theta_root_second = sp.sympify(
        root_data.get("theta_factor_second_derivative", "nan"), locals={"r": r}
    )
    if sp.simplify(stored_theta_root_second - theta_root_second) != 0:
        issues.append("stored theta-factor root-variable formula drifted")
    if (
        root_data.get("theta_factor_positive_numerator_derivative")
        != str(theta_numerator_derivative)
    ):
        issues.append("stored theta-factor positivity formula drifted")

    certificate = stored.get("certificate", {})
    fresh_certificate = target.explicit_countermodel_certificate()
    if certificate != fresh_certificate:
        issues.append("stored Arb certificate differs from fresh computation")
    for field in (
        "F_ball",
        "F_prime_ball",
        "F_second_ball",
        "laguerre_ball",
        "laguerre_upper",
    ):
        if field not in certificate:
            issues.append(f"certificate missing {field}")
    if certificate.get("delta") != "1/10":
        issues.append("explicit delta drifted")
    if certificate.get("epsilon") != "1/1000":
        issues.append("explicit epsilon drifted")
    if certificate.get("x") != "21/5":
        issues.append("witness x drifted")
    if certificate.get("correlation_frequency") != "42/5":
        issues.append("correlation frequency drifted")
    if certificate.get("passed") is not True:
        issues.append("certificate not marked passed")
    if parse_ball(certificate.get("laguerre_ball", "0")) >= 0:
        issues.append("theta-tail Laguerre ball did not certify negative")
    if parse_ball(certificate.get("F_ball", "0")).contains(0):
        issues.append("F witness ball unexpectedly contains zero")

    weighted = exact.get("weighted_correlation", {})
    for phrase in (
        "s^2",
        "L_1[F_(delta,epsilon)](x)=4*Fourier",
        "42/5",
        "<0",
    ):
        joined = " ".join(str(value) for value in weighted.values())
        if phrase not in joined:
            issues.append(f"weighted correlation gate missing {phrase}")
    decision = exact.get("nonpromotion_decision", "")
    for phrase in (
        "strict concavity of log(phi(sqrt(r)))",
        "do not imply",
        "K_(1)",
        "arithmetic or modular",
    ):
        if phrase not in decision:
            issues.append(f"nonpromotion decision missing {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    for phrase in (
        "not a proof of RH",
        "Countermodel Family",
        "Exact Endpoint Witness",
        "Weighted-Correlation Failure",
        "-743*pi*exp(-9/2)/51200",
        "delta=1/10",
        "epsilon=1/1000",
        "21/5",
        "42/5",
        "root-variable",
        "theta-type",
        "s^2",
        "https://arxiv.org/abs/1309.0055",
    ):
        if phrase not in note:
            issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF Newman weighted strong-log-concavity "
        "countermodel gate: 9 rows, 0 issues, 1 exact strong-curvature bound, "
        "1 exact root-variable concavity theorem, "
        "1 explicit theta-tail admissible kernel, "
        "1 Gaussian endpoint identity, 1 exact endpoint witness, "
        "1 Arb theta-tail witness, "
        "1 weighted-correlation countermodel, 1 Xi-specific handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
