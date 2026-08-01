#!/usr/bin/env python3
"""Independently validate the six-term ordinary-theta bridge tail theorem."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SCRIPT_ROOT = REPO_ROOT / "work" / "rh_compute" / "scripts"
MODULAR_SOURCE = RESULT_ROOT / "jensen_window_pf_newman_theta_modular_blend_gate.json"
BASE_STEM = "jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate"
BASE_RESULT = RESULT_ROOT / f"{BASE_STEM}.json"
BASE_BUILDER = SCRIPT_ROOT / f"{BASE_STEM}.py"
BASE_CHECKER = SCRIPT_ROOT / f"check_{BASE_STEM}.py"
CONTRACT_RESULT = (
    RESULT_ROOT
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)
PRECISION_BITS = 256
X_MIN = 38
X_MAX = 245
RETAINED_TERMS = 6
FIRST_OMITTED = 7
WITNESS_FIRST_OMITTED = (6, 7, 8)
RATIO_TARGET = Fraction(1, 100_000_000_000)
EXPECTED_IDS = [
    "ntfstbtg_01_forward_half_line_partition",
    "ntfstbtg_02_six_term_transform",
    "ntfstbtg_03_raw_forward_tail_majorant",
    "ntfstbtg_04_direct_c1_error",
    "ntfstbtg_05_bridge_endpoint_monotonicity",
    "ntfstbtg_06_six_term_gamma_scale_theorem",
    "ntfstbtg_07_finite_bridge_reduction",
    "ntfstbtg_08_retained_bridge_handoff",
]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def independent_tail(moment: int, first: int) -> arb:
    pi = arb.pi()
    theta = 1 - (arb(121) / 80) / (pi * first**2)
    rho = (arb(2) / first - pi * (2 * first + 1)).exp()
    if theta.lower() <= 0 or rho.upper() >= 1:
        raise RuntimeError(f"invalid tail denominator at K={first}")
    return (
        arb(5)
        * factorial(moment)
        * (arb(1) / 80).exp()
        * pi
        * first**2
        * (-pi * first**2).exp()
        / (4 * theta * (1 - rho))
    )


def independent_witness(first: int) -> dict[str, arb]:
    e0 = independent_tail(0, first)
    e1 = independent_tail(1, first)
    x = arb(X_MAX)
    gamma = (-arb.pi() * x / 8).exp()
    value_error = 16 * x**4 * e0
    derivative_error = 64 * x**3 * e0 + 16 * x**4 * e1
    return {
        "E0": e0,
        "E1": e1,
        "gamma_scale": gamma,
        "J_error": value_error,
        "J_prime_error": derivative_error,
        "J_error_over_gamma": value_error / gamma,
        "J_prime_error_over_gamma": derivative_error / gamma,
    }


def compare_ball(
    issues: list[str],
    label: str,
    stored: dict,
    expected: arb,
) -> None:
    try:
        observed = arb(stored["enclosure"])
    except Exception as exc:
        issues.append(f"{label} parse failed: {exc}")
        return
    if not observed.overlaps(expected):
        issues.append(f"{label} mismatch")
    if observed.rel_accuracy_bits() < 200:
        issues.append(f"{label} weak stored accuracy")


def symbolic_audit(issues: list[str]) -> None:
    x = sp.symbols("X")
    polynomial = sp.Poly(2 * x - 3, x)
    coefficient_norm = sum(abs(value) for value in polynomial.all_coeffs())
    if coefficient_norm != 5:
        issues.append("P_0 coefficient norm drifted")

    beta = sp.Rational(9, 4)
    g = beta - sp.Rational(3, 4)
    if g + sp.Rational(1, 80) != sp.Rational(121, 80):
        issues.append("forward theta denominator shift drifted")

    frequency, rate, power = sp.symbols(
        "frequency rate power",
        positive=True,
    )
    multiplier = sp.exp(rate * frequency) * frequency**power
    logarithmic_derivative = sp.simplify(
        sp.diff(multiplier, frequency) / multiplier
    )
    if sp.simplify(logarithmic_derivative - (rate + power / frequency)) != 0:
        issues.append("gamma-normalized multiplier derivative failed")


def source_audit(artifact: dict, issues: list[str]) -> None:
    audit = artifact.get("source_audit", {})
    expected_hashes = {
        "modular_source_sha256": digest(MODULAR_SOURCE),
        "base_result_sha256": digest(BASE_RESULT),
        "base_builder_sha256": digest(BASE_BUILDER),
        "base_checker_sha256": digest(BASE_CHECKER),
        "contract_result_sha256": digest(CONTRACT_RESULT),
    }
    for key, expected in expected_hashes.items():
        if audit.get(key) != expected:
            issues.append(f"source hash mismatch: {key}")

    modular = json.loads(MODULAR_SOURCE.read_text(encoding="utf-8"))
    theta = modular.get("exact", {}).get("theta_summand", {})
    expected_definition = (
        "phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))"
        "*exp(-pi*n^2*exp(4u))"
    )
    if theta.get("definition") != expected_definition:
        issues.append("ordinary theta summand definition mismatch")
    if theta.get("kernel_series") != "Phi(u)=sum_(n>=1)phi_n(u)=Phi(-u)":
        issues.append("ordinary theta series mismatch")
    positivity = (
        modular.get("exact", {})
        .get("strict_positivity", {})
        .get("positive_side")
    )
    if positivity != "phi_n(u)>0 for n>=1 and u>=0":
        issues.append("positive-half-line marker mismatch")

    contract = json.loads(CONTRACT_RESULT.read_text(encoding="utf-8"))
    direct = contract.get("exact", {}).get("direct_c1_contract", {})
    if direct.get("partial_value") != "J_(N,t)(x)=16*x^4*S_(N,t)(x)":
        issues.append("direct value contract mismatch")
    if direct.get("partial_derivative") != (
        "J_(N,t)'(x)=64*x^3*S_(N,t)(x)+16*x^4*S_(N,t)'(x)"
    ):
        issues.append("direct derivative contract mismatch")


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "t_cap_exact": "1/5",
        "x_interval": [X_MIN, X_MAX],
        "retained_terms": RETAINED_TERMS,
        "first_omitted": FIRST_OMITTED,
        "witness_first_omitted": list(WITNESS_FIRST_OMITTED),
        "ratio_target_exact": "1/100000000000",
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append("parameter mismatch")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if [row.get("readiness") for row in rows] != (
        ["proved"] * 7 + ["not_ready_to_apply"]
    ):
        issues.append("readiness mismatch")

    source_audit(artifact, issues)
    symbolic_audit(issues)

    witnesses = artifact.get("witnesses", [])
    if [row.get("first_omitted") for row in witnesses] != list(
        WITNESS_FIRST_OMITTED
    ):
        issues.append("witness sequence mismatch")
    rebuilt: dict[int, dict[str, arb]] = {}
    for stored in witnesses:
        first = int(stored.get("first_omitted", -1))
        if stored.get("retained_terms") != first - 1:
            issues.append(f"K={first} retained count mismatch")
        expected = independent_witness(first)
        rebuilt[first] = expected
        for key, value in expected.items():
            compare_ball(issues, f"K={first} {key}", stored.get(key, {}), value)

    target = arb(RATIO_TARGET.numerator) / RATIO_TARGET.denominator
    for key in ("J_error_over_gamma", "J_prime_error_over_gamma"):
        if rebuilt.get(7, {}).get(key, arb(1)).upper() >= target:
            issues.append(f"K=7 ratio target failed: {key}")
        if rebuilt.get(6, {}).get(key, arb(0)).lower() <= 1:
            issues.append(f"K=6 failure guard failed: {key}")
        if (
            rebuilt.get(8, {}).get(key, arb(1)).upper()
            >= rebuilt.get(7, {}).get(key, arb(0)).lower()
        ):
            issues.append(f"K=8 decrease guard failed: {key}")

    monotonicity = artifact.get("monotonicity", {})
    rate = arb.pi() / 8
    expected_diagnostics = {
        "spectral_rate_pi_over_8": rate,
        "value_x4_log_slope_at_x38": rate + arb(4) / X_MIN,
        "derivative_x3_log_slope_at_x38": rate + arb(3) / X_MIN,
        "derivative_x4_log_slope_at_x38": rate + arb(4) / X_MIN,
        "theta_K7": 1 - (arb(121) / 80) / (arb.pi() * 7**2),
        "rho_K7": (arb(2) / 7 - arb.pi() * 15).exp(),
    }
    for key, value in expected_diagnostics.items():
        compare_ball(issues, f"diagnostic {key}", monotonicity.get(key, {}), value)
        if value.lower() <= 0:
            issues.append(f"diagnostic is not positive: {key}")
    if monotonicity.get("monotonicity_formula") != (
        "d/dx log(exp(pi*x/8)*x^s)=pi/8+s/x>0, "
        "x>0, s in {3,4}"
    ):
        issues.append("monotonicity formula mismatch")

    tail = artifact.get("tail_majorant", {})
    if tail.get("kernel_polynomial") != "P_0(X)=2X-3":
        issues.append("tail kernel polynomial mismatch")
    if tail.get("kernel_coefficient_norm") != 5:
        issues.append("tail coefficient norm mismatch")
    if tail.get("beta") != "9/4" or tail.get("kernel_power") != 2:
        issues.append("tail template mismatch")

    theorem = artifact.get("theorem", {})
    if theorem.get("domain") != "0<=t<=1/5, 38<=x<=245":
        issues.append("theorem domain mismatch")
    if theorem.get("value_error") != (
        "|J_t-J_(6,t)^F|<10^-11*exp(-pi*x/8)"
    ):
        issues.append("value theorem mismatch")
    if theorem.get("first_derivative_error") != (
        "|J_t'-(J_(6,t)^F)'|<10^-11*exp(-pi*x/8)"
    ):
        issues.append("first-derivative theorem mismatch")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "38<=x<=245",
        "retained first-jet separation",
        "finite bridge",
        "cofinal retained",
        "strict Laguerre",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        return 1
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman theta forward six-term bridge tail gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['witnesses'])} witnesses, "
        "0 issues, direct gamma-scale tail theorem on 38<=x<=245, "
        "1 open retained-cover handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
