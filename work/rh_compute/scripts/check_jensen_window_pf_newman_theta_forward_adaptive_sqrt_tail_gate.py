#!/usr/bin/env python3
"""Independently validate the cofinal ordinary-theta square-root tail gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SCRIPT_ROOT = REPO_ROOT / "work" / "rh_compute" / "scripts"
FINITE_STEM = "jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate"
FINITE_RESULT = RESULT_ROOT / f"{FINITE_STEM}.json"
FINITE_BUILDER = SCRIPT_ROOT / f"{FINITE_STEM}.py"
FINITE_CHECKER = SCRIPT_ROOT / f"check_{FINITE_STEM}.py"
MODULAR_RESULT = (
    RESULT_ROOT / "jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate.json"
)
PRECISION_BITS = 256
X_THRESHOLD = 245
K_THRESHOLD = 7
EXPECTED_IDS = [
    "ntfastg_01_ordinary_positive_partition",
    "ntfastg_02_all_k_gaussian_tail",
    "ntfastg_03_adaptive_sqrt_count",
    "ntfastg_04_gamma_conversion",
    "ntfastg_05_ceiling_geometry",
    "ntfastg_06_envelope_monotonicity",
    "ntfastg_07_cofinal_sqrt_tail_theorem",
    "ntfastg_08_retained_sqrt_handoff",
]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_ball(container: dict, key: str, issues: list[str]) -> arb:
    try:
        value = arb(container[key]["enclosure"])
    except Exception as exc:
        issues.append(f"{key} parse failed: {exc}")
        return arb(0)
    if value.rel_accuracy_bits() < 190:
        issues.append(f"{key} has weak stored accuracy")
    return value


def compare(
    stored: dict,
    key: str,
    expected: arb,
    issues: list[str],
) -> None:
    observed = parse_ball(stored, key, issues)
    if not observed.overlaps(expected):
        issues.append(f"{key} mismatch")


def independent_values() -> dict[str, arb]:
    x = arb(X_THRESHOLD)
    pi = arb.pi()
    g = x / 8 + 2 * (1 + x).log()
    root_g = g.sqrt()
    k = arb(K_THRESHOLD)
    theta = 1 - arb(121) / (80 * pi * k**2)
    rho = (arb(2) / k - pi * (2 * k + 1)).exp()
    denominator = theta * (1 - rho)
    coefficient = 20 * (arb(1) / 80).exp() * pi / denominator
    g_prime = arb(1) / 8 + 2 / (1 + x)
    ceiling_factor = 2 * g + 2
    return {
        "g_at_x245": g,
        "sqrt_g_at_x245": root_g,
        "theta_K7": theta,
        "rho_K7": rho,
        "denominator_K7": denominator,
        "coefficient_C": coefficient,
        "log_slope_margin_at_x245": (
            2 * (1 + x).log() + 1 - 2 * x / (1 + x)
        ),
        "ceiling_factor_at_x245": ceiling_factor,
        "envelope_decay_slope_at_x245": (
            2 * pi / (1 + x) - 5 / x
        ),
        "tuning_rate_margin": pi - arb(1) / 42 - 3,
        "value_ratio_envelope_at_x245": (
            coefficient
            * x**4
            * ceiling_factor
            / (1 + x) ** (2 * pi)
        ),
        "derivative_ratio_envelope_at_x245": (
            coefficient
            * x**3
            * (x + 4)
            * ceiling_factor
            / (1 + x) ** (2 * pi)
        ),
        "count_to_saddle_limit": (pi / 2).sqrt(),
    }


def symbolic_audit(issues: list[str]) -> None:
    x, h = sp.symbols("x h", positive=True)
    pi = sp.pi
    g = x / 8 + 2 * sp.log(1 + x)
    if sp.simplify(
        -pi * (g + h)
        + pi * x / 8
        + pi * h
        + 2 * pi * sp.log(1 + x)
    ) != 0:
        issues.append("gamma conversion identity failed")
    v_log_derivative = sp.diff(
        sp.log(x**4 * (2 * g + 2) / (1 + x) ** (2 * pi)),
        x,
    )
    expected_value = (
        4 / x + sp.diff(g, x) / (g + 1) - 2 * pi / (1 + x)
    )
    if sp.simplify(v_log_derivative - expected_value) != 0:
        issues.append("value envelope derivative failed")
    d_log_derivative = sp.diff(
        sp.log(
            x**3 * (x + 4) * (2 * g + 2)
            / (1 + x) ** (2 * pi)
        ),
        x,
    )
    expected = (
        3 / x
        + 1 / (x + 4)
        + sp.diff(g, x) / (g + 1)
        - 2 * pi / (1 + x)
    )
    if sp.simplify(d_log_derivative - expected) != 0:
        issues.append("derivative envelope derivative failed")
    log_slope_numerator = sp.simplify(g + 1 - x * sp.diff(g, x))
    expected_numerator = (
        2 * sp.log(1 + x) + 1 - 2 * x / (1 + x)
    )
    if sp.simplify(log_slope_numerator - expected_numerator) != 0:
        issues.append("g slope domination identity failed")


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "t_cap_exact": "1/5",
        "x_threshold": X_THRESHOLD,
        "adaptive_first_omitted": (
            "ceil(sqrt(x/8+2*log(1+x)+h))"
        ),
        "adaptive_retained": (
            "ceil(sqrt(x/8+2*log(1+x)+h))-1"
        ),
        "tuning_parameter": "h>=0",
        "first_omitted_threshold": K_THRESHOLD,
        "ratio_target_exact": "exp(-3h)/50",
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

    audit = artifact.get("source_audit", {})
    expected_hashes = {
        "finite_result_sha256": digest(FINITE_RESULT),
        "finite_builder_sha256": digest(FINITE_BUILDER),
        "finite_checker_sha256": digest(FINITE_CHECKER),
        "modular_comparison_result_sha256": digest(MODULAR_RESULT),
    }
    for key, expected in expected_hashes.items():
        if audit.get(key) != expected:
            issues.append(f"source hash mismatch: {key}")
    finite = json.loads(FINITE_RESULT.read_text(encoding="utf-8"))
    if audit.get("ordinary_tail_formula") != finite.get(
        "tail_majorant", {}
    ).get("formula"):
        issues.append("ordinary tail formula mismatch")

    expected_values = independent_values()
    stored = artifact.get("diagnostics", {})
    for key, expected in expected_values.items():
        compare(stored, key, expected, issues)
    if stored.get("K_at_x245") != K_THRESHOLD:
        issues.append("K at x=245 mismatch")
    if stored.get("ratio_target_exact") != "1/50":
        issues.append("stored ratio target mismatch")

    g = expected_values["g_at_x245"]
    if not (g.lower() > 36 and g.upper() < 49):
        issues.append("K=7 endpoint count guard failed")
    for key in (
        "theta_K7",
        "denominator_K7",
        "log_slope_margin_at_x245",
        "envelope_decay_slope_at_x245",
        "tuning_rate_margin",
    ):
        if expected_values[key].lower() <= 0:
            issues.append(f"nonpositive independent guard: {key}")
    target = arb(1) / 50
    for key in (
        "value_ratio_envelope_at_x245",
        "derivative_ratio_envelope_at_x245",
    ):
        if expected_values[key].upper() >= target:
            issues.append(f"1/50 endpoint target failed: {key}")

    symbolic_audit(issues)
    theorem = artifact.get("theorem", {})
    expected_theorem = {
        "domain": "0<=t<=1/5, x>=245",
        "value_error": (
            "|J_t-J_(N_(F,h),t)^F|<exp(-3h-pi*x/8)/50"
        ),
        "first_derivative_error": (
            "|J_t'-(J_(N_(F,h),t)^F)'|"
            "<exp(-3h-pi*x/8)/50"
        ),
        "saddle_count_comparison": (
            "For fixed h, N_(F,h)(x)/sqrt(x/(4*pi))->sqrt(pi/2)"
        ),
        "arbitrary_tolerance_corollary": (
            "For eta>0, h>=max(0,log(1/(50eta))/3) gives both "
            "errors <eta*exp(-pi*x/8)"
        ),
    }
    for key, expected in expected_theorem.items():
        if theorem.get(key) != expected:
            issues.append(f"theorem mismatch: {key}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "x>=245",
        "square-root count",
        "retained first-jet",
        "tunable h",
        "t=0",
        "Lambda<=0",
        "RH",
        "PF-infinity",
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
        "validated Newman theta forward adaptive square-root tail gate: "
        f"{len(artifact['rows'])} rows, 0 issues, tunable cofinal "
        "exp(-3h)/50 value/derivative tail with saddle-scale retained count, "
        "1 open retained-separation handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
