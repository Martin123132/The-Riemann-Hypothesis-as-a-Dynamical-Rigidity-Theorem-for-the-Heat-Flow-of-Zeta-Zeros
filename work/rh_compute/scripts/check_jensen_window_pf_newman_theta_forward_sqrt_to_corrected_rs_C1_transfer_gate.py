#!/usr/bin/env python3
"""Independently validate the theta-to-corrected-RS C1 transfer gate."""

from __future__ import annotations

import argparse
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


STEM = (
    "jensen_window_pf_newman_theta_forward_"
    "sqrt_to_corrected_rs_C1_transfer_gate"
)
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SCRIPT_ROOT = REPO_ROOT / "work" / "rh_compute" / "scripts"
ADAPTIVE_STEM = (
    "jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate"
)
GLOBAL_STEM = (
    "jensen_window_pf_newman_polymath15_"
    "critical_C1_global_remainder_certificate"
)
CELL_STEM = (
    "jensen_window_pf_newman_polymath15_"
    "critical_C1_cell_remainder_certificate"
)
NORMALIZER_STEM = (
    "jensen_window_pf_newman_polymath15_normalized_laguerre_bridge"
)
PRECISION_BITS = 256
L_MIN = 50
X_COARSE = 245
EXPECTED_IDS = [
    "ntfsrsg_01_exact_normalizer_amplitude",
    "ntfsrsg_02_uniform_normalizer_lower_bound",
    "ntfsrsg_03_normalizer_first_derivative",
    "ntfsrsg_04_normalized_ordinary_tail",
    "ntfsrsg_05_direct_ordinary_transversality",
    "ntfsrsg_06_dual_approximation_identity",
    "ntfsrsg_07_corrected_rs_transfer_envelope",
    "ntfsrsg_08_scale_mismatch",
    "ntfsrsg_09_retained_arithmetic_handoff",
]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_paths(stem: str) -> tuple[Path, Path, Path]:
    return (
        RESULT_ROOT / f"{stem}.json",
        SCRIPT_ROOT / f"{stem}.py",
        SCRIPT_ROOT / f"check_{stem}.py",
    )


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
    pi = arb.pi()
    ell = arb(L_MIN)
    x = 4 * pi * ell.exp()
    coarse = arb(X_COARSE)
    coarse_ell = (coarse / (4 * pi)).log()
    scaled_factor = arb(1) / 2 + (1 + 4 / x) / ell
    x_power = ((arb(23) / 4) * x.log()).exp()
    theta_value_ratio = (3 * ell / 4).exp() / (25 * x_power)
    return {
        "x_at_L50": x,
        "normalizer_prefactor_margin_at_x245": (
            pi ** (arb(1) / 4)
            * (-arb(1) / (12 * coarse**2)).exp()
            - 1
        ),
        "Re_alpha_squared_floor_at_x245": (
            coarse_ell / 2 - 1 / coarse**2
        ) ** 2 - (pi / 4) ** 2,
        "Re_alpha_squared_floor_at_L50": (
            ell / 2 - 1 / x**2
        ) ** 2 - (pi / 4) ** 2,
        "log_A_derivative_bound_ratio_at_L50": (
            (ell / 4 + arb(1) / 2)
            * (1 + arb(1) / (2 * x))
            / (ell / 2)
        ),
        "scaled_theta_derivative_factor_at_L50": scaled_factor,
        "theta_value_envelope_over_exp_minus_3L_over_4_at_L50": (
            theta_value_ratio
        ),
        "theta_scaled_derivative_envelope_over_exp_minus_3L_over_4_at_L50": (
            scaled_factor * theta_value_ratio
        ),
        "composite_transfer_vector_constant": (
            arb(2501) ** 2 + arb(5001) ** 2
        ).sqrt(),
        "ordinary_to_corrected_squared_target_ratio_at_L50_h0": (
            (arb(2) / 625)
            / arb(32_000_000)
            / (((arb(23) / 2) * x.log()).exp())
            * (3 * ell / 2).exp()
        ),
    }


def symbolic_audit(issues: list[str]) -> None:
    x = sp.symbols("x", positive=True)
    delta_j = sp.Function("delta_j")(x)
    delta_h = delta_j / (16 * x**4)
    expected_delta_h_prime = (
        sp.diff(delta_j, x) / (16 * x**4)
        - delta_j / (4 * x**5)
    )
    if sp.simplify(sp.diff(delta_h, x) - expected_delta_h_prime) != 0:
        issues.append("characteristic-to-heat derivative identity failed")

    amplitude = sp.Function("amplitude")(x)
    error = sp.Function("error")(x)
    normalized = error / amplitude
    expected_normalized_prime = (
        sp.diff(error, x) / amplitude
        - sp.diff(amplitude, x) * error / amplitude**2
    )
    if sp.simplify(
        sp.diff(normalized, x) - expected_normalized_prime
    ) != 0:
        issues.append("normalized derivative identity failed")

    ordinary, corrected, r_rs, e_theta = sp.symbols(
        "ordinary corrected r_rs e_theta", real=True
    )
    if sp.expand(
        ordinary + e_theta - corrected - r_rs
        - (ordinary - corrected - r_rs + e_theta)
    ) != 0:
        issues.append("dual-approximation identity failed")


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    expected_domain = (
        "L=log(x/(4*pi))>=50, 0<t<=1/5, t*L<=25, h>=0"
    )
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "L_min": L_MIN,
        "t_cap_exact": "1/5",
        "scaled_time_cap": "t*L<=25",
        "tuning_parameter": "h>=0",
        "domain": expected_domain,
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append("parameter mismatch")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if [row.get("readiness") for row in rows] != (
        ["proved"] * 8 + ["not_ready_to_apply"]
    ):
        issues.append("readiness mismatch")

    audit = artifact.get("source_audit", {})
    for label, stem in (
        ("adaptive", ADAPTIVE_STEM),
        ("global", GLOBAL_STEM),
        ("cell", CELL_STEM),
        ("normalizer", NORMALIZER_STEM),
    ):
        paths = source_paths(stem)
        expected = {
            f"{label}_result_sha256": digest(paths[0]),
            f"{label}_builder_sha256": digest(paths[1]),
            f"{label}_checker_sha256": digest(paths[2]),
        }
        for key, value in expected.items():
            if audit.get(key) != value:
                issues.append(f"source hash mismatch: {key}")

    adaptive = json.loads(
        source_paths(ADAPTIVE_STEM)[0].read_text(encoding="utf-8")
    )
    global_gate = json.loads(
        source_paths(GLOBAL_STEM)[0].read_text(encoding="utf-8")
    )
    cell = json.loads(
        source_paths(CELL_STEM)[0].read_text(encoding="utf-8")
    )
    normalizer = json.loads(
        source_paths(NORMALIZER_STEM)[0].read_text(encoding="utf-8")
    )
    source_checks = {
        "adaptive_value_error": adaptive["theorem"]["value_error"],
        "adaptive_first_derivative_error": adaptive["theorem"][
            "first_derivative_error"
        ],
        "global_corrected_C1": global_gate["exact"]["global_c1"],
        "corrected_normalized_split": cell["exact"]["normalized_split"],
        "bulk_A_plus_B_envelope": cell["exact"]["eAB"],
        "published_M_t": normalizer["exact"]["published_input"]["M_t"],
        "published_alpha": normalizer["exact"]["published_input"]["alpha"],
    }
    for key, expected in source_checks.items():
        if audit.get(key) != expected:
            issues.append(f"source theorem mismatch: {key}")

    expected_values = independent_values()
    stored = artifact.get("diagnostics", {})
    for key, expected in expected_values.items():
        compare(stored, key, expected, issues)
    if stored.get("precision_bits") != PRECISION_BITS:
        issues.append("stored precision mismatch")
    if stored.get("L_min") != L_MIN:
        issues.append("stored L minimum mismatch")

    if expected_values[
        "normalizer_prefactor_margin_at_x245"
    ].lower() <= 0:
        issues.append("normalizer prefactor guard failed")
    if expected_values[
        "Re_alpha_squared_floor_at_x245"
    ].lower() <= 0:
        issues.append("coarse Re(alpha^2) guard failed")
    if expected_values["Re_alpha_squared_floor_at_L50"].lower() <= 0:
        issues.append("Re(alpha^2) guard failed")
    if expected_values[
        "log_A_derivative_bound_ratio_at_L50"
    ].upper() >= 1:
        issues.append("log-normalizer derivative guard failed")
    if expected_values[
        "scaled_theta_derivative_factor_at_L50"
    ].upper() >= arb(53) / 100:
        issues.append("scaled derivative factor guard failed")
    for key in (
        "theta_value_envelope_over_exp_minus_3L_over_4_at_L50",
        "theta_scaled_derivative_envelope_over_exp_minus_3L_over_4_at_L50",
        "ordinary_to_corrected_squared_target_ratio_at_L50_h0",
    ):
        if expected_values[key].upper() >= 1:
            issues.append(f"scale comparison failed: {key}")
    if expected_values[
        "composite_transfer_vector_constant"
    ].upper() >= 5600:
        issues.append("composite vector constant failed")

    symbolic_audit(issues)
    theorem = artifact.get("theorem", {})
    expected_theorem = {
        "domain": expected_domain,
        "ordinary_normalized_main": (
            "O_(h,t)=J_(N_(F,h),t)^F/(16*x^4*A_t)"
        ),
        "ordinary_value_error": (
            "E_F0=exp(-3h)/(25*x^(23/4))"
        ),
        "ordinary_scaled_first_derivative_error": (
            "E_F1=E_F0*(1/2+(1+4/x)/L)"
        ),
        "ordinary_direct_sufficient_target": (
            "T_L[O_h]>2*exp(-6h)/(625*x^(23/2))"
        ),
        "dual_approximation_identity": (
            "Z_t=O_(h,t)+e_F=J_hat_(N,t)+r_RS"
        ),
        "corrected_value_transfer": (
            "|O_(h,t)-J_hat_(N,t)|<2501*exp(-3L/4)"
        ),
        "corrected_scaled_derivative_transfer": (
            "|O_(h,t)'-J_hat_(N,t)'|/L"
            "<5001*exp(-3L/4)"
        ),
        "corrected_vector_transfer": (
            "||(O-J_hat),(O'-J_hat')/L||_2"
            "<5600*exp(-3L/4)"
        ),
    }
    for key, expected in expected_theorem.items():
        if theorem.get(key) != expected:
            issues.append(f"theorem mismatch: {key}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "L>=50",
        "0<t<=1/5",
        "tL<=25",
        "power-23/4",
        "Higher endpoint C_k terms alone",
        "not prove either retained lower bound",
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
    parser.add_argument(
        "--artifact", type=Path, default=DEFAULT_ARTIFACT
    )
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        return 1
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated theta square-root to corrected RS C1 transfer gate: "
        f"{len(artifact['rows'])} rows, 0 issues, exact dual "
        "approximation, direct ordinary threshold, and quantified "
        "certified-envelope mismatch"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
