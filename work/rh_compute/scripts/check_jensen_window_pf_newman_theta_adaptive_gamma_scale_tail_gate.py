#!/usr/bin/env python3
"""Independently validate the adaptive gamma-scale modular-tail theorem."""

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

import check_jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate as base_check


STEM = "jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate"
BASE_STEM = "jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
BASE_RESULT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{BASE_STEM}.json"
)
BASE_BUILDER = (
    REPO_ROOT / "work" / "rh_compute" / "scripts" / f"{BASE_STEM}.py"
)
BASE_CHECKER = (
    REPO_ROOT / "work" / "rh_compute" / "scripts" / f"check_{BASE_STEM}.py"
)
CONTRACT_RESULT = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)
PRECISION_BITS = 256
N_THRESHOLD = 63
X_THRESHOLD = 245
BOUNDARY_N = (62, 63, 64)
EXPECTED_IDS = [
    "ntagstg_01_adaptive_cell_partition",
    "ntagstg_02_within_cell_monotonicity",
    "ntagstg_03_forward_endpoint_monotonicity",
    "ntagstg_04_compact_endpoint_monotonicity",
    "ntagstg_05_large_endpoint_monotonicity",
    "ntagstg_06_positive_composition",
    "ntagstg_07_cofinal_gamma_scale_theorem",
    "ntagstg_08_retained_margin_handoff",
]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def power(base: arb, exponent: Fraction) -> arb:
    return (
        arb(exponent.numerator)
        / exponent.denominator
        * base.log()
    ).exp()


def independently_rebuild_boundary(
    retained: int,
    kernels: list[dict],
    switches: list[dict],
) -> dict[str, arb]:
    compiled = base_check.independently_compile(
        retained + 1,
        kernels,
        switches,
    )
    d0 = compiled["d0"]["total"]
    d1 = compiled["d1"]["total"]
    x = power(arb(retained), Fraction(4, 3)) - 1
    gamma = (-arb.pi() * x / 8).exp()
    value_error = 16 * d0 / x**5
    derivative_error = 64 * d0 / x**6 + 16 * d1 / x**5
    return {
        "x_right": x,
        "d0": d0,
        "d1": d1,
        "gamma_scale": gamma,
        "J_error": value_error,
        "J_prime_error": derivative_error,
        "J_error_over_gamma": value_error / gamma,
        "J_prime_error_over_gamma": derivative_error / gamma,
    }


def compare_boundary(
    issues: list[str],
    stored: dict,
    rebuilt: dict[str, arb],
) -> None:
    retained = stored.get("N")
    for key, expected in rebuilt.items():
        try:
            observed = arb(stored[key]["enclosure"])
        except Exception as exc:
            issues.append(f"N={retained} {key} parse failed: {exc}")
            continue
        if not observed.overlaps(expected):
            issues.append(f"N={retained} {key} mismatch")
        if observed.rel_accuracy_bits() < 200:
            issues.append(f"N={retained} {key} weak stored accuracy")


def independent_monotonicity() -> dict[str, arb]:
    pi = arb.pi()
    rate = pi / 8
    c_safe = (
        arb(3)
        / 16
        * power(arb(3), Fraction(2, 3))
        * power(pi, Fraction(2, 3))
    )
    n = arb(N_THRESHOLD)
    delta = (
        power(arb(N_THRESHOLD + 1), Fraction(4, 3))
        - power(n, Fraction(4, 3))
    )
    return {
        "reflected_rate_margin": c_safe - rate,
        "forward_log": (
            arb(20) / (N_THRESHOLD + 1)
            - pi * (2 * (N_THRESHOLD + 1) + 1)
            + rate * delta
        ),
        "compact_log": (
            arb(22) / (N_THRESHOLD + 1)
            - pi / arb(2).sqrt() * (2 * (N_THRESHOLD + 1) + 1)
            + rate * delta
        ),
        "large_margin": (
            arb(4)
            * (c_safe - rate)
            / 3
            * power(n, Fraction(1, 3))
            - arb(22) / n
        ),
        "value_slope": rate - arb(5) / X_THRESHOLD,
        "derivative_slope": rate - arb(6) / X_THRESHOLD,
    }


def compare_diagnostic(
    issues: list[str],
    stored: dict,
    key: str,
    expected: arb,
) -> None:
    try:
        observed = arb(stored[key]["enclosure"])
    except Exception as exc:
        issues.append(f"diagnostic {key} parse failed: {exc}")
        return
    if not observed.overlaps(expected):
        issues.append(f"diagnostic {key} mismatch")


def symbolic_monotonicity_audit(issues: list[str]) -> None:
    h, rate_prime = sp.symbols("h rate_prime", positive=True)
    expression = sp.factor(
        rate_prime * (1 + 1 / h) - (h + 1)
    )
    expected = sp.factor((h + 1) * (rate_prime / h - 1))
    if sp.simplify(expression - expected) != 0:
        issues.append("large-tail derivative factorization failed")

    n = sp.symbols("n", positive=True)
    c, a, degree = sp.symbols("c a degree", positive=True)
    phase_margin = (
        sp.Rational(4, 3) * (c - a) * n ** sp.Rational(1, 3)
        - degree / n
    )
    derivative = sp.diff(phase_margin, n)
    expected_derivative = (
        sp.Rational(4, 9) * (c - a) / n ** sp.Rational(2, 3)
        + degree / n**2
    )
    if sp.simplify(derivative - expected_derivative) != 0:
        issues.append("large-tail phase-margin monotonicity failed")


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "m_order": 9,
        "t_cap_exact": "1/5",
        "adaptive_count": "ceil((1+x)^(3/4))",
        "x_threshold": X_THRESHOLD,
        "n_threshold": N_THRESHOLD,
        "ratio_target_exact": "1/4",
        "boundary_n": list(BOUNDARY_N),
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
        "base_result_sha256": digest(BASE_RESULT),
        "base_builder_sha256": digest(BASE_BUILDER),
        "base_checker_sha256": digest(BASE_CHECKER),
        "contract_result_sha256": digest(CONTRACT_RESULT),
    }
    for key, expected in expected_hashes.items():
        if audit.get(key) != expected:
            issues.append(f"source hash mismatch: {key}")

    base_issues = base_check.validate(BASE_RESULT)
    issues.extend(f"base gate: {item}" for item in base_issues)
    kernels = base_check.independent_kernel_rows()
    switches = base_check.independent_switch_rows()
    stored_boundaries = artifact.get("boundary_witnesses", [])
    if [row.get("N") for row in stored_boundaries] != list(BOUNDARY_N):
        issues.append("boundary witness sequence mismatch")
    rebuilt_by_n: dict[int, dict[str, arb]] = {}
    for stored in stored_boundaries:
        retained = int(stored["N"])
        rebuilt = independently_rebuild_boundary(
            retained,
            kernels,
            switches,
        )
        rebuilt_by_n[retained] = rebuilt
        compare_boundary(issues, stored, rebuilt)

    target = arb(1) / 4
    for key in ("J_error_over_gamma", "J_prime_error_over_gamma"):
        if rebuilt_by_n.get(63, {}).get(key, arb(1)).upper() >= target:
            issues.append(f"N=63 target failed: {key}")
        if rebuilt_by_n.get(62, {}).get(key, arb(0)).lower() <= 1:
            issues.append(f"N=62 non-promotion guard failed: {key}")
        if (
            rebuilt_by_n.get(64, {}).get(key, arb(1)).upper()
            >= rebuilt_by_n.get(63, {}).get(key, arb(0)).lower()
        ):
            issues.append(f"N=64 endpoint decrease failed: {key}")

    monotonicity = independent_monotonicity()
    if monotonicity["reflected_rate_margin"].lower() <= 0:
        issues.append("reflected safe rate does not beat pi/8")
    if monotonicity["forward_log"].upper() >= 0:
        issues.append("forward normalized ratio check failed")
    if monotonicity["compact_log"].upper() >= 0:
        issues.append("compact normalized ratio check failed")
    if monotonicity["large_margin"].lower() <= 0:
        issues.append("large normalized tail check failed")
    if min(
        monotonicity["value_slope"].lower(),
        monotonicity["derivative_slope"].lower(),
    ) <= 0:
        issues.append("within-cell monotonicity check failed")

    stored_diagnostics = artifact.get("monotonicity", {})
    compare_diagnostic(
        issues,
        stored_diagnostics,
        "reflected_rate_margin",
        monotonicity["reflected_rate_margin"],
    )
    compare_diagnostic(
        issues,
        stored_diagnostics,
        "worst_forward_log_ratio_upper_at_N63",
        monotonicity["forward_log"],
    )
    compare_diagnostic(
        issues,
        stored_diagnostics,
        "worst_compact_log_ratio_upper_at_N63",
        monotonicity["compact_log"],
    )
    compare_diagnostic(
        issues,
        stored_diagnostics,
        "worst_large_derivative_margin_at_N63",
        monotonicity["large_margin"],
    )
    compare_diagnostic(
        issues,
        stored_diagnostics,
        "value_within_cell_log_slope_at_x245",
        monotonicity["value_slope"],
    )
    compare_diagnostic(
        issues,
        stored_diagnostics,
        "derivative_within_cell_log_slope_at_x245",
        monotonicity["derivative_slope"],
    )
    symbolic_monotonicity_audit(issues)

    theorem = artifact.get("theorem", {})
    if theorem.get("domain") != "0<=t<=1/5, x>=245":
        issues.append("theorem domain mismatch")
    if theorem.get("value_error") != (
        "|J_t-J_(N,t)|<exp(-pi*x/8)/4"
    ):
        issues.append("value theorem mismatch")
    if theorem.get("first_derivative_error") != (
        "|J_t'-J_(N,t)'|<exp(-pi*x/8)/4"
    ):
        issues.append("first-derivative theorem mismatch")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "x>=245",
        "retained J_N",
        "finite bridge 69<x<245",
        "terminating",
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
        "validated Newman theta adaptive gamma-scale tail gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['boundary_witnesses'])} boundary witnesses, "
        "0 issues, cofinal x>=245 omitted-tail theorem, "
        "1 open retained-margin handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
