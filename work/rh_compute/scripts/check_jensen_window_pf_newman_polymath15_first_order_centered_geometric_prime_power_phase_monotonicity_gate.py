#!/usr/bin/env python3
"""Validate the geometric prime-power phase-monotonicity gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_geometric_prime_power_"
    "phase_monotonicity_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "prime_power_heat": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_prime_power_heat_block_composition_guard.json"
    ),
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complete_prime_power_chain_occupation_gate.json"
    ),
    "ray_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_ray_bottom_logarithmic_flow_reduction.json"
    ),
}
EXPECTED_IDS = [
    "gppm_01_model",
    "gppm_02_factorization",
    "gppm_03_block_winding",
    "gppm_04_kernel",
    "gppm_05_phase_speed",
    "gppm_06_kernel_bounds",
    "gppm_07_monotone_margin",
    "gppm_08_projection_consequence",
    "gppm_09_length_monotonicity",
    "gppm_10_large_primes",
    "gppm_11_ternary_threshold",
    "gppm_12_dyadic_threshold",
    "gppm_13_complete_short_guard",
    "gppm_14_guard_count",
    "gppm_15_completion_verdict",
    "gppm_16_perturbation_target",
    "gppm_17_rejoin_guard",
    "gppm_18_handoff",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_phase_identity(issues: list[str]) -> None:
    a = sp.symbols("a", positive=True)
    theta = sp.symbols("theta", real=True)
    m = sp.symbols("M", positive=True, integer=True)
    def kernel(order: sp.Expr, radius: sp.Expr) -> sp.Expr:
        return order * (
            radius**2 - radius * sp.cos(order * theta)
        ) / (
            1
            - 2 * radius * sp.cos(order * theta)
            + radius**2
        )

    def factor_rate(order: sp.Expr, radius: sp.Expr) -> sp.Expr:
        factor = (
            1
            - radius * sp.cos(order * theta)
            - sp.I * radius * sp.sin(order * theta)
        )
        numerator = sp.expand_complex(
            sp.diff(factor, theta) * sp.conjugate(factor)
        )
        denominator = sp.expand_complex(factor * sp.conjugate(factor))
        return sp.simplify(sp.im(numerator) / denominator)

    numerator_rate = factor_rate(m, a**m)
    denominator_rate = factor_rate(sp.Integer(1), a)
    if sp.trigsimp(numerator_rate - kernel(m, a**m)) != 0:
        issues.append("numerator phase kernel failed")
    if sp.trigsimp(denominator_rate - kernel(1, a)) != 0:
        issues.append("denominator phase kernel failed")


def validate_thresholds(issues: list[str]) -> None:
    def mu(p: int, length: int) -> sp.Expr:
        a = 1 / sp.sqrt(p)
        return sp.simplify(
            1 / (1 + a)
            - length * a**length / (1 - a**length)
        )

    expected = {
        (5, 2): (3 - sp.sqrt(5)) / 4,
        (3, 4): (2 - sp.sqrt(3)) / 2,
        (2, 8): sp.Rational(22, 15) - sp.sqrt(2),
    }
    for key, value in expected.items():
        actual = mu(*key)
        if sp.simplify(actual - value) != 0:
            issues.append(f"threshold identity failed at {key}")
        if actual.is_positive is not True:
            issues.append(f"threshold positivity failed at {key}")

    a, length = sp.symbols("a M", positive=True)
    g_m = length * a**length / (1 - a**length)
    g_next = (length + 1) * a ** (length + 1) / (
        1 - a ** (length + 1)
    )
    numerator = sp.factor(
        sp.together(g_next - g_m).as_numer_denom()[0]
    )
    expected_numerator = a**length * (
        (length + 1) * a - length - a ** (length + 1)
    )
    if sp.simplify(numerator - expected_numerator) != 0:
        issues.append("tail monotonicity numerator failed")


def validate_short_guard(issues: list[str]) -> None:
    a, theta = sp.symbols("a theta", positive=True, real=True)
    x_value = -sp.sin(theta) - a * sp.sin(2 * theta)
    factorized = -sp.sin(theta) * (1 + 2 * a * sp.cos(theta))
    if sp.trigsimp(x_value - factorized) != 0:
        issues.append("short-chain projection factorization failed")

    derivative = sp.diff(x_value, theta)
    if sp.simplify(derivative.subs(theta, 0) + (1 + 2 * a)) != 0:
        issues.append("short-chain zero derivative at zero failed")
    if sp.simplify(derivative.subs(theta, sp.pi) - (1 - 2 * a)) != 0:
        issues.append("short-chain zero derivative at pi failed")

    cosine = -1 / (2 * a)
    sine_squared = 1 - cosine**2
    extra_rate = sp.simplify(
        derivative.subs(
            {
                sp.cos(theta): cosine,
                sp.cos(2 * theta): 2 * cosine**2 - 1,
            }
        )
    )
    if sp.simplify(extra_rate - 2 * a * sine_squared) != 0:
        issues.append("short-chain extra-zero orientation failed")

    for prime in (2, 3):
        if not (sp.Rational(1, 2) < 1 / sp.sqrt(prime) < 1):
            issues.append(f"short-chain domain failed for p={prime}")
    if not (1 / sp.sqrt(5) < sp.Rational(1, 2)):
        issues.append("large-prime short-chain separation failed")


def validate_sources(
    stored: dict, sources: dict[str, dict], issues: list[str]
) -> None:
    expected_hashes = {
        key: file_hash(path)
        for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored.get("source_audit", {}).get("source_sha256") != expected_hashes:
        issues.append("source hashes drifted")

    heat = sources["prime_power_heat"].get("exact", {})
    if "winding-three sum" not in heat.get("composition_guard", ""):
        issues.append("source composition guard drifted")
    chain = sources["complete_chain"]
    if "external phase B_m" not in chain.get(
        "threshold_guard", {}
    ).get("external_phase_flip", ""):
        issues.append("source external phase drifted")
    ray = sources["ray_flow"].get("exact", {})
    if "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)" not in (
        ray.get("ray_proxy", {}).get("defect_flux", "")
    ):
        issues.append("source ray defect drifted")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "gate result", issues)
    sources = {
        key: load_json(path, f"{key} source", issues)
        for key, path in SOURCES.items()
    }
    if issues:
        return issues

    if stored.get("kind") != STEM:
        issues.append("result kind drifted")
    if stored.get("date") != "2026-07-29":
        issues.append("result date drifted")
    if [row.get("id") for row in stored.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    expected_summary = {
        "rows": 18,
        "exact_phase_derivative_identities": 1,
        "uniform_monotone_prime_families": 3,
        "complete_short_chain_recrossing_guards": 1,
        "nonpromotion_guards": 2,
        "open_C1_perturbation_targets": 1,
        "proved_joined_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    exact = stored.get("exact", {})
    if exact.get("thresholds", {}).get("dyadic") != (
        "p=2, M>=8: mu_M>=22/15-sqrt(2)>0"
    ):
        issues.append("stored dyadic threshold drifted")
    if exact.get("short_chain_guard", {}).get("first_jet_winding") != (
        "wind(X_a+i*X_a')=-2"
    ):
        issues.append("stored short-chain winding drifted")
    if "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)" not in (
        exact.get("route", {}).get("open_handoff", "")
    ):
        issues.append("stored joined-defect handoff drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove C1 stability",
        "joined p-free",
        "Xi Abel gap",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for required in (
            "# Geometric Prime-Power Phase-Monotonicity Gate",
            "## Phase-Speed Margin",
            "## Complete Short-Chain Guard",
            "0 joined Abel gaps",
            "0 successor winding bounds",
        ):
            if required not in note:
                issues.append(f"note marker missing: {required}")

    validate_phase_identity(issues)
    validate_thresholds(issues)
    validate_short_guard(issues)
    validate_sources(stored, sources, issues)
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"ERROR: {issue}")
    if issues:
        return 1
    print(
        "validated geometric prime-power phase-monotonicity gate: "
        "18 rows, 1 phase-derivative identity, 3 uniform monotone "
        "prime families, 1 complete short-chain recrossing guard, "
        "2 nonpromotion guards, 1 open C1 perturbation target, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
