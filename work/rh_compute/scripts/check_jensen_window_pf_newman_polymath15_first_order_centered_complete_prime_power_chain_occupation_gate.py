#!/usr/bin/env python3
"""Validate the complete-prime-power-chain occupation gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "complete_prime_power_chain_occupation_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "occupation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "signed_occupation_transport_reduction.json"
    ),
    "prime_power": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_block_composition_guard.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contact_signed_transport_reduction.json"
    ),
    "adjacent_chart": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
    ),
}


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def validate_sources(payload: dict) -> None:
    recorded = payload.get("source_hashes")
    require(isinstance(recorded, dict), "source_hashes missing")
    require(set(recorded) == set(SOURCES), "source key drift")
    for key, path in SOURCES.items():
        require(path.is_file(), f"missing source {key}")
        require(recorded[key] == file_hash(path), f"source hash drift: {key}")


def validate_chain_algebra() -> None:
    t, log_m, ell, sigma, tau = sp.symbols(
        "t log_m ell sigma tau", real=True
    )
    k = sp.symbols("k", integer=True, nonnegative=True)
    s_star = sigma + sp.I * tau
    lhs = (
        t * (log_m + k * ell) ** 2 / 4
        - s_star * (log_m + k * ell)
    )
    rhs = (
        t * log_m**2 / 4
        - s_star * log_m
        + k * ell * (t * log_m / 2 - s_star)
        + t * ell**2 * k**2 / 4
    )
    require(sp.expand(lhs - rhs) == 0, "chain factorization failed")

    c_s, b, kappa_0, u_0 = sp.symbols(
        "c_s b kappa_0 u_0", real=True
    )
    z_re, z_im, K_re, K_im = sp.symbols(
        "z_re z_im K_re K_im", real=True
    )
    direct = (
        c_s * (kappa_0 * z_re - ell * K_re)
        - b * (u_0 * z_im - ell * K_im)
    )
    projected = sp.re(
        (c_s * kappa_0 + sp.I * b * u_0)
        * (z_re + sp.I * z_im)
        - ell
        * (c_s + sp.I * b)
        * (K_re + sp.I * K_im)
    ).expand(complex=True)
    require(
        sp.simplify(direct - projected) == 0,
        "Euler moment failed",
    )


def validate_threshold_guard() -> None:
    ell = sp.symbols("ell", positive=True)
    ratio = sp.sqrt(2) / 2
    theta = sp.pi / 3
    masses = (
        sp.cos(theta),
        ratio * sp.cos(2 * theta),
        ratio**2 * sp.cos(3 * theta),
    )
    slopes = (
        sp.Rational(5, 2) * ell * sp.tan(theta),
        sp.Rational(3, 2) * ell * sp.tan(2 * theta),
        sp.Rational(1, 2) * ell * sp.tan(3 * theta),
    )
    expected_masses = (
        sp.Rational(1, 2),
        -sp.sqrt(2) / 4,
        -sp.Rational(1, 2),
    )
    expected_slopes = (
        5 * sp.sqrt(3) * ell / 2,
        -3 * sp.sqrt(3) * ell / 2,
        0,
    )
    for actual, expected in zip(masses, expected_masses):
        require(
            sp.simplify(actual - expected) == 0,
            "threshold mass failed",
        )
    for actual, expected in zip(slopes, expected_slopes):
        require(
            sp.simplify(actual - expected) == 0,
            "threshold slope failed",
        )
    ordered_profile = (
        sp.simplify(sum(masses)),
        sp.simplify(masses[0] + masses[2]),
        sp.simplify(masses[0]),
    )
    require(
        ordered_profile
        == (-sp.sqrt(2) / 4, 0, sp.Rational(1, 2)),
        "negative-zero-positive profile failed",
    )
    moment = sp.simplify(
        sum(mass * slope for mass, slope in zip(masses, slopes))
    )
    require(
        moment == ell * (10 * sp.sqrt(3) + 3 * sp.sqrt(6)) / 8,
        "threshold moment failed",
    )

    q = sp.Rational(2, 5)
    theta_0 = sp.acos(q)
    theta_1 = (sp.pi + theta_0) / 2
    theta_2 = sp.pi
    require(
        sp.simplify(sp.cos(theta_0) + q * sp.cos(theta_2)) == 0,
        "general endpoint cancellation failed",
    )
    require(
        float(sp.N(sp.cos(theta_1), 50)) < 0,
        "middle mass sign failed",
    )
    require(
        float(sp.N(sp.tan(theta_0), 50)) > 0
        and float(sp.N(sp.tan(theta_1), 50)) < -1,
        "general slope signs failed",
    )


def validate_outer_decomposition() -> None:
    for prime in (2, 3, 5, 7, 11):
        for n_max in range(1, 4097):
            middle = n_max // prime
            formula = n_max - 2 * middle + middle // prime
            actual = sum(
                1
                for base in range(middle + 1, n_max + 1)
                if base % prime != 0
            )
            require(formula == actual, "singleton formula failed")
    for n_max in range(1, 4097):
        middle = n_max // 2
        count = n_max - 2 * middle + middle // 2
        require(count >= n_max // 4, "dyadic singleton bound failed")

    q = sp.symbols("q")
    prefix = sp.symbols("S_0:9")
    telescoping = sp.expand(
        sum(
            q**k * (prefix[k] - q * prefix[k + 1])
            for k in range(8)
        )
    )
    require(
        telescoping == prefix[0] - q**8 * prefix[8],
        "p-free telescoping failed",
    )

    for n_max in range(1, 257):
        for prime in (2, 3, 5):
            represented = []
            for base in range(1, n_max + 1):
                if base % prime == 0:
                    continue
                value = base
                while value <= n_max:
                    represented.append(value)
                    value *= prime
            require(
                sorted(represented) == list(range(1, n_max + 1)),
                "p-free valuation partition failed",
            )


def validate_payload(payload: dict) -> None:
    require(
        payload.get("kind")
        == "complete_prime_power_chain_occupation_gate",
        "kind mismatch",
    )
    require(
        payload.get("status")
        == "exact_reduction_and_route_obstruction_not_an_rh_proof",
        "status mismatch",
    )
    rows = payload.get("rows")
    require(isinstance(rows, list) and len(rows) == 21, "row count mismatch")
    identifiers = [row.get("id") for row in rows]
    require(len(set(identifiers)) == 21, "duplicate row id")
    require(
        identifiers[0] == "cppco_01_coefficient_family"
        and identifiers[-1] == "cppco_21_proof_boundary",
        "row ordering mismatch",
    )
    text = json.dumps(payload, sort_keys=True)
    required = (
        "there is no Mobius factor",
        "B_m*a_(m,k)*w_(m,p)^k",
        "negative, zero on an open interval, and positive",
        "S_2(N)>=floor(N/4)",
        "p-free external phase",
        "first unresolved compact boundary beyond Q_207",
        "not an actual Xi point",
        "Clay-prize conclusion",
    )
    for marker in required:
        require(marker in text, f"payload marker missing: {marker}")
    summary = payload.get("summary", {})
    require(summary.get("row_count") == 21, "summary row count drift")
    require(
        summary.get("threshold_sign_obstructions") == 1,
        "threshold obstruction count drift",
    )
    require(
        summary.get("open_joined_boundary_theorems") == 1,
        "open theorem count drift",
    )


def validate_note() -> None:
    require(NOTE.is_file(), "note missing")
    text = NOTE.read_text(encoding="utf-8")
    markers = (
        "# Complete Prime-Power-Chain Occupation Gate",
        "## Coefficient Contract",
        "## Exact Chain Algebra",
        "## Complete-Chain Threshold Guard",
        "## P-Free Rejoining",
        "chainwise threshold positivity is false",
        "This is not a proof of RH",
    )
    for marker in markers:
        require(marker in text, f"note marker missing: {marker}")
    require(text.count("```") % 2 == 0, "unbalanced note fences")


def main() -> int:
    require(RESULT.is_file(), "result missing")
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    validate_sources(payload)
    validate_chain_algebra()
    validate_threshold_guard()
    validate_outer_decomposition()
    validate_payload(payload)
    validate_note()
    print(
        "validated complete prime-power-chain occupation gate: "
        "21 rows, exact chain factorization and Euler moment, "
        "negative-zero-positive threshold guard, external-phase "
        "reversal, order-N singleton theorem, p-free telescoping, "
        "joined-boundary handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
