#!/usr/bin/env python3
"""Validate the prime-power heat-starlikeness base certificate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

from prime_power_heat_bernstein import build_base_certificates


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_prime_power_heat_"
    "starlikeness_base_certificate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "phase_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_prime_power_"
        "logarithmic_phase_flow_gate.json"
    ),
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complete_prime_power_"
        "chain_occupation_gate.json"
    ),
    "geometric_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_geometric_prime_power_"
        "phase_monotonicity_gate.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_direct_projection_"
        "regime_reduction.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_normalized_prefix_"
        "phase_flux_reduction.json"
    ),
}
EXPECTED_IDS = [
    "pphsb_01_domain",
    "pphsb_02_ideal_coordinates",
    "pphsb_03_current",
    "pphsb_04_chebyshev",
    "pphsb_05_bernstein",
    "pphsb_06_radicals",
    "pphsb_07_dyadic_box",
    "pphsb_08_dyadic_margin",
    "pphsb_09_ternary_box",
    "pphsb_10_ternary_margin",
    "pphsb_11_large_prime",
    "pphsb_12_quinary_audit",
    "pphsb_13_offset",
    "pphsb_14_coefficient_error",
    "pphsb_15_current_error",
    "pphsb_16_ray_error",
    "pphsb_17_actual_margins",
    "pphsb_18_boundary",
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


def exp_bounds_positive(
    value: Fraction, terms: int = 100
) -> tuple[Fraction, Fraction]:
    if value < 0:
        raise ValueError("positive exponential enclosure required")
    total = Fraction(1)
    term = Fraction(1)
    for index in range(1, terms + 1):
        term *= value / index
        total += term
    first_omitted = term * value / (terms + 1)
    ratio = value / (terms + 2)
    if ratio >= 1:
        raise ValueError("Taylor tail ratio is not contractive")
    return total, total + first_omitted / (1 - ratio)


def validate_current_identity(issues: list[str]) -> None:
    z = sp.symbols("z", nonzero=True)
    coefficients = sp.symbols("c0:8", positive=True, real=True)
    q_value = sum(
        coefficient * z**index
        for index, coefficient in enumerate(coefficients)
    )
    q_conjugate = sum(
        coefficient * z ** (-index)
        for index, coefficient in enumerate(coefficients)
    )
    h_value = sum(
        index * coefficient * z**index
        for index, coefficient in enumerate(coefficients)
    )
    h_conjugate = sum(
        index * coefficient * z ** (-index)
        for index, coefficient in enumerate(coefficients)
    )
    direct = sp.expand(
        q_value * q_conjugate
        + (h_value * q_conjugate + q_value * h_conjugate) / 2
    )
    paired = sum(
        (index + 1) * coefficients[index] ** 2
        for index in range(8)
    )
    paired += sum(
        sp.Rational(left + right + 2, 2)
        * coefficients[left]
        * coefficients[right]
        * (
            z ** (right - left)
            + z ** (left - right)
        )
        for left in range(8)
        for right in range(left + 1, 8)
    )
    if sp.expand(direct - paired) != 0:
        issues.append("Chebyshev current pairing identity failed")


def validate_parameter_boxes(issues: list[str]) -> None:
    specifications = (
        (2, Fraction(7, 10), Fraction(47, 50), Fraction(22, 25)),
        (3, Fraction(11, 10), Fraction(17, 20), Fraction(73, 100)),
        (5, Fraction(161, 100), Fraction(18, 25), Fraction(13, 25)),
    )
    cutoff_lower = 2**25
    for prime, log_upper, q_lower, s_lower in specifications:
        exp_log_lower, _ = exp_bounds_positive(log_upper)
        if not exp_log_lower > prime:
            issues.append(f"log({prime}) rational upper bound failed")

        q_exponent = log_upper * log_upper / 8
        _, q_exp_upper = exp_bounds_positive(q_exponent)
        if not q_exp_upper < 1 / q_lower:
            issues.append(f"p={prime} q-box containment failed")

        s_exponent = (
            log_upper
            * (log_upper + Fraction(1, cutoff_lower))
            / 4
        )
        _, s_exp_upper = exp_bounds_positive(s_exponent)
        if not s_exp_upper < 1 / s_lower:
            issues.append(f"p={prime} s-box containment failed")

    # log(2)>2/3 follows from the first positive atanh-series terms:
    # 2(1/3+(1/3)^3/3)>2/3.
    if not 2 * (Fraction(1, 3) + Fraction(1, 81)) > Fraction(
        2, 3
    ):
        issues.append("log(2) lower bound failed")


def validate_large_prime_family(issues: list[str]) -> None:
    # On 0<=r<=1/sqrt(5), f(r)=1-3r+2r^2 is decreasing because
    # 1/sqrt(5)<1/2. Its endpoint exceeds 1/20 precisely because
    # 1/sqrt(5)<9/20, whose square is the exact 80<81 comparison.
    if not Fraction(1, 5) < Fraction(1, 4):
        issues.append("large-prime derivative interval failed")
    if not Fraction(1, 5) < Fraction(81, 400):
        issues.append("large-prime endpoint margin failed")

    r = sp.symbols("r", nonnegative=True, real=True)
    theta = sp.symbols("theta", real=True)
    z = sp.exp(sp.I * theta)
    q_value = 1 + r * z
    h_value = r * z
    current = sp.simplify(
        sp.expand_complex(
            q_value * sp.conjugate(q_value)
            + sp.re(h_value * sp.conjugate(q_value))
        )
    )
    expected = 1 + 2 * r**2 + 3 * r * sp.cos(theta)
    if sp.simplify(sp.trigsimp(current - expected)) != 0:
        issues.append("two-level current identity failed")


def validate_actual_transfer(issues: list[str]) -> None:
    x_lower = 12 * 2**50
    d_upper = Fraction(2189, x_lower)
    if not d_upper < Fraction(1, 2):
        issues.append("normalized correction denominator bound failed")

    # For p=2, M=8 and p=3, M=4, the rational log upper bounds give
    # t*k*log(p)<7/2 and <33/20, respectively. For p>=5, M=2,
    # p<=N and log(N)<(L+1)/2 give t*log(p)<51/4<13.
    # In all cases eta=t*k*log(p)*delta_a/2<2/x.
    if not Fraction(7, 2) < 16:
        issues.append("dyadic saddle-offset exponent budget failed")
    if not Fraction(33, 20) < 16:
        issues.append("ternary saddle-offset exponent budget failed")
    if not Fraction(51, 4) < 13:
        issues.append("large-prime saddle-offset exponent budget failed")

    # Since exp(eta)-1<=eta/(1-eta), eta<2/x and x>4 imply
    # exp(eta)-1<4/x and exp(eta)<2. Together with
    # 2D/(1-D)<4D this gives |R_k-1|<18000/x.
    if not Fraction(2, x_lower - 2) < Fraction(4, x_lower):
        issues.append("saddle-offset exponential enclosure failed")
    r_upper = Fraction(18000, x_lower)
    if not r_upper < 1:
        issues.append("relative coefficient error is not contractive")
    if not Fraction(4 + 8 * 2189, x_lower) < r_upper:
        issues.append("relative coefficient perturbation budget failed")

    frozen_error = 576 * r_upper + 288 * r_upper**2
    if not frozen_error < 900 * r_upper:
        issues.append("frozen-current Lipschitz budget failed")

    ray_error = (
        Fraction(16_200_448, x_lower)
        + Fraction(6_486_528, x_lower**2)
    )
    if not ray_error < Fraction(1, 1000):
        issues.append("fixed-ray perturbation budget failed")

    margins = (
        Fraction(1, 300) - Fraction(1, 1000),
        Fraction(1, 25) - Fraction(1, 1000),
        Fraction(1, 20) - Fraction(1, 1000),
    )
    if margins != (
        Fraction(7, 3000),
        Fraction(39, 1000),
        Fraction(49, 1000),
    ):
        issues.append("actual current margins failed")


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

    phase = sources["phase_flow"].get("exact", {})
    if "J_ray=|Q|^2+Re(H*conj(Q))" not in phase.get(
        "phase_flow", {}
    ).get("ray_current", ""):
        issues.append("source ray current drifted")
    projection = sources["direct_projection"].get("exact", {})
    if "0<=delta_a<1/(4x)" not in projection.get(
        "saddle_bounds", ""
    ):
        issues.append("source saddle offset bound drifted")
    prefix = sources["normalized_prefix"].get("exact", {})
    if "|epsilon_n|<8446/x^2" not in prefix.get(
        "q_ge_1_inward", ""
    ):
        issues.append("source correction-current bound drifted")


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
    if stored.get("date") != "2026-07-30":
        issues.append("result date drifted")
    if [row.get("id") for row in stored.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    expected_summary = {
        "rows": 18,
        "exact_bernstein_base_certificates": 3,
        "analytic_large_prime_base_families": 1,
        "actual_ray_positive_base_families": 3,
        "length_propagation_theorems": 0,
        "proved_joined_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    rebuilt = [
        certificate.to_dict()
        for certificate in build_base_certificates()
    ]
    rebuilt_json = json.loads(json.dumps(rebuilt))
    certificates = stored.get("exact", {}).get(
        "bernstein_method", {}
    ).get("certificates", [])
    if certificates != rebuilt_json:
        issues.append("stored exact Bernstein certificates drifted")
    if len(rebuilt) != 3 or not all(
        certificate.get("positive") for certificate in rebuilt
    ):
        issues.append("one or more exact Bernstein certificates failed")

    exact = stored.get("exact", {})
    if "sum_(k<l)(k+l+2)" not in exact.get(
        "ideal_family", {}
    ).get("chebyshev_expansion", ""):
        issues.append("stored current polynomial drifted")
    if (
        "p^(-1/2)" not in exact.get("large_prime_family", {}).get(
            "reduction", ""
        )
        or "p>=5" not in stored.get("rows", [])[10].get("claim", "")
    ):
        issues.append("stored large-prime family drifted")
    if "<1/1000" not in exact.get("actual_transfer", {}).get(
        "ray_error", ""
    ):
        issues.append("stored ray-error bound drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove length propagation",
        "joined p-free",
        "Xi Abel gap",
        "successor winding bound",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for required in (
            "# Prime-Power Heat-Starlikeness Base Certificate",
            "## Exact Current",
            "## Exact Bernstein Certificates",
            "## Actual Fixed-Ray Transfer",
            "3 actual ray-positive base families",
            "0 length-propagation theorems",
        ):
            if required not in note:
                issues.append(f"note marker missing: {required}")

    validate_current_identity(issues)
    validate_parameter_boxes(issues)
    validate_large_prime_family(issues)
    validate_actual_transfer(issues)
    validate_sources(stored, sources, issues)
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"ERROR: {issue}")
    if issues:
        return 1
    print(
        "validated prime-power heat-starlikeness base certificate: "
        "18 rows, 3 exact Bernstein base certificates, "
        "1 analytic p>=5 two-level family, "
        "3 actual ray-positive base families, "
        "0 length-propagation theorems, 0 joined Abel gaps, "
        "0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
