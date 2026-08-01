#!/usr/bin/env python3
"""Validate the prime-power logarithmic phase-flow gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import isqrt
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_prime_power_"
    "logarithmic_phase_flow_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complete_prime_power_"
        "chain_occupation_gate.json"
    ),
    "ray_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_ray_bottom_logarithmic_flow_reduction.json"
    ),
    "geometric_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_geometric_prime_power_"
        "phase_monotonicity_gate.json"
    ),
    "prime_power_heat": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_prime_power_heat_block_"
        "composition_guard.json"
    ),
    "critical_frame": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_normalized_prefix_phase_flux_reduction.json"
    ),
}
EXPECTED_IDS = [
    "pplpf_01_chart",
    "pplpf_02_factorization",
    "pplpf_03_currents",
    "pplpf_04_moments",
    "pplpf_05_clock",
    "pplpf_06_reduced_derivative",
    "pplpf_07_shifted_speed",
    "pplpf_08_frozen_current",
    "pplpf_09_external_phase",
    "pplpf_10_heat_profile",
    "pplpf_11_geometric_boundary",
    "pplpf_12_witness_chart",
    "pplpf_13_phase_sweep",
    "pplpf_14_heat_separation",
    "pplpf_15_actual_stability",
    "pplpf_16_C1_rejection",
    "pplpf_17_replacement",
    "pplpf_18_modulus_payoff",
    "pplpf_19_family_split",
    "pplpf_20_rejoin",
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
    value: Fraction, terms: int
) -> tuple[Fraction, Fraction]:
    """Return exact Taylor lower/upper bounds for exp(value), value>=0."""
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
    tail_upper = first_omitted / (1 - ratio)
    return total, total + tail_upper


def log_two_bounds(terms: int = 30) -> tuple[Fraction, Fraction]:
    """Enclose log(2) by the atanh series at 1/3."""
    z = Fraction(1, 3)
    lower = Fraction(0)
    for index in range(terms + 1):
        lower += 2 * z ** (2 * index + 1) / (2 * index + 1)
    tail = (
        2
        * z ** (2 * terms + 3)
        / ((2 * terms + 3) * (1 - z * z))
    )
    return lower, lower + tail


def sqrt_two_bounds(digits: int = 60) -> tuple[Fraction, Fraction]:
    scale = 10**digits
    lower_integer = isqrt(2 * scale * scale)
    return (
        Fraction(lower_integer, scale),
        Fraction(lower_integer + 1, scale),
    )


def validate_flow_algebra(issues: list[str]) -> None:
    x, h, b, c = sp.symbols("x h b c", real=True, nonzero=True)
    h_re, h_im, e_re, e_im = sp.symbols(
        "H_r H_i E_r E_i", real=True
    )
    q_re, q_im = sp.symbols("Q_r Q_i", real=True)
    H = h_re + sp.I * h_im
    E = e_re + sp.I * e_im
    Q = q_re + sp.I * q_im
    tau = -x * h * b
    rho = -c / b
    saddle = c + sp.I * b

    reduced = sp.expand(
        tau * ((sp.I - rho) * H + E / (h * (-b)))
        - (x * E - x * h * saddle * H)
    )
    if sp.simplify(reduced) != 0:
        issues.append("reduced logarithmic derivative identity failed")

    hq = sp.expand_complex(H * sp.conjugate(Q))
    eq = sp.expand_complex(E * sp.conjugate(Q))
    normalized = sp.expand_complex(
        ((sp.I - rho) * H + E / (h * (-b)))
        * sp.conjugate(Q)
    )
    expected_imaginary = (
        sp.re(hq)
        - rho * sp.im(hq)
        + sp.im(eq) / (h * (-b))
    )
    if sp.simplify(sp.im(normalized) - expected_imaginary) != 0:
        issues.append("division-free ray-current identity failed")


def validate_heat_profile(issues: list[str]) -> None:
    k, length = sp.symbols("k M", integer=True, nonnegative=True)
    h, t, u, delta = sp.symbols(
        "h t u delta", positive=True, real=True
    )
    u_zero = u + (length - 1) * h
    saddle_b = sp.Rational(1, 2) - t * delta / 2
    original = (
        -k * h * saddle_b
        - t * k * h * u_zero / 2
        + t * k**2 * h**2 / 4
    )
    log_q = -t * h**2 / 4
    v = 2 * (u - delta) / h
    compact = -k * h / 2 + log_q * k * (
        2 * length - 2 - k + v
    )
    if sp.simplify(original - compact) != 0:
        issues.append("dimensionless heat-profile identity failed")


def witness_coarse_bounds(issues: list[str]) -> dict[str, Fraction]:
    base = 562_538_277
    cutoff = 72_004_899_456
    if base % 2 != 1:
        issues.append("witness base is not dyadic-free")
    if cutoff != 128 * base:
        issues.append("witness chain does not have exactly eight levels")

    left_scale = Fraction(cutoff * cutoff) - Fraction(1, 64)
    _, exp_50_upper = exp_bounds_positive(Fraction(50), 300)
    if not left_scale > exp_50_upper:
        issues.append("witness left boundary is not rigorously above L=50")
    if not cutoff * cutoff < 2**100:
        issues.append("witness does not rigorously satisfy L<100")

    log_lower, log_upper = log_two_bounds()
    if not log_lower > Fraction(2, 3):
        issues.append("log(2) phase-sweep lower bound failed")
    if not log_lower * (2 * cutoff + 1) > 1:
        issues.append("witness cell does not sweep one full base phase")

    sqrt_lower, sqrt_upper = sqrt_two_bounds()
    geometric_lower: list[Fraction] = []
    geometric_upper: list[Fraction] = []
    for index in range(8):
        if index % 2 == 0:
            value = Fraction(1, 2 ** (index // 2))
            geometric_lower.append(value)
            geometric_upper.append(value)
        else:
            factor = Fraction(1, 2 ** ((index - 1) // 2))
            geometric_lower.append(factor / sqrt_upper)
            geometric_upper.append(factor / sqrt_lower)

    geometric_mass_upper = sum(geometric_upper)
    geometric_moment_lower = sum(
        index * geometric_lower[index] for index in range(8)
    )
    geometric_speed_lower = (
        1 + geometric_moment_lower / geometric_mass_upper
    )
    if not geometric_speed_lower > Fraction(14, 5):
        issues.append("geometric aligned-speed lower enclosure failed")

    heated_lower: list[Fraction] = []
    heated_upper: list[Fraction] = []
    for index in range(8):
        if index == 0:
            heated_lower.append(Fraction(1))
            heated_upper.append(Fraction(1))
            continue
        linear = Fraction(index, 2)
        quadratic = Fraction(14 * index - index * index, 16)
        exponent_lower = (
            linear * log_lower
            + quadratic * log_lower * log_lower
        )
        exponent_upper = (
            linear * log_upper
            + quadratic * log_upper * log_upper
        )
        exp_lower, _ = exp_bounds_positive(exponent_lower, 120)
        _, exp_upper = exp_bounds_positive(exponent_upper, 120)
        heated_lower.append(1 / exp_upper)
        heated_upper.append(1 / exp_lower)

    heated_mass_lower = sum(heated_lower)
    heated_moment_upper = sum(
        index * heated_upper[index] for index in range(8)
    )
    heated_speed_upper = 1 + heated_moment_upper / heated_mass_lower
    if not heated_speed_upper < Fraction(11, 5):
        issues.append("heated aligned-speed upper enclosure failed")
    if not (
        geometric_speed_lower - heated_speed_upper > Fraction(3, 5)
    ):
        issues.append("heat/geometric coarse separation failed")

    margin_upper = Fraction(22, 15) - Fraction(7, 5)
    if not margin_upper == Fraction(1, 15):
        issues.append("dyadic margin rational upper bound failed")

    # x=4*pi*left_scale > 12*left_scale because pi>3.
    d_upper = Fraction(2189, 12) / left_scale
    coefficient_relative_error = (
        Fraction(4, cutoff) + 4 * d_upper
    )
    if not coefficient_relative_error < Fraction(1, 10**10):
        issues.append("actual coefficient perturbation bound failed")

    # The ratio perturbation is <3r. Radial and epsilon currents are
    # bounded coarsely using -b>1/2 and log(2)>2/3.
    radial_upper = Fraction(13, 10) / (48 * left_scale)
    epsilon_upper = Fraction(8446, 144) / left_scale**2
    correction_current_upper = 6 * epsilon_upper
    actual_speed_error = (
        3 * coefficient_relative_error
        + radial_upper
        + correction_current_upper
    )
    if not actual_speed_error < Fraction(1, 10**8):
        issues.append("actual logarithmic-flow perturbation bound failed")

    if not (
        geometric_speed_lower
        - heated_speed_upper
        - actual_speed_error
        > Fraction(1, 2)
    ):
        issues.append("actual C1 counterexample separation failed")

    return {
        "geometric_speed_lower": geometric_speed_lower,
        "heated_speed_upper": heated_speed_upper,
        "actual_speed_error": actual_speed_error,
        "margin_upper": margin_upper,
    }


def validate_numeric_audit(stored: dict, issues: list[str]) -> None:
    audit = stored.get("exact", {}).get("numeric_audit", {})
    if audit.get("prime_free_base") != 562_538_277:
        issues.append("stored witness base drifted")
    if audit.get("cutoff") != 72_004_899_456:
        issues.append("stored witness cutoff drifted")
    if audit.get("arithmetic") != (
        "100-digit evaluation; the checker independently proves "
        "the displayed coarse inequalities with rational enclosures"
    ):
        issues.append("stored arithmetic description drifted")

    mp.mp.dps = 80
    geometric = mp.mpf(audit.get("geometric_aligned_speed", "nan"))
    heated = mp.mpf(audit.get("heated_aligned_speed", "nan"))
    margin = mp.mpf(audit.get("geometric_margin_mu_8", "nan"))
    if not geometric > mp.mpf(14) / 5:
        issues.append("stored geometric diagnostic violates coarse bound")
    if not heated < mp.mpf(11) / 5:
        issues.append("stored heat diagnostic violates coarse bound")
    if not geometric - heated > mp.mpf(3) / 5:
        issues.append("stored diagnostic separation violates coarse bound")
    if not 0 < margin < mp.mpf(1) / 15:
        issues.append("stored dyadic diagnostic margin drifted")


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

    chain = sources["complete_chain"].get("chain", {})
    if "a_(m,k)*w_(m,p)^k" not in chain.get(
        "exact_factorization", ""
    ):
        issues.append("source chain factorization drifted")
    ray = sources["ray_flow"].get("exact", {})
    if "x*partial_x" not in ray.get("logarithmic_derivative", {}).get(
        "operator", ""
    ):
        issues.append("source logarithmic operator drifted")
    geometric = sources["geometric_phase"].get("exact", {})
    if "22/15-sqrt(2)" not in geometric.get("thresholds", {}).get(
        "dyadic", ""
    ):
        issues.append("source geometric margin drifted")
    heat = sources["prime_power_heat"].get("exact", {})
    if "8756/x" not in heat.get("correction_phase", ""):
        issues.append("source correction-phase bound drifted")


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
        "rows": 20,
        "exact_logarithmic_phase_identities": 7,
        "dimensionless_heat_profiles": 1,
        "complete_dyadic_counterexamples": 1,
        "rejected_absolute_C1_targets": 1,
        "replacement_starlikeness_targets": 1,
        "proved_joined_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    exact = stored.get("exact", {})
    if "J_ray=|Q|^2+Re(H*conj(Q))" not in exact.get(
        "phase_flow", {}
    ).get("ray_current", ""):
        issues.append("stored ray current drifted")
    if "q^[k*(2M-2-k+v)]" not in exact.get(
        "heat_profile", {}
    ).get("coefficient", ""):
        issues.append("stored heat profile drifted")
    if ">1/2" not in exact.get("dyadic_counterexample", {}).get(
        "rejection", ""
    ):
        issues.append("stored C1 rejection drifted")
    if "J_heat=|Q_heat|^2" not in exact.get(
        "replacement_target", {}
    ).get("current", ""):
        issues.append("stored replacement target drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "rejects the previous absolute geometric C1",
        "does not prove that target uniformly",
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
            "# Prime-Power Logarithmic Phase-Flow Gate",
            "## Exact Chain Flow",
            "## Complete Counterexample",
            "## Replacement Target",
            "0 joined Abel gaps",
            "0 successor winding bounds",
        ):
            if required not in note:
                issues.append(f"note marker missing: {required}")

    validate_flow_algebra(issues)
    validate_heat_profile(issues)
    witness_coarse_bounds(issues)
    validate_numeric_audit(stored, issues)
    validate_sources(stored, sources, issues)
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"ERROR: {issue}")
    if issues:
        return 1
    print(
        "validated prime-power logarithmic phase-flow gate: "
        "20 rows, 7 exact logarithmic-phase identities, "
        "1 dimensionless heat profile, 1 complete dyadic counterexample, "
        "1 rejected absolute C1 target, 1 replacement starlikeness target, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
