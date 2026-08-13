#!/usr/bin/env python3
"""Validate the physical correction-coefficient envelope gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "physical_plin_correction_coefficient_envelope_gate"
)
DEFAULT_ARTIFACT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "critical_ray": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "critical_ray_finite_height_real_edge_gate.json"
    ),
    "absolute_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "real_residual": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_residual_reduction.json"
    ),
    "two_carrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "physical_plin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_amplitude_gate.json"
    ),
}

H_DENOMINATOR = 72_000_000_000
SMALLNESS_DENOMINATOR = 16_892


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/", 1)
    return Fraction(int(numerator), int(denominator))


def require(condition: bool, message: str, issues: list[str]) -> None:
    if not condition:
        issues.append(message)


def require_zero(expression: sp.Expr, message: str, issues: list[str]) -> None:
    if sp.cancel(sp.expand(expression)) != 0:
        issues.append(message)


def validate_sources(artifact: dict, issues: list[str]) -> None:
    source_audit = artifact.get("source_audit", {})
    require(
        set(source_audit) == set(SOURCE_PATHS),
        "source-audit key set drifted",
        issues,
    )
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source {key}", issues)
        if not path.is_file():
            continue
        payloads[key] = load_json(path)
        stored = source_audit.get(key, {})
        require(
            stored.get("path")
            == str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            f"source path drifted for {key}",
            issues,
        )
        require(
            stored.get("sha256") == file_hash(path),
            f"source hash drifted for {key}",
            issues,
        )

    if set(payloads) != set(SOURCE_PATHS):
        return
    critical = payloads["critical_ray"]
    bounds = critical.get("majorant_certificate", {}).get(
        "source_bounds", {}
    )
    require(bounds.get("alpha_prime") == "|alpha'|<h^2", "alpha' source drifted", issues)
    require(bounds.get("alpha_second") == "|alpha''|<h^4", "alpha'' source drifted", issues)
    require(
        bounds.get("correction") == "|d|<h^2/4, |d_x|<h^4/32",
        "critical correction source drifted",
        issues,
    )
    effective_box = critical.get("exact", {}).get("domain", {}).get(
        "effective_box", ""
    )
    require(
        "h<exp(-25)<1/72000000000" in effective_box,
        "critical h cap drifted",
        issues,
    )
    chain = critical.get("majorant_certificate", {}).get(
        "normalized_majorant_chain", {}
    )
    expected_caps = {
        "chi_absolute": Fraction(2),
        "alpha_prime_over_h2": Fraction(1),
        "alpha_second_over_h4": Fraction(1),
        "d_over_h2": Fraction(1, 4),
        "d_x_over_h4": Fraction(1, 32),
    }
    for key, cap in expected_caps.items():
        row = chain.get(key, {})
        try:
            upper = parse_fraction(row.get("upper", "0/0"))
            stored_cap = parse_fraction(row.get("cap", "0/0"))
        except (ValueError, ZeroDivisionError):
            issues.append(f"invalid critical fraction at {key}")
            continue
        require(stored_cap == cap, f"critical cap drifted at {key}", issues)
        require(upper < stored_cap, f"critical inequality failed at {key}", issues)

    anchor = payloads["absolute_anchor"].get("exact", {}).get(
        "first_coefficient", ""
    )
    require(
        "|d_1|<2189/x<1/2" in anchor,
        "n=1 denominator source drifted",
        issues,
    )
    derivative = (
        payloads["real_residual"]
        .get("exact", {})
        .get("main_and_derivative_bounds", {})
        .get("correction_derivative", "")
    )
    require(
        "|d_(n,x)|<4223/x^2" in derivative,
        "all-carrier derivative source drifted",
        issues,
    )
    scale = (
        payloads["two_carrier"]
        .get("finite_sum_certificate", {})
        .get("physical_chi_bound", "")
    )
    require("x>a^2=h^-2" in scale, "x-to-h scale source drifted", issues)
    require("16892h^2<1" in scale, "smallness source drifted", issues)
    counts = payloads["physical_plin"].get("counts", {})
    require(
        counts.get("exact_centered_complex_coefficients") == 6,
        "parent coefficient count drifted",
        issues,
    )
    require(
        counts.get("maximum_polynomial_degree") == 5,
        "parent degree drifted",
        issues,
    )


def validate_symbolic(artifact: dict, issues: list[str]) -> None:
    ell, A, y = sp.symbols("ell A y")
    alpha, alpha_x = sp.symbols("alpha alpha_x")
    B, B_x, a_0, a_0_x = sp.symbols("B B_x a_0 a_0_x", nonzero=True)

    numerator = a_0 - 2 * B * alpha * ell + B * ell**2
    numerator_x = (
        a_0_x
        - 2 * B_x * alpha * ell
        - 2 * B * alpha_x * ell
        + B_x * ell**2
    )
    C = numerator / a_0
    D = numerator_x / a_0 - numerator * a_0_x / a_0**2

    centered_C = sp.Poly(sp.expand(C.subs(ell, A + y)), y)
    centered_D = sp.Poly(sp.expand(D.subs(ell, A + y)), y)
    rho_2 = B / a_0
    rho_2_x = B_x / a_0 - B * a_0_x / a_0**2
    chi = alpha - A

    require_zero(
        centered_C.coeff_monomial(y**2) - rho_2,
        "independent c_2 identity failed",
        issues,
    )
    require_zero(
        centered_C.coeff_monomial(y) + 2 * chi * rho_2,
        "independent c_1 identity failed",
        issues,
    )
    require_zero(
        centered_C.coeff_monomial(1) - numerator.subs(ell, A) / a_0,
        "independent c_0 ratio failed",
        issues,
    )
    require_zero(
        centered_D.coeff_monomial(y**2) - rho_2_x,
        "independent d_2 identity failed",
        issues,
    )
    require_zero(
        centered_D.coeff_monomial(y)
        + 2 * alpha_x * rho_2
        + 2 * chi * rho_2_x,
        "independent d_1 identity failed",
        issues,
    )
    require_zero(
        centered_D.coeff_monomial(1)
        - (
            numerator_x.subs(ell, A) / a_0
            - numerator.subs(ell, A) * a_0_x / a_0**2
        ),
        "independent d_0 ratio failed",
        issues,
    )

    certificate = artifact.get("symbolic_certificate", {})
    require(
        certificate.get("symbolic_differences")
        == {
            "c_0_ratio": "0",
            "c_1_centering": "0",
            "d_0_ratio": "0",
            "d_1_centering": "0",
        },
        "stored symbolic differences drifted",
        issues,
    )
    derivative_text = certificate.get("normalized_coefficients", "")
    for token in (
        "B_x=-i*alpha''(s)t^2/16",
        "rho_(2,x)=B_x/a_0-B*frak_d_(1,x)/a_0^2",
    ):
        require(token in derivative_text, f"missing derivative token {token}", issues)


def validate_bounds(artifact: dict, issues: list[str]) -> None:
    h_0 = Fraction(1, H_DENOMINATOR)
    require(
        SMALLNESS_DENOMINATOR * h_0**2 < 1,
        "16892 h-cap arithmetic failed",
        issues,
    )
    b = Fraction(1, 32)
    b_x = Fraction(1, 64)
    rho_2 = 2 * b
    rho_2_x = (
        2 * b_x
        + 4 * b * Fraction(4_223, SMALLNESS_DENOMINATOR)
    )
    require(rho_2 == Fraction(1, 16), "rho_2 arithmetic failed", issues)
    require(rho_2_x == Fraction(1, 16), "rho_2_x arithmetic failed", issues)

    c_0_raw = 2 * (Fraction(1, 4) + 2_189)
    require(c_0_raw == Fraction(8_757, 2), "c_0 raw constant failed", issues)
    require(c_0_raw < 4_379, "c_0 rounding failed", issues)
    require(3 * 4_379 < SMALLNESS_DENOMINATOR, "c_0 modulus guard failed", issues)
    c_1 = 2 * 2 * rho_2
    require(c_1 == Fraction(1, 4), "c_1 arithmetic failed", issues)

    d_0_raw = (
        2 * Fraction(1, 32)
        + 4
        * (1 + Fraction(1, 4 * SMALLNESS_DENOMINATOR))
        * 4_223
    )
    require(d_0_raw == Fraction(270_277, 16), "d_0 raw constant failed", issues)
    require(d_0_raw < 16_893, "d_0 rounding failed", issues)
    d_1 = 2 * Fraction(1, 2) * rho_2 + 2 * 2 * rho_2_x
    require(d_1 == Fraction(5, 16), "d_1 arithmetic failed", issues)

    stored = artifact.get("bound_certificate", {})
    diagnostics = stored.get("rational_diagnostics", {})
    expected = {
        "h_upper": Fraction(1, H_DENOMINATOR),
        "16892_h_upper_squared": SMALLNESS_DENOMINATOR * h_0**2,
        "rho_2_over_h2": Fraction(1, 16),
        "rho_2_x_over_h4_using_16892h2_lt_1": Fraction(1, 16),
        "c_0_minus_1_raw_over_h2": Fraction(8_757, 2),
        "c_0_minus_1_rounded_over_h2": Fraction(4_379),
        "d_0_raw_over_h4_using_16892h2_lt_1": Fraction(270_277, 16),
        "d_0_rounded_over_h4": Fraction(16_893),
        "d_1_over_h4": Fraction(5, 16),
    }
    for key, value in expected.items():
        try:
            observed = parse_fraction(diagnostics.get(key, "0/0"))
        except (ValueError, ZeroDivisionError):
            issues.append(f"invalid stored diagnostic {key}")
            continue
        require(observed == value, f"stored diagnostic drifted at {key}", issues)

    quotient = stored.get("quotient_bounds", {})
    value_bounds = stored.get("centered_value_bounds", {})
    derivative_bounds = stored.get("centered_derivative_bounds", {})
    require(quotient.get("rho_2") == "|rho_2|<h^2/16", "stored rho_2 bound drifted", issues)
    require(
        "<h^4/16" in quotient.get("rho_2_x", ""),
        "stored rho_2_x bound drifted",
        issues,
    )
    require(
        value_bounds
        == {
            "c_0": "|c_0-1|<4379h^2<1/3; hence 2/3<|c_0|<4/3",
            "c_1": "|c_1|<h^2/4",
            "c_2": "|c_2|<h^2/16",
        },
        "stored C envelopes drifted",
        issues,
    )
    require(
        derivative_bounds
        == {
            "d_0": "|d_0|<16893h^4",
            "d_1": "|d_1|<5h^4/16",
            "d_2": "|d_2|<h^4/16",
        },
        "stored D envelopes drifted",
        issues,
    )


def validate_rows_and_note(
    artifact: dict, note_path: Path, issues: list[str]
) -> None:
    counts = artifact.get("counts", {})
    expected_counts = {
        "rows": 16,
        "sources": 5,
        "quotient_bounds": 2,
        "centered_coefficients": 6,
        "coefficient_envelopes": 6,
        "nonvanishing_constant_channels": 1,
        "physical_plin_bounds": 0,
        "grouped_interior_bounds": 0,
        "signed_flow_bounds": 0,
    }
    require(counts == expected_counts, "count block drifted", issues)
    rows = artifact.get("rows", [])
    expected_ids = [f"pccenv_{index:02d}_{suffix}" for index, suffix in (
        (1, "domain"),
        (2, "source_separation"),
        (3, "exact_correction"),
        (4, "exact_quotient"),
        (5, "denominator"),
        (6, "B_bounds"),
        (7, "rho_bounds"),
        (8, "value_centering"),
        (9, "derivative_centering"),
        (10, "c0_ratio"),
        (11, "value_envelopes"),
        (12, "d0_ratio"),
        (13, "derivative_envelopes"),
        (14, "constant_channel"),
        (15, "handoff"),
        (16, "boundary"),
    )]
    require([row.get("id") for row in rows] == expected_ids, "row ids drifted", issues)
    require(len({row.get("id") for row in rows}) == 16, "duplicate row ids", issues)
    require(
        artifact.get("status")
        == (
            "exact rho_2 quotient and derivative bounds plus six certified "
            "terminal-centered correction-coefficient envelopes"
        ),
        "status drifted",
        issues,
    )
    boundary = artifact.get("proof_boundary", "")
    for token in ("no P_lin norm", "Phi_B bound", "RH"):
        require(token in boundary, f"proof boundary missing {token}", issues)

    require(note_path.is_file(), "missing note", issues)
    if not note_path.is_file():
        return
    note = note_path.read_text(encoding="utf-8")
    required_note = (
        "|rho_2|<h^2/16",
        "|rho_(2,x)|<h^4/32+(4223/8)h^6<h^4/16",
        "|c_0-1|<4379h^2<1/3",
        "|d_0|<16893h^4",
        "p_5=i*rho_2*(X_T+V_p)*(s_*')^2",
        "not a proof",
    )
    for token in required_note:
        require(token in note, f"note missing {token}", issues)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    issues: list[str] = []
    if not args.artifact.is_file():
        print(f"missing artifact: {args.artifact}")
        return 1
    artifact = load_json(args.artifact)
    require(artifact.get("kind") == STEM, "kind drifted", issues)
    validate_sources(artifact, issues)
    validate_symbolic(artifact, issues)
    validate_bounds(artifact, issues)
    validate_rows_and_note(artifact, args.note, issues)

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    counts = artifact["counts"]
    print(
        "validated physical correction-coefficient envelope gate: "
        f"{counts['rows']} rows, {counts['quotient_bounds']} quotient bounds, "
        f"{counts['centered_coefficients']} centered coefficients, "
        f"{counts['coefficient_envelopes']} coefficient envelopes, "
        f"{counts['nonvanishing_constant_channels']} nonvanishing constant channel, "
        f"{counts['physical_plin_bounds']} P_lin bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
