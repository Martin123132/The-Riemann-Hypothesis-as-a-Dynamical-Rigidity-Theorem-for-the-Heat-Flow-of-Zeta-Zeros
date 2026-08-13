#!/usr/bin/env python3
"""Validate the physical P_lin joined-coefficient envelope gate."""

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
    "physical_plin_joined_coefficient_envelope_gate"
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
    "prefix_flux": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
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
    "correction_envelope": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_correction_coefficient_envelope_gate.json"
    ),
}


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
    audit = artifact.get("source_audit", {})
    require(set(audit) == set(SOURCE_PATHS), "source key set drifted", issues)
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source {key}", issues)
        if not path.is_file():
            continue
        payloads[key] = load_json(path)
        stored = audit.get(key, {})
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
    effective_box = critical.get("exact", {}).get("domain", {}).get(
        "effective_box", ""
    )
    require("0<t<=1/2" in effective_box, "critical t cap drifted", issues)
    require(
        "h<exp(-25)<1/72000000000" in effective_box,
        "critical h cap drifted",
        issues,
    )
    saddle = critical.get("exact", {}).get("terminal", {}).get("alpha", "")
    require(
        "s_*'=-i/2-it alpha'/4" in saddle,
        "s_*' identity drifted",
        issues,
    )
    require(
        "s_*''=-t alpha''/8" in saddle,
        "s_*'' identity drifted",
        issues,
    )
    source_bounds = critical.get("majorant_certificate", {}).get(
        "source_bounds", {}
    )
    require(
        source_bounds.get("alpha_prime") == "|alpha'|<h^2",
        "alpha-prime bound drifted",
        issues,
    )
    require(
        source_bounds.get("alpha_second") == "|alpha''|<h^4",
        "alpha-second bound drifted",
        issues,
    )

    flux = payloads["prefix_flux"].get("exact", {}).get(
        "normalized_flux_derivative", ""
    )
    require(
        "u_(N,x)=1/(8*pi*a^2)=1/(4*T_0)" in flux,
        "u_(N,x) source drifted",
        issues,
    )

    radial = (
        payloads["two_carrier"]
        .get("finite_sum_certificate", {})
        .get("physical_radial_defect_bound", "")
    )
    for token in (
        "u_N=-log(1-h*theta)",
        "|s_*'+i/2|<h^2/6000",
        "0<=u_N-h*theta<h^2",
        "|i*b*u_N-chi_N|<4h^2",
    ):
        require(token in radial, f"physical radial source drifted: {token}", issues)

    plin = payloads["physical_plin"]
    require(
        plin.get("counts", {}).get("exact_centered_complex_coefficients") == 6,
        "physical P_lin coefficient count drifted",
        issues,
    )
    require(
        "p_5=i*rho_2*(X_T+V_p)*(s_*')^2"
        in plin.get("coefficient_certificate", {}).get("leading_coefficient", ""),
        "physical p_5 source drifted",
        issues,
    )

    correction = payloads["correction_envelope"].get(
        "bound_certificate", {}
    )
    require(
        correction.get("centered_value_bounds", {})
        == {
            "c_0": "|c_0-1|<4379h^2<1/3; hence 2/3<|c_0|<4/3",
            "c_1": "|c_1|<h^2/4",
            "c_2": "|c_2|<h^2/16",
        },
        "C source envelope drifted",
        issues,
    )
    require(
        correction.get("centered_derivative_bounds", {})
        == {
            "d_0": "|d_0|<16893h^4",
            "d_1": "|d_1|<5h^4/16",
            "d_2": "|d_2|<h^4/16",
        },
        "D source envelope drifted",
        issues,
    )


def validate_symbolic(artifact: dict, issues: list[str]) -> None:
    y = sp.symbols("y")
    r_0, r_1, t_0, t_1, delta, u_x = sp.symbols(
        "r_0 r_1 t_0 t_1 delta u_x"
    )
    f_v, f_n, f_a, f_q = sp.symbols("f_v f_n f_a f_q", real=True)
    R = r_0 + r_1 * y
    R_x = t_0 + t_1 * y
    U = sp.expand(
        f_n * (R * (R + delta) + R_x)
        + f_q * (R + delta)
        + f_a * R
        + f_v
    )
    V = sp.expand(f_n * R + f_q)
    u_poly = sp.Poly(U, y)
    v_poly = sp.Poly(V, y)
    u = [u_poly.coeff_monomial(y**j) for j in range(3)]
    v = [v_poly.coeff_monomial(y**j) for j in range(2)]
    expected_u = [
        f_n * (r_0 * (r_0 + delta) + t_0)
        + f_q * (r_0 + delta)
        + f_a * r_0
        + f_v,
        f_n * (r_1 * (2 * r_0 + delta) + t_1)
        + (f_q + f_a) * r_1,
        f_n * r_1**2,
    ]
    expected_v = [f_n * r_0 + f_q, f_n * r_1]
    for index, (actual, expected) in enumerate(zip(u, expected_u)):
        require_zero(actual - expected, f"U recurrence {index} failed", issues)
    for index, (actual, expected) in enumerate(zip(v, expected_v)):
        require_zero(actual - expected, f"V recurrence {index} failed", issues)

    ideal = {
        r_0: 0,
        r_1: sp.I / 2,
        t_0: -sp.I * u_x / 2,
        t_1: 0,
        delta: 0,
    }
    expected_ideal_u = [
        f_v - sp.I * u_x * f_n / 2,
        sp.I * (f_q + f_a) / 2,
        -f_n / 4,
    ]
    expected_ideal_v = [f_q, sp.I * f_n / 2]
    for index, expected in enumerate(expected_ideal_u):
        require_zero(
            u[index].subs(ideal) - expected,
            f"ideal U coefficient {index} failed",
            issues,
        )
    for index, expected in enumerate(expected_ideal_v):
        require_zero(
            v[index].subs(ideal) - expected,
            f"ideal V coefficient {index} failed",
            issues,
        )

    alpha = sp.symbols("alpha_V alpha_N alpha_A alpha_Q", real=True)
    beta = sp.symbols("beta_V beta_N beta_A beta_Q", real=True)

    def substitute_row(
        coefficients: list[sp.Expr], row: tuple[sp.Symbol, ...]
    ) -> list[sp.Expr]:
        substitutions = dict(zip((f_v, f_n, f_a, f_q), row))
        return [sp.expand(item.subs(substitutions)) for item in coefficients]

    u_alpha = substitute_row(u, alpha)
    u_beta = substitute_row(u, beta)
    v_alpha = substitute_row(v, alpha)
    v_beta = substitute_row(v, beta)
    gamma = [
        u_alpha[0],
        u_alpha[1] + sp.I * u_beta[0],
        u_alpha[2] + sp.I * u_beta[1],
        sp.I * u_beta[2],
    ]
    eta = [
        v_alpha[0],
        v_alpha[1] + sp.I * v_beta[0],
        sp.I * v_beta[1],
    ]
    expected_gamma = [
        alpha[0] - sp.I * u_x * alpha[1] / 2,
        sp.I * (alpha[3] + alpha[2]) / 2
        + sp.I * beta[0]
        + u_x * beta[1] / 2,
        -alpha[1] / 4 - (beta[3] + beta[2]) / 2,
        -sp.I * beta[1] / 4,
    ]
    expected_eta = [
        alpha[3],
        sp.I * (alpha[1] / 2 + beta[3]),
        -beta[1] / 2,
    ]
    for index, expected in enumerate(expected_gamma):
        require_zero(
            gamma[index].subs(ideal) - expected,
            f"ideal gamma {index} failed",
            issues,
        )
    for index, expected in enumerate(expected_eta):
        require_zero(
            eta[index].subs(ideal) - expected,
            f"ideal eta {index} failed",
            issues,
        )

    c = sp.symbols("c_0:3")
    d = sp.symbols("d_0:3")
    C = sum(c[j] * y**j for j in range(3))
    D = sum(d[j] * y**j for j in range(3))
    P = sp.expand(
        C * sum(gamma[j] * y**j for j in range(4))
        + D * sum(eta[j] * y**j for j in range(3))
    )
    p = [sp.Poly(P, y).coeff_monomial(y**j) for j in range(6)]
    expected_p = [
        c[0] * gamma[0] + d[0] * eta[0],
        c[0] * gamma[1]
        + c[1] * gamma[0]
        + d[0] * eta[1]
        + d[1] * eta[0],
        c[0] * gamma[2]
        + c[1] * gamma[1]
        + c[2] * gamma[0]
        + d[0] * eta[2]
        + d[1] * eta[1]
        + d[2] * eta[0],
        c[0] * gamma[3]
        + c[1] * gamma[2]
        + c[2] * gamma[1]
        + d[1] * eta[2]
        + d[2] * eta[1],
        c[1] * gamma[3] + c[2] * gamma[2] + d[2] * eta[2],
        c[2] * gamma[3],
    ]
    for index, expected in enumerate(expected_p):
        require_zero(p[index] - expected, f"p_{index} convolution failed", issues)

    rho_2, s_prime = sp.symbols("rho_2 s_prime")
    require_zero(
        p[5].subs({c[2]: rho_2, r_1: -s_prime})
        - sp.I * rho_2 * beta[1] * s_prime**2,
        "p_5 fibre factor failed",
        issues,
    )

    symbolic = artifact.get("symbolic_certificate", {})
    require(
        symbolic.get("symbolic_zero_checks")
        == {
            "row_recurrences": 5,
            "ideal_row_coefficients": 5,
            "ideal_gamma_eta": 7,
            "p_convolution": 6,
            "physical_p5": 1,
        },
        "stored symbolic counts drifted",
        issues,
    )
    require(
        symbolic.get("p_5")
        == "p_5=i*rho_2*beta_N*(s_*')^2=i*rho_2*(X_T+V_p)*(s_*')^2",
        "stored p_5 identity drifted",
        issues,
    )


def validate_bounds(artifact: dict, issues: list[str]) -> None:
    h_cap = Fraction(1, 72_000_000_000)
    require(h_cap < Fraction(1, 2), "h cap is not below one half", issues)
    r_1_abs = Fraction(1, 2) + Fraction(1, 6000)
    require(r_1_abs == Fraction(3001, 6000), "r_1 cap failed", issues)

    t_0_raw = Fraction(1, 64) + Fraction(3001, 144000)
    require(t_0_raw < Fraction(1, 16), "t_0 rational bound failed", issues)

    g_0_raw = (
        Fraction(1, 16)
        + Fraction(1, 144000)
        + Fraction(1, 1500)
        + Fraction(1, 36_000_000)
    )
    require(g_0_raw < Fraction(1, 15), "G_0 rational bound failed", issues)
    g_1_raw = r_1_abs * Fraction(6001, 1500) + Fraction(1, 16)
    require(g_1_raw < Fraction(17, 8), "G_1 rational bound failed", issues)
    g_2_raw = Fraction(1, 6000) + Fraction(1, 36_000_000)
    require(g_2_raw < Fraction(1, 3000), "G_2 rational bound failed", issues)
    p_5_raw = Fraction(1, 16) * r_1_abs**2
    require(p_5_raw < Fraction(1, 63), "p_5 rational bound failed", issues)

    diagnostics = artifact.get("bound_certificate", {}).get(
        "rational_diagnostics", {}
    )
    expected = {
        "h_upper": h_cap,
        "r_1_absolute_cap": r_1_abs,
        "t_0_raw_at_h_half": t_0_raw,
        "G_0_raw_at_h_half": g_0_raw,
        "G_1_raw": g_1_raw,
        "G_2_raw": g_2_raw,
        "p_5_raw_over_h2_betaN": p_5_raw,
    }
    for key, value in expected.items():
        try:
            stored = parse_fraction(diagnostics.get(key, "0/0"))
        except (ValueError, ZeroDivisionError):
            issues.append(f"invalid diagnostic fraction {key}")
            continue
        require(stored == value, f"diagnostic drifted at {key}", issues)

    bounds = artifact.get("bound_certificate", {})
    require(
        bounds.get("rate_box", {})
        == {
            "r_0": "|r_0|<h^3/3000",
            "r_1": "|r_1-i/2|<h^2/6000 and |r_1|<3001/6000",
            "t_0": "|t_0|<h^2/16",
            "t_1": "|t_1|<h^4/16",
            "delta": "|delta|<4h^2",
        },
        "stored rate box drifted",
        issues,
    )
    require(
        bounds.get("quadratic_rate_box", {})
        == {
            "G_0": "|G_0+i*u_(N,x)/2|<h^4/15",
            "G_1": "|G_1|<17h^2/8",
            "G_2": "|G_2+1/4|<h^2/3000",
        },
        "stored quadratic rate box drifted",
        issues,
    )
    p_bounds = bounds.get("p_deviation_envelopes", {})
    for key in ("Pi_0", "Pi_1", "Pi_2", "Pi_3", "Pi_4", "Pi_5"):
        require(key in p_bounds, f"missing {key}", issues)
    require(
        p_bounds.get("Pi_5")
        == "p_5=i*rho_2*beta_N*(s_*')^2 and |p_5|<h^2|beta_N|/63",
        "stored top bound drifted",
        issues,
    )


def validate_structure(
    artifact: dict, artifact_path: Path, note_path: Path, issues: list[str]
) -> None:
    require(artifact.get("kind") == STEM, "artifact kind drifted", issues)
    counts = artifact.get("counts", {})
    expected_counts = {
        "rows": 17,
        "sources": 5,
        "carrier_rate_bounds": 5,
        "quadratic_rate_bounds": 3,
        "row_error_envelopes": 5,
        "ideal_joined_coefficients": 7,
        "centered_plin_coefficients_propagated": 6,
        "exact_top_fibre_factors": 1,
        "numerical_observation_row_bounds": 0,
        "grouped_interior_bounds": 0,
        "signed_flow_bounds": 0,
    }
    require(counts == expected_counts, "artifact counts drifted", issues)
    rows = artifact.get("rows", [])
    require(len(rows) == 17, "row count drifted", issues)
    require(
        len({row.get("id") for row in rows}) == 17,
        "row ids are not unique",
        issues,
    )
    require(note_path.is_file(), "missing companion note", issues)
    if note_path.is_file():
        note = note_path.read_text(encoding="utf-8")
        for token in (
            "# Physical P_lin Joined-Coefficient Envelope Gate",
            "|p_5|<h^2|beta_N|/63",
            "No independent numerical bound for an alpha or beta component has been assumed.",
            "The p_j above multiply y=lambda-log(a).",
            "This gate proves no numerical retained-observation bound",
        ):
            require(token in note, f"note token missing: {token}", issues)
    require(artifact_path.is_file(), "artifact path missing", issues)
    boundary = artifact.get("proof_boundary", "")
    for token in (
        "no numerical retained-observation bound",
        "basis-translated tail bound",
        "grouped Morse estimate",
        "prize-level conclusion",
    ):
        require(token in boundary, f"proof boundary missing: {token}", issues)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    if not args.artifact.is_file():
        print(f"missing artifact: {args.artifact}")
        return 1
    artifact = load_json(args.artifact)
    issues: list[str] = []
    validate_sources(artifact, issues)
    validate_symbolic(artifact, issues)
    validate_bounds(artifact, issues)
    validate_structure(artifact, args.artifact, args.note, issues)

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1

    counts = artifact["counts"]
    print(
        "validated physical P_lin joined-coefficient envelope gate: "
        f"{counts['rows']} rows, {counts['carrier_rate_bounds']} rate bounds, "
        f"{counts['quadratic_rate_bounds']} quadratic-rate bounds, "
        f"{counts['centered_plin_coefficients_propagated']} propagated p_j, "
        f"{counts['numerical_observation_row_bounds']} numerical row bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
