#!/usr/bin/env python3
"""Validate the physical P_lin coefficient and Morse-amplitude gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_amplitude_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"plma_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "source",
            "center",
            "defect",
            "correction",
            "rows",
            "channels",
            "uv",
            "gamma",
            "coefficients",
            "leading",
            "ideal",
            "cubic",
            "contact",
            "morse",
            "c",
            "cprime",
            "saddle",
            "modes",
            "zero",
            "route",
            "target",
            "quadratic",
            "pi",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 24,
    "physical_base_polynomials": 4,
    "physical_correction_channels": 2,
    "exact_centered_complex_coefficients": 6,
    "maximum_polynomial_degree": 5,
    "degree_five_total_value_factors": 1,
    "correction_free_maximum_degree": 3,
    "correction_free_contact_maximum_degree": 2,
    "exact_morse_cprime_formulas": 2,
    "enumerated_physical_modes": 0,
    "coefficient_norm_bounds": 0,
    "grouped_interior_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}
SUMMARY = (
    "validated physical P_lin/Morse-amplitude gate: 24 rows, "
    "2 correction channels, 6 centered coefficients, degree 5, "
    "1 total-value leading factor, 2 c-prime formulas, "
    "0 grouped interior bounds, 0 signed flow bounds"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def require_zero(expression: sp.Expr, label: str, issues: list[str]) -> None:
    if sp.expand(expression) != 0:
        issues.append(label)


def close(
    left: mp.mpf | mp.mpc,
    right: mp.mpf | mp.mpc,
    tolerance: str = "1e-60",
) -> bool:
    return abs(left - right) <= mp.mpf(tolerance) * max(1, abs(left), abs(right))


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    saved = payload.get("source_audit", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        entry = saved.get(key, {})
        relative = str(path.relative_to(REPO_ROOT)).replace("\\", "/")
        if entry.get("path") != relative:
            issues.append(f"source path drifted: {key}")
        if entry.get("sha256") != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    expected = {
        "observation_image": {
            "visible_residual_observations": 8,
            "maximum_polynomial_degree": 5,
        },
        "two_carrier_kernel": {
            "carrier_polynomials": 4,
            "common_rate_cancellations": 2,
        },
        "morse_fresnel_endpoint": {"exact_transition_integrals": 6},
    }
    for source, checks in expected.items():
        path = gate.SOURCE_PATHS[source]
        if not path.is_file():
            continue
        counts = json.loads(path.read_text(encoding="utf-8")).get("counts", {})
        for key, value in checks.items():
            if counts.get(key) != value:
                issues.append(f"source count drifted: {source}.{key}")


def independent_centering_audit(issues: list[str]) -> None:
    x, log_a, u = sp.symbols("x log_a u", real=True)
    c, b, c_x, b_x, u_x = sp.symbols("c b c_x b_x u_x", real=True)
    s_1 = c + sp.I * b
    s_2 = c_x + sp.I * b_x
    lam = log_a + x
    log_n = log_a - u
    original_r = s_1 * (log_n - lam) + sp.I * b * u
    original_rx = s_2 * (log_n - lam) + sp.I * (b_x * u + b * u_x)
    require_zero(original_r - (-s_1 * x - c * u), "independent R centering failed", issues)
    require_zero(
        original_rx - (-s_2 * x - c_x * u + sp.I * b * u_x),
        "independent R_x centering failed",
        issues,
    )


def independent_polynomial_audit(issues: list[str]) -> None:
    x = sp.symbols("x", real=True)
    c = [sp.Rational(7, 5), sp.Rational(-3, 11) + sp.I / 13, sp.Rational(2, 17)]
    d = [sp.Rational(-5, 19), sp.Rational(3, 23), sp.Rational(1, 29) - sp.I / 31]
    r = [sp.Rational(2, 7) + sp.I / 5, -sp.Rational(4, 9) + sp.I / 8]
    rx = [sp.Rational(-3, 10) + sp.I / 12, sp.Rational(5, 14) - sp.I / 15]
    delta = sp.Rational(1, 6) - sp.I / 16
    c_poly = sum(c[index] * x**index for index in range(3))
    d_poly = sum(d[index] * x**index for index in range(3))
    r_poly = r[0] + r[1] * x
    rx_poly = rx[0] + rx[1] * x
    q_poly = (r_poly + delta) * c_poly + d_poly
    a_poly = r_poly * c_poly
    n_poly = sp.expand(r_poly * q_poly + rx_poly * c_poly)

    alpha = [sp.Rational(2, 3), -sp.Rational(3, 5), sp.Rational(5, 7), -sp.Rational(7, 11)]
    beta = [sp.Rational(-11, 13), sp.Rational(13, 17), -sp.Rational(17, 19), sp.Rational(19, 23)]

    def uv(values: list[sp.Expr]) -> tuple[list[sp.Expr], list[sp.Expr]]:
        f_v, f_n, f_a, f_q = values
        g = [
            r[0] * (r[0] + delta) + rx[0],
            r[1] * (2 * r[0] + delta) + rx[1],
            r[1] ** 2,
        ]
        u_values = [
            f_n * g[0] + f_q * (r[0] + delta) + f_a * r[0] + f_v,
            f_n * g[1] + (f_q + f_a) * r[1],
            f_n * g[2],
        ]
        v_values = [f_n * r[0] + f_q, f_n * r[1]]
        return u_values, v_values

    u_a, v_a = uv(alpha)
    u_b, v_b = uv(beta)
    gamma = [u_a[0], u_a[1] + sp.I * u_b[0], u_a[2] + sp.I * u_b[1], sp.I * u_b[2]]
    eta = [v_a[0], v_a[1] + sp.I * v_b[0], sp.I * v_b[1]]
    recurrence = sp.expand(
        c_poly * sum(gamma[index] * x**index for index in range(4))
        + d_poly * sum(eta[index] * x**index for index in range(3))
    )
    direct_alpha = alpha[0] * c_poly + alpha[1] * n_poly + alpha[2] * a_poly + alpha[3] * q_poly
    direct_beta = beta[0] * c_poly + beta[1] * n_poly + beta[2] * a_poly + beta[3] * q_poly
    require_zero(
        recurrence - direct_alpha - sp.I * x * direct_beta,
        "independent C/D coefficient recurrence failed",
        issues,
    )

    p_5 = sp.Poly(recurrence, x).nth(5)
    expected_p_5 = sp.I * c[2] * beta[1] * r[1] ** 2
    require_zero(p_5 - expected_p_5, "independent p_5 factor failed", issues)


def independent_ideal_audit(payload: dict, issues: list[str]) -> None:
    x, u_x = sp.symbols("x u_x", real=True)
    n_p, v_p, q_p, a_p = sp.symbols("N_p V_p Q_p A_p", real=True)
    n_d, q_d = sp.symbols("N_p_xi Q_p_xi", real=True)
    x_t, a_t, x_tx, a_tx = sp.symbols("X_T A_T X_T_x A_T_x", real=True)
    c_poly = sp.Integer(1)
    a_poly = sp.I * x / 2
    q_poly = a_poly
    n_poly = -x**2 / 4 - sp.I * u_x / 2
    direct = sp.expand(
        n_d * c_poly
        + 2 * a_p * n_poly
        - q_d * a_poly
        - q_d * q_poly
        + sp.I
        * x
        * (
            (n_p + a_tx) * c_poly
            + (x_t + v_p) * n_poly
            - (q_p + x_tx) * a_poly
            - (a_t + a_p) * q_poly
        )
    ).subs(q_p, a_p)
    expected = sp.expand(
        n_d
        - sp.I * a_p * u_x
        + (sp.I * (n_p + a_tx - q_d) + (x_t + v_p) * u_x / 2) * x
        + (a_p + a_t + x_tx) * x**2 / 2
        - sp.I * (x_t + v_p) * x**3 / 4
    )
    require_zero(direct - expected, "independent ideal cubic failed", issues)
    contact = sp.expand(expected.subs({v_p: -x_t, a_p: -a_t}))
    if sp.degree(contact, x) != 2:
        issues.append("independent ideal contact degree failed")
    saved = payload.get("ideal_certificate", {})
    if saved.get("ideal_polynomial_expanded") != str(expected):
        issues.append("saved ideal polynomial drifted")
    if saved.get("contact_polynomial_expanded") != str(contact):
        issues.append("saved contact polynomial drifted")


def morse_z(v: mp.mpf) -> mp.mpf:
    return mp.sign(v - 1) * mp.sqrt(2 * (v - 1 - mp.log(v)))


def morse_j(v: mp.mpf) -> mp.mpf:
    z = morse_z(v)
    return v * z / (v - 1)


def independent_morse_audit(issues: list[str]) -> None:
    mp.mp.dps = 90
    alpha = mp.mpf("53.125")
    mode = 11
    sigma = mp.mpf("0.503")
    heat = mp.mpf(1) / 7000
    coefficients = [
        mp.mpc(2, -1) / 3,
        mp.mpc(-3, 2) / 5,
        mp.mpc(5, 1) / 7,
        mp.mpc(-7, -2) / 11,
        mp.mpc(11, 3) / 13,
        mp.mpc(-13, 5) / 17,
    ]
    u_0 = alpha / mode
    lambda_0 = mp.log(u_0)

    def polynomial(lam: mp.mpf) -> mp.mpc:
        return mp.fsum(coefficients[index] * lam**index for index in range(6))

    def polynomial_1(lam: mp.mpf) -> mp.mpc:
        return mp.fsum(index * coefficients[index] * lam ** (index - 1) for index in range(1, 6))

    def polynomial_2(lam: mp.mpf) -> mp.mpc:
        return mp.fsum(index * (index - 1) * coefficients[index] * lam ** (index - 2) for index in range(2, 6))

    def exponential(lam: mp.mpf) -> mp.mpf:
        return mp.exp(heat * lam**2 / 4 - sigma * lam)

    def w(v: mp.mpf) -> mp.mpc:
        lam = mp.log(u_0 * v)
        return exponential(lam) * polynomial(lam) * morse_j(v)

    w_1 = exponential(lambda_0) * polynomial(lambda_0)
    for v in (mp.mpf("0.43"), mp.mpf("0.82"), mp.mpf("1.27"), mp.mpf("2.4")):
        z = morse_z(v)
        jacobian = morse_j(v)
        direct = (
            z * jacobian * mp.diff(w, v) - w(v) + w_1
        ) / (mode * mp.sqrt(alpha) * z**2)
        lam = mp.log(u_0 * v)
        g = heat * lam / 2 - sigma
        j_prime = mp.diff(morse_j, v)
        closed = (
            exponential(lam)
            * (
                z * jacobian**2 * (polynomial_1(lam) + g * polynomial(lam)) / v
                + (z * jacobian * j_prime - jacobian) * polynomial(lam)
            )
            + w_1
        ) / (mode * mp.sqrt(alpha) * z**2)
        if not close(direct, closed, "1e-70"):
            issues.append(f"independent off-saddle Morse formula failed at v={v}")

    g_0 = heat * lambda_0 / 2 - sigma
    saddle = exponential(lambda_0) / (2 * mode * mp.sqrt(alpha)) * (
        polynomial_2(lambda_0)
        + (2 * g_0 + 1) * polynomial_1(lambda_0)
        + (g_0**2 + g_0 + heat / 2 + mp.mpf(1) / 6) * polynomial(lambda_0)
    )
    epsilon = mp.mpf("1e-10")
    near_values = []
    for v in (1 - epsilon, 1 + epsilon):
        z = morse_z(v)
        jacobian = morse_j(v)
        near_values.append(
            (z * jacobian * mp.diff(w, v) - w(v) + w_1)
            / (mode * mp.sqrt(alpha) * z**2)
        )
    symmetric_limit = mp.fsum(near_values) / 2
    if not close(symmetric_limit, saddle, "1e-18"):
        issues.append("independent removable saddle formula failed")


def independent_saddle_series_audit(issues: list[str]) -> None:
    z = sp.symbols("z", real=True)
    a_0, a_1, a_2 = sp.symbols("A0 A1 A2")
    delta_v = z + z**2 / 3
    delta_lambda = sp.expand(delta_v - delta_v**2 / 2)
    jacobian = sp.expand(1 + sp.Rational(2, 3) * delta_v - sp.Rational(5, 36) * delta_v**2)
    amplitude = sp.expand(a_0 + a_1 * delta_lambda + a_2 * delta_lambda**2 / 2)
    product = sp.series(amplitude * jacobian, z, 0, 3).removeO().expand()
    expected_second = a_2 / 2 + a_1 / 2 + a_0 / 12
    require_zero(
        product.coeff(z, 2) - expected_second,
        "independent saddle series coefficient failed",
        issues,
    )


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len(ids) != len(set(ids)):
        issues.append("duplicate row ids")
    readiness = [row.get("readiness") for row in rows]
    if readiness.count("proved") != 23 or readiness.count("open") != 1:
        issues.append("readiness counts drifted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "no coefficient norm",
        "grouped c-prime sum",
        "no",
        "RH",
        "prize-level conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    try:
        rebuilt = gate.build_payload()
    except Exception as exc:  # pragma: no cover
        issues.append(f"builder replay failed: {exc}")
    else:
        if payload != rebuilt:
            issues.append("saved result differs from a fresh builder replay")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Physical P_lin And Morse-Amplitude Gate",
        "p_5=i*rho_2",
        "c'_(P,r)(0)",
        "0 grouped interior bounds",
        "this is not a proof",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_centering_audit(issues)
        independent_polynomial_audit(issues)
        independent_ideal_audit(payload, issues)
        independent_morse_audit(issues)
        independent_saddle_series_audit(issues)
        structural_audit(payload, issues)
    note_audit(issues)
    if issues:
        for issue in issues:
            print(f"FAIL: {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
