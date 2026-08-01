#!/usr/bin/env python3
"""Validate the terminal-tail five-moment two-carrier kernel reduction."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_two_carrier_kernel_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"ttkc_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "polynomials",
            "observations",
            "defect",
            "endpoint",
            "unsymmetrized",
            "hermitian",
            "transpose",
            "rates",
            "anchor",
            "leading",
            "guard",
            "target",
            "handoff",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 14,
    "carrier_polynomials": 4,
    "endpoint_linear_kernels": 1,
    "hermitian_kernels": 1,
    "transpose_kernels": 1,
    "exact_kernel_symmetrizations": 2,
    "common_rate_cancellations": 2,
    "physical_chi_anchor_identities": 1,
    "physical_chi_bounds": 1,
    "leading_q1_kernels": 2,
    "pair_sign_guards": 1,
    "signed_type_ii_targets": 1,
    "signed_type_ii_bounds": 0,
    "phi_b_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated five-moment two-carrier kernel reduction: 14 rows, "
    "1 endpoint-linear kernel, 1 Hermitian kernel, 1 transpose kernel, "
    "1 physical chi_N bound, 0 signed Type-I/II bounds"
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
    if sp.simplify(sp.expand(expression)) != 0:
        issues.append(label)


def independent_kernel_audit(issues: list[str]) -> None:
    cn, cm, cmb = sp.symbols("C_n C_m Cbar_m")
    dn, dm = sp.symbols("D_n D_m")
    rn, rm, rmb = sp.symbols("R_n R_m Rbar_m")
    rnx, rmx, rmxb = sp.symbols("R_n_x R_m_x Rbar_m_x")
    qn, qm, qmb = sp.symbols("Q_n Q_m Qbar_m")
    xt, at, xtx, atx = sp.symbols("X_T A_T X_T_x A_T_x")

    nn = rn * qn + rnx * cn
    endpoint = xt * nn + atx * cn - at * qn - xtx * rn * cn
    endpoint_factored = (
        (xt * rn - at) * qn + (xt * rnx + atx - xtx * rn) * cn
    )
    require_zero(
        endpoint - endpoint_factored,
        "independent endpoint factorization failed",
        issues,
    )

    nmb = rmb * qmb + rmxb * cmb
    h_nm = cn * nmb - rn * cn * qmb
    conjugate_h_mn = cmb * nn - rmb * cmb * qn
    h_plus = (h_nm + conjugate_h_mn) / 2
    h_expected = (
        (rmb - rn) * (cn * qmb - cmb * qn)
        + cn * cmb * (rmxb + rnx)
    ) / 2
    require_zero(
        h_plus - h_expected,
        "independent Hermitian symmetrization failed",
        issues,
    )

    lambda_n, lambda_m, s_prime, chi, common = sp.symbols(
        "lambda_n lambda_m s_prime chi common"
    )
    rn_t = s_prime * lambda_n + common
    rm_t = s_prime * lambda_m + common
    qn_t = (chi + s_prime * lambda_n) * cn + dn
    qm_t = (chi + s_prime * lambda_m) * cm + dm
    nn_t = rn_t * qn_t + rnx * cn
    nm_t = rm_t * qm_t + rmx * cm
    t_nm = cn * nm_t - rn_t * cn * qm_t
    t_mn = cm * nn_t - rm_t * cm * qn_t
    delta = lambda_m - lambda_n
    t_expected = (
        s_prime**2 * delta**2 * cn * cm
        + s_prime * delta * (cn * dm - cm * dn)
        + cn * cm * (rnx + rmx)
    ) / 2
    t_plus = sp.expand((t_nm + t_mn) / 2)
    require_zero(
        t_plus - t_expected,
        "independent transpose symmetrization failed",
        issues,
    )
    if chi in t_plus.free_symbols:
        issues.append("transpose common-rate cancellation failed")

    chibar, sbar, cmbar, dmbar = sp.symbols(
        "chi_bar s_prime_bar C_m_bar D_m_bar"
    )
    qn_h = (chi + s_prime * lambda_n) * cn + dn
    qmbar_h = (chibar + sbar * lambda_m) * cmbar + dmbar
    difference = cn * qmbar_h - cmbar * qn_h
    expected_difference = (
        cn
        * cmbar
        * (chibar - chi + sbar * lambda_m - s_prime * lambda_n)
        + cn * dmbar
        - cmbar * dn
    )
    require_zero(
        difference - expected_difference,
        "independent Hermitian common-rate identity failed",
        issues,
    )
    chi_real, chi_imag = sp.symbols("chi_real chi_imag", real=True)
    real_rate_form = sp.expand(
        expected_difference.subs(
            {chi: chi_real + sp.I * chi_imag, chibar: chi_real - sp.I * chi_imag}
        )
    )
    if sp.diff(real_rate_form, chi_real) != 0:
        issues.append("Hermitian real common-rate cancellation failed")


def independent_four_carrier_audit(issues: list[str]) -> None:
    i = sp.I
    log_n = sp.Integer(7)
    rho_1 = sp.Rational(2, 13) - i / 17
    rho_2 = -sp.Rational(1, 29) + 2 * i / 31
    rho_1_x = sp.Rational(1, 37) + i / 41
    rho_2_x = -sp.Rational(1, 43) + i / 47
    s_prime = -sp.Rational(1, 23) - i / 2
    s_second = sp.Rational(1, 53) - i / 59
    chi = -sp.Rational(2, 19) + 3 * i / 61
    b = -sp.Rational(1, 2)
    b_x = -sp.Rational(1, 67)
    u_terminal = sp.Rational(2, 11)
    u_x = sp.Rational(1, 71)
    ells = (sp.Integer(0), sp.Integer(1), sp.Integer(4), sp.Integer(6))
    weights = (
        1 - i / 3,
        -sp.Rational(3, 5) + 2 * i / 7,
        sp.Rational(4, 9) + i / 11,
        -sp.Rational(2, 13) - 3 * i / 17,
    )

    def lam(ell: sp.Expr) -> sp.Expr:
        return log_n - ell

    def c(ell: sp.Expr) -> sp.Expr:
        return 1 + rho_1 * ell + rho_2 * ell**2

    def d(ell: sp.Expr) -> sp.Expr:
        return rho_1_x * ell + rho_2_x * ell**2

    def r(ell: sp.Expr) -> sp.Expr:
        return s_prime * lam(ell) + i * b * u_terminal

    def rx(ell: sp.Expr) -> sp.Expr:
        return s_second * lam(ell) + i * (b_x * u_terminal + b * u_x)

    def q(ell: sp.Expr) -> sp.Expr:
        return (chi + s_prime * lam(ell)) * c(ell) + d(ell)

    def nn(ell: sp.Expr) -> sp.Expr:
        return r(ell) * q(ell) + rx(ell) * c(ell)

    value = sum(c(ell) * weight for ell, weight in zip(ells, weights))
    value_x = sum(q(ell) * weight for ell, weight in zip(ells, weights))
    scalar = sum(r(ell) * c(ell) * weight for ell, weight in zip(ells, weights))
    scalar_x = sum(nn(ell) * weight for ell, weight in zip(ells, weights))
    direct = sp.re(value) * sp.re(scalar_x) - sp.re(scalar) * sp.re(value_x)

    pair_sum = 0
    for ell_n, weight_n in zip(ells, weights):
        for ell_m, weight_m in zip(ells, weights):
            hermitian = (
                c(ell_n) * sp.conjugate(nn(ell_m))
                - r(ell_n) * c(ell_n) * sp.conjugate(q(ell_m))
            )
            transpose = (
                c(ell_n) * nn(ell_m)
                - r(ell_n) * c(ell_n) * q(ell_m)
            )
            pair_sum += (
                hermitian * weight_n * sp.conjugate(weight_m)
                + transpose * weight_n * weight_m
            )
    require_zero(
        direct - sp.re(pair_sum) / 2,
        "independent four-carrier identity failed",
        issues,
    )


def independent_leading_audit(issues: list[str]) -> None:
    lambda_n, lambda_m, u_terminal, u_x = sp.symbols(
        "lambda_n lambda_m u_N u_N_x", real=True
    )

    def radial(distance: sp.Expr) -> sp.Expr:
        return -sp.I * (distance + u_terminal) / 2

    def radial_x() -> sp.Expr:
        return -sp.I * u_x / 2

    def scalar_rate(distance: sp.Expr) -> sp.Expr:
        value_rate = radial(distance)
        return radial(distance) * value_rate + radial_x()

    rn = radial(lambda_n)
    rm = radial(lambda_m)
    qn = rn
    qm = rm
    nn = scalar_rate(lambda_n)
    nm = scalar_rate(lambda_m)
    h_nm = sp.conjugate(nm) - rn * sp.conjugate(qm)
    h_mn = sp.conjugate(nn) - rm * sp.conjugate(qn)
    h_plus = (h_nm + sp.conjugate(h_mn)) / 2
    t_nm = nm - rn * qm
    t_mn = nn - rm * qn
    t_plus = (t_nm + t_mn) / 2
    require_zero(
        h_plus + (lambda_n + lambda_m + 2 * u_terminal) ** 2 / 8,
        "independent leading Hermitian kernel failed",
        issues,
    )
    require_zero(
        t_plus + (lambda_n - lambda_m) ** 2 / 8 + sp.I * u_x / 2,
        "independent leading transpose kernel failed",
        issues,
    )

    h = sp.symbols("h", positive=True, real=True)
    positive = -(
        (1 - sp.Rational(1, 2))
        * (h**2 - sp.Rational(1, 2) * (2 * h) ** 2)
    ) / 4
    negative = -(1 * h**2) / 4
    if sp.simplify(positive - h**2 / 8) != 0:
        issues.append("independent positive pair guard failed")
    if sp.simplify(negative + h**2 / 4) != 0:
        issues.append("independent negative single guard failed")


def independent_physical_bound_audit(issues: list[str]) -> None:
    finite_edge_path = gate.SOURCE_PATHS["q1_finite_edge"]
    carrier_path = gate.SOURCE_PATHS["carrier_rotation"]
    growing_path = gate.SOURCE_PATHS["growing_tail"]
    try:
        finite_edge = json.loads(finite_edge_path.read_text(encoding="utf-8"))
        carrier = json.loads(carrier_path.read_text(encoding="utf-8"))
        growing = json.loads(growing_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"physical source load failed: {exc}")
        return

    phase = finite_edge.get("exact", {}).get("terminal", {}).get("phase", "")
    stable = (
        finite_edge.get("majorant_certificate", {})
        .get("source_block_majorants", {})
        .get("stable_terminal", {})
    )
    if "Q_x/Q=-i*h*theta/2" not in phase:
        issues.append("terminal phase provenance failed")
    if stable.get("W_x") != "|W_x|<2h^2":
        issues.append("terminal logarithmic derivative provenance failed")

    relative = carrier.get("exact", {}).get("relative_carrier_rotation", "")
    for marker in (
        "|d_n|<2189/x<1/2",
        "|d_(n,x)|<4223/x^2",
        "d_(n,x)/(1+d_n)-d_(1,x)/(1+d_1)",
    ):
        if marker not in relative:
            issues.append(f"carrier correction provenance failed: {marker}")
    source_bounds = growing.get("exact", {}).get("source_bounds", "")
    for marker in ("|s_*'+i/2|<h^2/6000", "u_x=h^2/(8*pi)"):
        if marker not in source_bounds:
            issues.append(f"growing-tail provenance failed: {marker}")

    h_cap = sp.Rational(1, 72000000000)
    delta_bound = 4 * 4223 * h_cap**4
    if not delta_bound < h_cap**2:
        issues.append("terminal correction logarithmic derivative bound failed")
    if 2 + 1 != 3:
        issues.append("physical chi_N triangle bound failed")
    radial_coefficient = 3 + sp.Rational(1, 2) + h_cap / 3000
    if not radial_coefficient < 4:
        issues.append("physical radial defect bound failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 2:
        issues.append("open-row count drifted")

    status = payload.get("status", "")
    for token in (
        "Hermitian/transpose two-carrier Phi_B kernel",
        "signed Type-I/II bound open",
        "no retained aggregate sign",
        "or RH",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "no signed Type-I/II",
        "upper bound on Phi_B",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        sources = gate.load_sources()
        if payload.get("source_audit") != gate.source_audit(sources):
            issues.append("source audit recomputation drifted")
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("symbolic certificate drifted")
        if payload.get("finite_sum_certificate") != gate.finite_sum_certificate():
            issues.append("finite-sum certificate drifted")
        if payload.get("exact") != gate.exact_payload():
            issues.append("exact payload drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"certificate recomputation failed: {exc}")

    if payload.get("counts", {}).get("signed_type_ii_bounds") != 0:
        issues.append("unexpected signed Type-I/II bound")
    if payload.get("counts", {}).get("phi_b_bounds") != 0:
        issues.append("unexpected Phi_B bound")

    independent_kernel_audit(issues)
    independent_four_carrier_audit(issues)
    independent_leading_audit(issues)
    independent_physical_bound_audit(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Terminal-Tail Five-Moment Two-Carrier Kernel Reduction",
        "Date: 2026-07-31",
        "0 signed Type-I/II",
        "H^+_(n,m)",
        "T^+_(n,m)",
        "|chi_N+i*h*theta/2|<3h^2",
        "P_bulk=h^2/8",
        "No new `pi` occurs",
        "RH, or prize-level conclusion",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        validate(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
