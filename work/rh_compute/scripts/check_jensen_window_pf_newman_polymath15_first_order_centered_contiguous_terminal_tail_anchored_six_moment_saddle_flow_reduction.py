#!/usr/bin/env python3
"""Validate the endpoint-composed six-moment saddle-flow reduction."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_saddle_flow_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"sfr_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "family",
            "derivative",
            "symmetry",
            "transfer",
            "leading",
            "pair",
            "diagonal",
            "degree",
            "six",
            "mangoldt",
            "hyperbola",
            "guard",
            "target",
            "boundary",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 15,
    "full_current_flow_identities": 1,
    "signed_transfer_integrals": 1,
    "flow_kernel_symmetries": 2,
    "leading_pair_flow_identities": 1,
    "correction_free_moments": 6,
    "max_moment_order": 5,
    "new_mangoldt_orders": 1,
    "order_five_prime_edge_guards": 1,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated six-moment saddle-flow reduction: 15 rows, "
    "1 full-current flow, 2 flow symmetries, 1 leading pair identity, "
    "6 moments, 1 new Mangoldt order, 0 signed flow bounds"
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


def independent_source_audit(payload: dict, issues: list[str]) -> None:
    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")
    try:
        sources = {
            key: json.loads(path.read_text(encoding="utf-8"))
            for key, path in gate.SOURCE_PATHS.items()
        }
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"independent source load failed: {exc}")
        return
    if sources["quadratic_barrier"].get("counts", {}).get(
        "signed_physical_bounds"
    ) != 0:
        issues.append("independent barrier boundary audit failed")
    kernel = sources["two_carrier"].get("symbolic_certificate", {})
    if "h^2*Phi_B" not in kernel.get("full_two_carrier_identity", ""):
        issues.append("independent full-kernel provenance failed")
    mangoldt = sources["mangoldt"].get("symbolic_certificate", {})
    if "Lambda(d)+Lambda(m)" not in mangoldt.get("mangoldt_moments", ""):
        issues.append("independent Mangoldt provenance failed")


def independent_leading_pair_audit(issues: list[str]) -> None:
    count = 4
    u = sp.symbols(f"u0:{count}", real=True)
    x = sp.symbols(f"x0:{count}", real=True)
    y = sp.symbols(f"y0:{count}", real=True)
    u_x = sp.symbols("u_x", real=True)
    f = sum(x)
    g = sum(y)
    p = sum(u[n] * x[n] for n in range(count))
    q = sum(u[n] * y[n] for n in range(count))
    r = sum(u[n] ** 2 * x[n] for n in range(count))
    cubic = sum(u[n] ** 3 * y[n] for n in range(count))
    direct = sp.expand(q * r - f * cubic + 2 * u_x * (q * g - f * p))
    pair = sp.Integer(0)
    for n in range(count):
        for m in range(count):
            delta = u[n] - u[m]
            sigma = u[n] + u[m]
            pair += -sp.Rational(1, 4) * (
                delta**2 * sigma * (x[n] * y[m] + y[n] * x[m])
                + delta
                * sigma**2
                * (y[n] * x[m] - x[n] * y[m])
            ) - u_x * sigma * (x[n] * x[m] - y[n] * y[m])
    require_zero(
        direct - pair,
        "independent leading pair-flow identity failed",
        issues,
    )


def independent_full_flow_numeric_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    u = [mp.mpf(1) / 7, mp.mpf(2) / 5, mp.mpf(5) / 6]
    amplitudes = [mp.mpf(7) / 5, mp.mpf(11) / 9, mp.mpf(13) / 8]
    phase = mp.e ** (mp.j * mp.mpf(2) / 9)
    linear = [
        mp.mpc(2, 1) / 7,
        mp.mpc(-3, 2) / 11,
        mp.mpc(5, -1) / 13,
    ]
    raw_h = [
        [mp.mpc((n + 2) * (m + 1), n - m) / 31 for m in range(3)]
        for n in range(3)
    ]
    hermitian = [
        [
            (raw_h[n][m] + mp.conj(raw_h[m][n])) / 2
            for m in range(3)
        ]
        for n in range(3)
    ]
    raw_t = [
        [mp.mpc(n + m + 1, 2 * n - m) / 29 for m in range(3)]
        for n in range(3)
    ]
    transpose = [
        [(raw_t[n][m] + raw_t[m][n]) / 2 for m in range(3)]
        for n in range(3)
    ]
    t0 = mp.mpf(17) / 5
    epsilon = mp.mpf(1) / 137

    def carriers(xi: mp.mpf) -> list[mp.mpc]:
        return [
            phase * amplitudes[n] * mp.e ** (-mp.j * xi * u[n])
            for n in range(3)
        ]

    def psi(xi: mp.mpf) -> mp.mpf:
        z = carriers(xi)
        value = sum(linear[n] * z[n] for n in range(3))
        value += mp.mpf(1) / 2 * sum(
            hermitian[n][m] * z[n] * mp.conj(z[m])
            + transpose[n][m] * z[n] * z[m]
            for n in range(3)
            for m in range(3)
        )
        return mp.re(value)

    def derivative(xi: mp.mpf) -> mp.mpf:
        z = carriers(xi)
        value = sum(-mp.j * u[n] * linear[n] * z[n] for n in range(3))
        value += mp.mpf(1) / 2 * sum(
            -mp.j
            * (u[n] - u[m])
            * hermitian[n][m]
            * z[n]
            * mp.conj(z[m])
            - mp.j
            * (u[n] + u[m])
            * transpose[n][m]
            * z[n]
            * z[m]
            for n in range(3)
            for m in range(3)
        )
        return mp.re(value)

    lhs = psi(t0 - epsilon) - psi(t0)
    rhs = -epsilon * mp.quad(
        lambda theta: derivative(t0 - theta * epsilon), [0, 1]
    )
    if abs(lhs - rhs) > mp.mpf("1e-60"):
        issues.append("independent full-current line integral failed")


def independent_degree_audit(payload: dict, issues: list[str]) -> None:
    ell_n, ell_m, log_a = sp.symbols("ell_n ell_m log_a")
    c = lambda v: 1 + 2 * v + 3 * v**2
    d = lambda v: 5 + 7 * v + 11 * v**2
    r = lambda v: 13 + 17 * v
    q = lambda v: (19 + 23 * v) * c(v) + d(v)
    r_x = lambda v: 29 + 31 * v
    n_poly = lambda v: sp.expand(r(v) * q(v) + r_x(v) * c(v))
    linear = sp.expand(
        37 * n_poly(ell_n)
        + 41 * c(ell_n)
        - 43 * q(ell_n)
        - 47 * r(ell_n) * c(ell_n)
    )
    pair = lambda a, b: sp.expand(
        c(a) * n_poly(b) - r(a) * c(a) * q(b)
    )
    symmetric = sp.expand(
        (pair(ell_n, ell_m) + pair(ell_m, ell_n)) / 2
    )
    degrees = (
        sp.degree((log_a - ell_n) * linear, ell_n),
        sp.degree((ell_m - ell_n) * symmetric, ell_n),
        sp.degree((ell_m - ell_n) * symmetric, ell_m),
        sp.degree((2 * log_a - ell_n - ell_m) * symmetric, ell_n),
        sp.degree((2 * log_a - ell_n - ell_m) * symmetric, ell_m),
    )
    if max(degrees) != 5 or min(degrees) < 5:
        issues.append("independent degree-five closure failed")
    saved = payload.get("degree_certificate", {})
    if saved.get("max_degree") != 5 or saved.get("moment_count") != 6:
        issues.append("saved six-moment degree audit drifted")


def independent_order_five_audit(issues: list[str]) -> None:
    log_d, log_m, t_symbol, s_symbol = sp.symbols("log_d log_m t_symbol s_symbol")
    combined = t_symbol * (log_d + log_m) ** 2 / 4 - s_symbol * (
        log_d + log_m
    )
    separated = (
        t_symbol * log_d**2 / 4
        - s_symbol * log_d
        + t_symbol * log_m**2 / 4
        - s_symbol * log_m
        + t_symbol * log_d * log_m / 2
    )
    require_zero(
        combined - separated,
        "independent auxiliary heat factor failed",
        issues,
    )
    rows = 0
    for integer in range(2, 64):
        factors = sp.factorint(integer)
        symbols = {
            prime: sp.Symbol(f"P_{prime}") for prime in sorted(factors)
        }
        log_n = sum(
            exponent * symbols[prime]
            for prime, exponent in factors.items()
        )

        def mangoldt(divisor: int) -> sp.Expr:
            decomposition = sp.factorint(divisor)
            if len(decomposition) != 1:
                return sp.Integer(0)
            prime = next(iter(decomposition))
            return symbols.get(prime, sp.Symbol(f"P_{prime}"))

        divisor_sum = sum(mangoldt(divisor) for divisor in sp.divisors(integer))
        require_zero(
            log_n**4 * divisor_sum - log_n**5,
            f"independent order-five identity failed at {integer}",
            issues,
        )
        rows += 1
    if rows != 62:
        issues.append("independent order-five row count failed")


def structural_audit(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("kind drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids drifted")
    if len({row.get("id") for row in rows}) != len(rows):
        issues.append("duplicate row ids")
    if payload.get("counts", {}).get("signed_flow_bounds") != 0:
        issues.append("signed flow bound was overpromoted")
    if payload.get("counts", {}).get("phi_b_bounds") != 0:
        issues.append("Phi_B bound was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not prove a signed flow estimate",
        "Phi_B upper bound",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")


def note_audit(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Endpoint-Composed Six-Moment Saddle Flow",
        "Psi'(xi)=Re sum_n[-i*u_n*L_n]z_n",
        "K_diag=-2u_x sum_n u_n Re(z_n^2)",
        "G_5(xi)=",
        "0 signed flow bounds",
        "not a proof of RH",
        f"python work/rh_compute/scripts/check_{STEM}.py",
    ):
        if marker not in text:
            issues.append(f"note marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        independent_source_audit(payload, issues)
        independent_leading_pair_audit(issues)
        independent_full_flow_numeric_audit(issues)
        independent_degree_audit(payload, issues)
        independent_order_five_audit(issues)
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
