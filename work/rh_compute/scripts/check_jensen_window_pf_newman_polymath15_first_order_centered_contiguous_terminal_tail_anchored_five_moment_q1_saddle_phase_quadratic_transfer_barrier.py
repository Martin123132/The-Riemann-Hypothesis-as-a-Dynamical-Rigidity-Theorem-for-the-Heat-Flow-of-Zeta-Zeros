#!/usr/bin/env python3
"""Validate the q=1 saddle-phase quadratic transfer barrier."""

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

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_quadratic_transfer_barrier as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"qpb_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "current",
            "collapse",
            "zero",
            "flow",
            "derivative",
            "line",
            "endpoint",
            "terminal",
            "block",
            "mass",
            "majorant",
            "barrier",
            "reject",
            "handoff",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 15,
    "discriminant_collapses": 1,
    "moment_frequency_flow_identities": 5,
    "signed_line_integrals": 1,
    "endpoint_perturbation_identities": 1,
    "terminal_amplitude_lower_bounds": 1,
    "positive_mass_block_bounds": 2,
    "absolute_mass_transfer_barriers": 1,
    "signed_physical_bounds": 0,
    "phi_b_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated q=1 quadratic transfer barrier: 15 rows, "
    "1 discriminant collapse, 1 signed line integral, "
    "1 terminal amplitude lower bound, 1 absolute-mass barrier, "
    "0 signed physical bounds"
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

    saddle = sources["saddle_variance"]
    if "Omega=-Im(s_*)" not in saddle.get(
        "physical_q1_certificate", {}
    ).get("frequency", ""):
        issues.append("independent saddle-frequency provenance failed")
    if saddle.get("counts", {}).get("signed_physical_bounds") != 0:
        issues.append("independent upstream boundary audit failed")
    stable = sources["finite_edge"].get("majorant_certificate", {}).get(
        "source_block_majorants", {}
    ).get("stable_terminal", {})
    if stable.get("W") != "|W|<6h":
        issues.append("independent terminal-log provenance failed")
    identity = sources["two_carrier"].get(
        "finite_sum_certificate", {}
    ).get("physical_chi_identity", "")
    if "C_Nw_N=-Q(p)exp(W_theta)" not in identity:
        issues.append("independent terminal carrier provenance failed")
    correction = sources["mangoldt"].get(
        "symbolic_certificate", {}
    ).get("quadratic_correction", "")
    if "c_n=(1+d_n)/(1+d_1)" not in correction:
        issues.append("independent correction provenance failed")


def independent_discriminant_audit(issues: list[str]) -> None:
    f, g, p, q, r, s, u = sp.symbols("f g p q r s u", real=True)
    m0 = f + sp.I * g
    m1 = p + sp.I * q
    m2 = r + sp.I * s
    modulus = sp.expand(m0 * sp.conjugate(m0))
    a_coord = sp.re(m1 * sp.conjugate(m0)).expand(complex=True)
    c_coord = sp.expand((m2 * m0 - m1**2) * sp.conjugate(m0) ** 2)
    d_hat = sp.expand(
        f**2 * sp.re(c_coord)
        + a_coord**2 * modulus
        - (2 * u * modulus**2 + sp.im(c_coord)) * f * g
    )
    four_p = -f * r - q**2 + 2 * u * f * g
    require_zero(
        d_hat + four_p * modulus**2,
        "independent D_hat collapse failed",
        issues,
    )


def independent_flow_audit(issues: list[str]) -> None:
    f, g, p, q, r, y, u = sp.symbols("f g p q r y u", real=True)
    current = (-f * r - q**2 + 2 * u * f * g) / 4
    derivative = (
        sp.diff(current, f) * q
        + sp.diff(current, g) * (-p)
        + sp.diff(current, q) * (-r)
        + sp.diff(current, r) * y
    )
    expected = (q * r - f * y + 2 * u * (q * g - f * p)) / 4
    require_zero(
        derivative - expected,
        "independent frequency derivative failed",
        issues,
    )

    df, dg, dq, dr = sp.symbols("df dg dq dr", real=True)
    shifted = (
        -(f + df) * (r + dr)
        - (q + dq) ** 2
        + 2 * u * (f + df) * (g + dg)
    ) / 4
    expected_delta = (
        -df * r
        - f * dr
        - df * dr
        - 2 * q * dq
        - dq**2
        + 2 * u * (df * g + f * dg + df * dg)
    )
    require_zero(
        4 * (shifted - current) - expected_delta,
        "independent endpoint perturbation failed",
        issues,
    )


def independent_finite_line_audit(issues: list[str]) -> None:
    mp.mp.dps = 70
    distances = [mp.mpf(1) / 7, mp.mpf(2) / 5, mp.mpf(3) / 4]
    amplitudes = [mp.mpf(7) / 5, mp.mpf(11) / 7, mp.mpf(13) / 9]
    omega = mp.e ** (mp.j * mp.mpf(3) / 11)
    t0 = mp.mpf(19) / 7
    epsilon = mp.mpf(1) / 131
    omega_phase = t0 - epsilon
    u_x = mp.mpf(1) / 97

    def moments(xi: mp.mpf) -> list[mp.mpc]:
        return [
            omega
            * sum(
                amplitude
                * distance**order
                * mp.e ** (-mp.j * xi * distance)
                for amplitude, distance in zip(amplitudes, distances)
            )
            for order in range(4)
        ]

    def current(xi: mp.mpf) -> mp.mpf:
        m0, m1, m2, _ = moments(xi)
        return (
            -mp.re(m0) * mp.re(m2)
            - mp.im(m1) ** 2
            + 2 * u_x * mp.re(m0) * mp.im(m0)
        ) / 4

    def integrand(theta: mp.mpf) -> mp.mpf:
        m0, m1, m2, m3 = moments(t0 - theta * epsilon)
        f, g = mp.re(m0), mp.im(m0)
        p, q = mp.re(m1), mp.im(m1)
        r, y = mp.re(m2), mp.im(m3)
        return q * r - f * y + 2 * u_x * (q * g - f * p)

    lhs = 4 * (current(omega_phase) - current(t0))
    rhs = -epsilon * mp.quad(integrand, [0, 1])
    if abs(lhs - rhs) > mp.mpf("1e-60"):
        issues.append("independent numerical signed line integral failed")

    positive = [
        sum(
            amplitude * distance**order
            for amplitude, distance in zip(amplitudes, distances)
        )
        for order in range(4)
    ]
    majorant = epsilon * (
        positive[0] * positive[3]
        + positive[1] * positive[2]
        + 4 * u_x * positive[0] * positive[1]
    ) / 4
    if abs(current(omega_phase) - current(t0)) > majorant:
        issues.append("independent numerical absolute majorant failed")


def independent_scale_audit(payload: dict, issues: list[str]) -> None:
    h = sp.Rational(1, 72_000_000_000)
    if not 8 * h < sp.Rational(69, 100):
        issues.append("independent A_a exponential bound failed")
    ratio = (
        sp.Integer(2) ** 50
        * sp.Rational(69, 100) ** 3
        / (96 * sp.Rational(22, 7) * 50**2)
    )
    if not ratio > 490_000_000:
        issues.append("independent transfer-gap rational bound failed")
    saved = payload.get("positive_mass_barrier_certificate", {}).get(
        "rational_certificate", {}
    )
    if saved.get("ratio_lower") != str(ratio):
        issues.append("saved transfer-gap fraction drifted")
    if saved.get("majorant_coefficient") != "1/38400":
        issues.append("saved majorant coefficient drifted")


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
    if any(row.get("readiness") not in {"proved", "guard_validated", "open"} for row in rows):
        issues.append("invalid row readiness")
    if payload.get("counts", {}).get("signed_physical_bounds") != 0:
        issues.append("signed physical bound was overpromoted")
    if payload.get("counts", {}).get("phi_b_bounds") != 0:
        issues.append("Phi_B bound was overpromoted")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not prove the true saddle transfer is large",
        "upper-bound Phi_B",
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
        "# Physical q=1 Quadratic Transfer Barrier",
        "D_hat=(f*r+q^2-2u_x*f*g)|M_0|^4",
        "B_abs/(h^2/400)>490000000",
        "does not lower-bound the true signed transfer error",
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
        independent_discriminant_audit(issues)
        independent_flow_audit(issues)
        independent_finite_line_audit(issues)
        independent_scale_audit(payload, issues)
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
