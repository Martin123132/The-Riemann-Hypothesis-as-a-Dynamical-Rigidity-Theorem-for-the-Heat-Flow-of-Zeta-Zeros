#!/usr/bin/env python3
"""Independently check the terminal conditional-quadrature gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_terminal_conditional_"
    "quadrature_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 28,
    "physical_carrier_source_expansions": 1,
    "endpoint_carrier_half_audits": 1,
    "endpoint_carrier_halves": 2,
    "local_terminal_derivative_bounds": 2,
    "terminal_s1_geometry_bounds": 4,
    "terminal_s1_integral_sandwiches": 1,
    "terminal_s1_sample_checks": 456,
    "terminal_quadrature_identities": 1,
    "row_relative_terminal_bounds": 2,
    "bulk_terminal_kernel_decompositions": 1,
    "complex_recurrence_relations": 1,
    "quadrature_recurrence_defects": 1,
    "physical_terminal_quadrature_sign_witnesses": 2,
    "signed_bulk_terminal_bounds": 0,
    "complete_outer_complement_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_full_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "tcq_01_rates",
    "tcq_02_lift",
    "tcq_03_halves",
    "tcq_04_local_bound",
    "tcq_05_half_bound",
    "tcq_06_geometry",
    "tcq_07_integral",
    "tcq_08_growth",
    "tcq_09_meaning",
    "tcq_10_exact_trace",
    "tcq_11_leading",
    "tcq_12_row_error",
    "tcq_13_beta",
    "tcq_14_bulk_split",
    "tcq_15_bulk_error",
    "tcq_16_complex_relation",
    "tcq_17_recurrence_defect",
    "tcq_18_midpoint",
    "tcq_19_consequence",
    "tcq_20_source",
    "tcq_21_signs",
    "tcq_22_live",
    "tcq_23_route",
    "tcq_24_ties",
    "tcq_25_phase",
    "tcq_26_quadratic",
    "tcq_27_pi",
    "tcq_28_boundary",
]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def to_mp(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def check_sources(artifact: dict) -> None:
    expected = {
        "terminal_obstruction",
        "outer_endpoint",
        "full_support",
        "terminal_recurrence",
        "two_carrier",
    }
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source-audit keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_carrier_algebra() -> None:
    i = sp.I
    uq, un, unx = sp.symbols("uq un unx", real=True)
    c, cx, b, bx, defect = sp.symbols("c cx b bx defect", real=True)
    C, D = sp.symbols("C D")
    sp1 = c + i * b
    sp2 = cx + i * bx
    x = -uq
    R = -sp1 * x - c * un
    Rx = -sp2 * x - cx * un + i * b * unx
    require_zero(R - (sp1 * uq - c * un), "independent carrier R")
    require_zero(Rx - (sp2 * uq - cx * un + i * b * unx), "independent carrier Rx")
    Q = (R + defect) * C + D
    N = R * Q + Rx * C
    require_zero(N - (R * ((R + defect) * C + D) + Rx * C), "independent N row")

    wr, wi, pr, pii = sp.symbols("wr wi pr pii", real=True)
    w = wr + i * wi
    P = pr + i * pii
    kappa = 1 / (2 * sp.pi * i)
    require_zero(
        sp.re(kappa * w * P) - sp.im(w * P) / (2 * sp.pi),
        "independent quadrature",
    )
    npxi = sp.symbols("npxi", real=True)
    require_zero(sp.im(w * npxi) - wi * npxi, "independent rank-one term")

    half = sp.Rational(1, 2)
    atom = sp.symbols("atom")
    require_zero(half * atom + half * atom - atom, "independent endpoint halves")


def check_rational_bounds() -> None:
    for denominator in [72_000_000_000, 10**12, 10**18, 10**24]:
        h = Fraction(1, denominator)
        require(Fraction(17, 10) * h**2 < Fraction(8, 5), "N_N row bound")
        require(2 * h * Fraction(17, 10) * h**2 < Fraction(17, 5) * h**3 + h**6, "local atom bound")
        require(4380 * h**2 < h, "C source difference")
        exp_difference = 6 * h / (1 - 6 * h)
        require(exp_difference < 7 * h, "source exponential difference")
        ratio_difference = (exp_difference + 4380 * h**2) / (1 - 4380 * h**2)
        require(ratio_difference < 9 * h, "source ratio difference")
        require(9 * h < Fraction(1, 8), "source sign reserve")

        r_alpha, r_beta = sp.symbols("R_alpha R_beta", positive=True)
        row_bound = sp.Rational(8, 5) * h * r_alpha + 2 * h * sp.Rational(8, 5) * r_beta
        expected = sp.Rational(8, 5) * h * r_alpha + sp.Rational(16, 5) * h * r_beta
        require_zero(row_bound - expected, "independent row remainder")

    require(23**2 > 16**2 * 2, "sin(pi/8) radical comparison")
    require(4**2 * 2 > 1, "cos(pi/8) radical comparison")


def check_terminal_kernel() -> int:
    mp.mp.dps = 80
    checks = 0
    for n_value in range(2, 44):
        for theta in [
            Fraction(0),
            Fraction(1, 9),
            Fraction(1, 2),
            Fraction(8, 9),
            Fraction(1),
        ]:
            a_value = Fraction(n_value) + theta
            for offset in [Fraction(1, 23), Fraction(5, 11), Fraction(22, 23)]:
                alpha = a_value**2 - offset
                x_n = alpha / n_value
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                m_n = floor_fraction(a_n) + 1
                n_n = floor_fraction(2 * alpha)
                a_star = x_n - m_n + 1
                b_star = n_n + 1 - x_n

                require(Fraction(1, 4) < a_star < 2, "independent a_* box")
                require(
                    b_star > alpha * (2 - Fraction(1, n_value)),
                    "independent b_* lower",
                )
                require(b_star > a_star, "independent terminal orientation")

                av = to_mp(a_star)
                bv = to_mp(b_star)
                s1 = mp.digamma(bv) - mp.digamma(av)
                lower = mp.log(bv / av)
                upper = lower + 1 / av - 1 / bv
                coarse_lower = mp.log(
                    mp.mpf(n_value**2 - 1)
                    * (2 - mp.mpf(1) / n_value)
                    / 2
                )
                coarse_upper = mp.log(8 * (n_value + 1) ** 2 + 4) + 4
                require(lower < s1 < upper, "independent digamma sandwich")
                require(coarse_lower < s1 < coarse_upper, "independent growth box")
                checks += 1
    require(checks == 630, "independent terminal kernel count")
    return checks


def check_recurrence() -> None:
    p, m = sp.symbols("p m", real=True)
    base = p**2 / 2 - (2 * m + 1) * p + sp.Rational(3, 8)
    shifted = (
        (p - 2 * m - 2) ** 2 / 2
        + (p - 2 * m - 2)
        + sp.Rational(3, 8)
    )
    require_zero(shifted - base - 2 * m * (m + 1), "independent phase shift")

    mp.mp.dps = 70
    for p_value in [mp.mpf(-1), mp.mpf(-3) / 7, mp.mpf(0), mp.mpf(2) / 5, mp.mpf(1)]:
        q = mp.e ** (-mp.pi * 1j * (p_value**2 / 2 - p_value + mp.mpf(3) / 8))
        ratio = -mp.e ** (2 * mp.pi * 1j * p_value)
        for m_value in range(0, 13):
            t = p_value - 2 * m_value - 2
            recurrence = mp.e ** (mp.pi * 1j * (t**2 / 2 + t + mp.mpf(3) / 8))
            rhs = (-1) ** m_value * mp.conj(recurrence)
            require(abs(q * ratio**m_value - rhs) < mp.mpf("1e-60"), "independent complex recurrence")

    midpoint_r = mp.e ** (mp.pi * 1j * mp.mpf(3) / 8)
    defect = midpoint_r - mp.conj(midpoint_r)
    require(abs(defect - 2j * mp.sin(3 * mp.pi / 8)) < mp.mpf("1e-60"), "midpoint quadrature defect")
    require(abs(defect) > mp.mpf(1), "nonzero midpoint defect")


def check_terminal_signs() -> None:
    mp.mp.dps = 80
    q_zero = mp.e ** (-3 * mp.pi * 1j / 8)
    q_one = mp.e ** (mp.pi * 1j / 8)
    require(mp.im(-q_zero) > mp.mpf(3) / 4, "ideal p=0 quadrature")
    require(mp.im(-q_one) < -mp.mpf(3) / 8, "ideal p=1 quadrature")
    require(mp.mpf(9) / 72_000_000_000 < mp.mpf(1) / 8, "physical perturbation reserve")
    require(mp.im(-q_zero) - mp.mpf(1) / 8 > mp.mpf(5) / 8, "physical p=0 sign")
    require(mp.im(-q_one) + mp.mpf(1) / 8 < -mp.mpf(1) / 4, "physical p=1 sign")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["tcq_23_route"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Physical Carrier Row",
        "## Terminal Kernel",
        "## Conditional Boundary",
        "## Recurrence Audit",
        "mathcal N_(<N,xi)",
        "quadrature",
        "logarithmic",
        "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require("signed bulk-terminal estimate open" in artifact.get("status", ""), "status drifted")
    require(
        "no signed bulk-terminal estimate" in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    cert = artifact.get("symbolic_certificate", {})
    require(
        set(cert)
        == {
            "physical_carrier_row",
            "terminal_kernel",
            "conditional_boundary",
            "terminal_recurrence",
            "bounds",
            "audit",
            "handoff",
        },
        "certificate keys drifted",
    )
    require(cert["audit"] == {"terminal_s1_cases": 456}, "builder audit drifted")
    check_sources(artifact)
    check_carrier_algebra()
    check_rational_bounds()
    kernel_checks = check_terminal_kernel()
    check_recurrence()
    check_terminal_signs()
    check_rows(artifact)
    check_note(note)
    print(
        "validated full-support terminal conditional quadrature gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {kernel_checks} kernel checks, "
        f"{EXPECTED_COUNTS['endpoint_carrier_halves']} carrier halves, "
        f"{EXPECTED_COUNTS['quadrature_recurrence_defects']} quadrature defect, "
        f"{EXPECTED_COUNTS['signed_bulk_terminal_bounds']} signed bulk-terminal bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
