#!/usr/bin/env python3
"""Validate the Suzuki Jordan-error kernel reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_jordan_error_kernel_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md"
)

REQUIRED_IDS = {
    f"jek_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "dirichlet_convolution"),
        (2, "cumulative_mass"),
        (3, "residue_main_term"),
        (4, "error_definition"),
        (5, "weight_derivative"),
        (6, "abel_formula"),
        (7, "archimedean_zero"),
        (8, "derivative_moment_zero"),
        (9, "main_term_annihilation"),
        (10, "error_kernel"),
        (11, "power_sum_decomposition"),
        (12, "elementary_error_bound"),
        (13, "absolute_kernel_bound"),
        (14, "smoothing_barrier"),
        (15, "positive_main_term_rejected"),
        (16, "open_cancellation_target"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Suzuki Jordan-Error Kernel Reduction",
    "A_omega=1/((1+omega)*zeta(1+2omega))",
    "W_(omega,k)(u)=-u*g_(omega,k)'(u)",
    "integral_0^1 W_(omega,k)(u)u^omega du=0",
    "entire positive residue main term",
    "exactly zero to",
    "E_omega(x)",
    "sum_(d>x)mu(d)d^(-(1+2omega))",
    "x^(1/2-omega)*(1+log x)",
    "Increasing `k` does not improve the power",
    "do not prove the required",
    "eventual sign, RH, or Lambda<=0",
)


def mobius(n: int) -> int:
    value = sp.mobius(n)
    return int(value)


def coefficient(omega: mp.mpf, n: int) -> mp.mpf:
    total = mp.mpf("0")
    for d in sp.divisors(n):
        total += mobius(int(d)) * mp.power(d, -omega) * mp.power(
            n // int(d), omega
        )
    return total


def independent_error_decomposition_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    for x in (mp.mpf(7), mp.mpf("12.5"), mp.mpf(31)):
        upper = int(mp.floor(x))
        cumulative = mp.fsum(
            coefficient(omega, n) for n in range(1, upper + 1)
        )
        a_omega = 1 / ((1 + omega) * mp.zeta(1 + 2 * omega))
        direct_error = cumulative - a_omega * mp.power(x, 1 + omega)

        local = mp.mpf("0")
        reciprocal_partial = mp.mpf("0")
        for d in range(1, upper + 1):
            mu_d = mobius(d)
            reciprocal_partial += mu_d * mp.power(
                d, -(1 + 2 * omega)
            )
            y = x / d
            power_sum = mp.fsum(
                mp.power(m, omega)
                for m in range(1, int(mp.floor(y)) + 1)
            )
            remainder = power_sum - mp.power(y, 1 + omega) / (
                1 + omega
            )
            local += mu_d * mp.power(d, -omega) * remainder

        tail = 1 / mp.zeta(1 + 2 * omega) - reciprocal_partial
        decomposed = (
            local
            - mp.power(x, 1 + omega) / (1 + omega) * tail
        )
        if abs(direct_error - decomposed) > mp.mpf("1e-55"):
            issues.append(
                f"error decomposition failed at x={x}: "
                f"{abs(direct_error-decomposed)}"
            )
    return issues


def independent_moment_checks() -> list[str]:
    issues: list[str] = []
    s, omega = sp.symbols("s omega")
    polynomial_factor = (
        (s - omega)
        * (s - omega - 1)
        / ((s + omega) * (s + omega - 1))
    )
    if sp.simplify(polynomial_factor.subs(s, 1 + omega)) != 0:
        issues.append("archimedean factor does not vanish at s=1+omega")

    test_omega = sp.Rational(1, 4)
    exponent = test_omega - 1
    q = sp.Rational(1, 2) - test_omega
    for k in range(1, 5):
        leading_g = -sp.Symbol("a", positive=True) * q ** (-(k - 1))
        leading_w = sp.simplify(-(exponent) * leading_g)
        if leading_w.is_negative is not True:
            issues.append(f"W near-zero sign failed at k={k}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != (
        "jensen_window_pf_suzuki_jordan_error_kernel_reduction"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    ids = {item.get("id") for item in rows}
    if len(rows) != 16:
        issues.append(f"expected 16 rows, found {len(rows)}")
    if ids != REQUIRED_IDS:
        issues.append(
            f"row ids differ: missing={sorted(REQUIRED_IDS-ids)}, "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )
    if len(ids) != len(rows):
        issues.append("duplicate row id")

    for item in rows:
        if not item.get("proof_boundary"):
            issues.append(f"{item.get('id')}: missing proof boundary")

    for row_id in (
        "jek_14_smoothing_barrier",
        "jek_15_positive_main_term_rejected",
    ):
        item = next(
            (candidate for candidate in rows if candidate.get("id") == row_id),
            {},
        )
        if item.get("status") != "guard_validated":
            issues.append(f"{row_id}: not guard validated")

    target = next(
        (
            candidate
            for candidate in rows
            if candidate.get("id") == "jek_16_open_cancellation_target"
        ),
        {},
    )
    if target.get("status") != "open_target":
        issues.append("open cancellation target not marked open")

    audit = payload.get("audit", {})
    expected = {
        "row_count": 16,
        "exact_reduction_count": 11,
        "route_obstruction_count": 2,
        "open_arithmetic_gate_count": 1,
        "main_term_survives": False,
        "eventual_sign_proved": False,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, value in expected.items():
        if audit.get(key) != value:
            issues.append(
                f"audit {key}: expected {value!r}, "
                f"found {audit.get(key)!r}"
            )

    for marker in REQUIRED_NOTE:
        if marker not in note:
            issues.append(f"note missing: {marker}")

    issues.extend(independent_error_decomposition_check())
    issues.extend(independent_moment_checks())

    lowered = note.lower()
    for forbidden in (
        "we have proved rh",
        "this proves rh",
        "positive main term dominates",
        "absolute bound proves eventual sign",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        print(
            "failed Suzuki Jordan-error kernel reduction: "
            f"{len(issues)} issues"
        )
        return 1

    print(
        "validated Suzuki Jordan-error kernel reduction: "
        "16 rows, 0 issues, 8 exact identities/bounds, "
        "2 route guards, 1 open cancellation gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
